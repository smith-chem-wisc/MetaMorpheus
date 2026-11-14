using EngineLayer;
using EngineLayer.Indexing;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Reading a cached peptide index back restores each peptide's proteins, modifications and digestion settings after
    /// deserializing it. That used to run on one thread and was about half of the read. It now runs in parallel. These tests hold
    /// it to the one-thread loop, kept verbatim below: every peptide comes back with the same fields, the same protein and
    /// modification objects and the same masses, and a cache that cannot be restored fails with the same exception.
    /// </summary>
    [TestFixture]
    public static class ReadPeptideIndexParallelTests
    {
        private const BindingFlags PrivateStatic = BindingFlags.NonPublic | BindingFlags.Static;

        private static readonly MethodInfo ReadPeptideIndex = typeof(MetaMorpheusTask).GetMethod("ReadPeptideIndex", PrivateStatic,
            new[] { typeof(string), typeof(List<Protein>), typeof(int) });

        private static readonly MethodInfo WritePeptideIndex = typeof(MetaMorpheusTask).GetMethod("WritePeptideIndex", PrivateStatic,
            new[] { typeof(List<IBioPolymerWithSetMods>), typeof(string) });

        private static readonly MethodInfo WriteIndexEngineParams = typeof(MetaMorpheusTask).GetMethod("WriteIndexEngineParams", PrivateStatic);

        /// <summary>
        /// The small yeast database with decoys, carbamidomethyl fixed and two variable modifications, so many peptides carry
        /// several modifications, fixed and variable together. It is written as a cache folder.
        /// </summary>
        private static (string Folder, List<Protein> Proteins, List<IBioPolymerWithSetMods> Built) WrittenCache()
        {
            var parameters = new CommonParameters(maxThreadsToUsePerFile: 4,
                digestionParams: new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7, maxModsForPeptides: 3));
            var proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta"),
                true, DecoyType.Reverse, false, out _);
            var variable = GlobalVariables.AllModsKnown.Where(m => (m.ModificationType == "Common Variable" && m.IdWithMotif == "Oxidation on M")
                || (m.ModificationType == "Common Biological" && m.IdWithMotif == "Phosphorylation on S")).ToList();
            var fixedMods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Fixed" && m.IdWithMotif == "Carbamidomethyl on C").ToList();
            var indexEngine = new IndexingEngine(proteins, variable, fixedMods, null, null, null, 0, DecoyType.Reverse, parameters,
                new List<(string, CommonParameters)> { ("", parameters) }, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var built = ((IndexingResults)indexEngine.Run()).PeptideIndex;

            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "ReadPeptideIndexParallel-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(folder);
            WriteIndexEngineParams.Invoke(null, new object[] { indexEngine, Path.Combine(folder, MetaMorpheusTask.IndexEngineParamsFileName) });
            WritePeptideIndex.Invoke(null, new object[] { built, Path.Combine(folder, MetaMorpheusTask.PeptideIndexFileName) });
            return (folder, proteins, built);
        }

        private static List<PeptideWithSetModifications> Read(string folder, List<Protein> proteins, int maxThreads)
        {
            try
            {
                return (List<PeptideWithSetModifications>)ReadPeptideIndex.Invoke(null,
                    new object[] { Path.Combine(folder, MetaMorpheusTask.PeptideIndexFileName), proteins, maxThreads });
            }
            catch (TargetInvocationException e)
            {
                System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(e.InnerException).Throw();
                throw;
            }
        }

        [Test]
        public static void TheReflectedMethodsExist()
        {
            Assert.That(ReadPeptideIndex, Is.Not.Null, "ReadPeptideIndex(string, List<Protein>, int) not found; update this test");
            Assert.That(WritePeptideIndex, Is.Not.Null);
            Assert.That(WriteIndexEngineParams, Is.Not.Null);
        }

        /// <summary>
        /// One thread, eight threads and the old loop restore every peptide identically, field by field.
        /// </summary>
        [Test]
        public static void ParallelRestoreGivesEveryPeptideExactlyWhatTheOneThreadLoopGave()
        {
            var (folder, proteins, built) = WrittenCache();
            try
            {
                List<PeptideWithSetModifications> reference = ReferenceRead(Path.Combine(folder, MetaMorpheusTask.PeptideIndexFileName), proteins);
                Assume.That(reference.Count, Is.EqualTo(built.Count));
                Assume.That(reference.Count(p => p.AllModsOneIsNterminus.Count >= 2), Is.GreaterThan(100), "several modifications on many peptides");

                foreach (int threads in new[] { 1, 8, 8, 8 })
                {
                    List<PeptideWithSetModifications> actual = Read(folder, proteins, threads);
                    Assert.That(actual.Count, Is.EqualTo(reference.Count));
                    var failures = new List<string>();
                    for (int i = 0; i < reference.Count && failures.Count < 10; i++)
                    {
                        string difference = Difference(reference[i], actual[i]);
                        if (difference != null)
                        {
                            failures.Add($"{threads} threads, peptide {i} {reference[i].FullSequence}: {difference}");
                        }
                    }
                    Assert.That(failures, Is.Empty, string.Join("; ", failures));
                }
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A cache whose peptides name proteins the database no longer has cannot be restored. The read stops with the exception
        /// the one-thread loop threw: the one for the first such peptide in the index, whichever thread reaches a failure first.
        /// </summary>
        [Test]
        public static void ACacheThatCannotBeRestoredFailsWithTheExceptionTheOneThreadLoopThrew()
        {
            var (folder, proteins, built) = WrittenCache();
            try
            {
                // drop a protein from every few, so failing peptides are spread through the index and across threads' ranges
                var remaining = proteins.Where((_, i) => i % 7 != 3).ToList();
                var referenceFailure = Assert.Catch(() => ReferenceRead(Path.Combine(folder, MetaMorpheusTask.PeptideIndexFileName), remaining));
                Assume.That(referenceFailure, Is.InstanceOf<MzLibException>());

                for (int attempt = 0; attempt < 5; attempt++)
                {
                    var failure = Assert.Catch(() => Read(folder, remaining, 8));
                    Assert.That(failure.GetType(), Is.EqualTo(referenceFailure.GetType()));
                    Assert.That(failure.Message, Is.EqualTo(referenceFailure.Message), "the first failing peptide, as the one-thread loop reported it");
                }
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        private static string Difference(PeptideWithSetModifications x, PeptideWithSetModifications y)
        {
            if (x.FullSequence != y.FullSequence || x.BaseSequence != y.BaseSequence) return "sequence";
            if (!ReferenceEquals(x.Parent, y.Parent)) return "protein";
            if (BitConverter.DoubleToInt64Bits(x.MonoisotopicMass) != BitConverter.DoubleToInt64Bits(y.MonoisotopicMass)) return "mass";
            if (x.OneBasedStartResidue != y.OneBasedStartResidue || x.OneBasedEndResidue != y.OneBasedEndResidue) return "position";
            if (x.MissedCleavages != y.MissedCleavages || x.NumFixedMods != y.NumFixedMods) return "missed cleavages or fixed mods";
            if (x.CleavageSpecificityForFdrCategory != y.CleavageSpecificityForFdrCategory || x.PeptideDescription != y.PeptideDescription) return "specificity";
            if (x.PairedTargetDecoySequence != y.PairedTargetDecoySequence) return "paired sequence";
            if (DigestionParamsText(x.DigestionParams) != DigestionParamsText(y.DigestionParams)) return "digestion parameters";
            if (!x.AllModsOneIsNterminus.Keys.SequenceEqual(y.AllModsOneIsNterminus.Keys)) return "modification positions or their order";
            if (x.AllModsOneIsNterminus.Any(kv => !ReferenceEquals(kv.Value, y.AllModsOneIsNterminus[kv.Key]))) return "modification objects";
            if (!x.Equals(y)) return "Equals";
            return null;
        }

        /// <summary>
        /// Every peptide of one read shares one DigestionParams object, so each is rendered once rather than once per peptide.
        /// </summary>
        private static readonly System.Runtime.CompilerServices.ConditionalWeakTable<object, string> RenderedDigestionParams = new();

        private static string DigestionParamsText(object digestionParams)
            => digestionParams == null ? "none" : RenderedDigestionParams.GetValue(digestionParams, p => Toml.WriteString(p, MetaMorpheusTask.tomlConfig));

        /// <summary>
        /// MetaMorpheusTask.ReadPeptideIndex as it was before the parallel restore, verbatim apart from reading the digestion
        /// parameters directly rather than through its one-line private helper.
        /// </summary>
        private static List<PeptideWithSetModifications> ReferenceRead(string peptideIndexFileName, List<Protein> allKnownProteins)
        {
            var messageTypes = new[] { typeof(List<PeptideWithSetModifications>) };
            var ser = new NetSerializer.Serializer(messageTypes);
            List<PeptideWithSetModifications> peptideIndex;
            using (var file = File.OpenRead(peptideIndexFileName))
            {
                peptideIndex = (List<PeptideWithSetModifications>)ser.Deserialize(file);
            }

            // populate dictionaries of known proteins for deserialization
            Dictionary<string, Protein> proteinDictionary = new Dictionary<string, Protein>();
            foreach (Protein protein in allKnownProteins)
            {
                if (!proteinDictionary.ContainsKey(protein.Accession))
                {
                    proteinDictionary.Add(protein.Accession, protein);
                }
                else if (proteinDictionary[protein.Accession].BaseSequence != protein.BaseSequence)
                {
                    throw new MetaMorpheusException($"The protein database contained multiple proteins with accession {protein.Accession} ! This is not allowed for index-based searches (modern, non-specific, crosslink searches)");
                }
            }

            // get digestion info from file
            var storedDigestParams = Toml.ReadFile<DigestionParams>(Path.Combine(Path.GetDirectoryName(peptideIndexFileName), "DigestionParameters.toml"), MetaMorpheusTask.tomlConfig);

            // get non-serialized information for the peptides (proteins, mod info)
            foreach (var peptide in peptideIndex)
            {
                peptide.SetNonSerializedPeptideInfo(GlobalVariables.AllModsKnownDictionary, proteinDictionary, storedDigestParams);
            }

            return peptideIndex;
        }
    }
}
