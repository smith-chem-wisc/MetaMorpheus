using EngineLayer;
using Nett;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The cached peptide index stores each peptide's full sequence, and nothing else of its modifications: reading the
    /// index back rebuilds them from it. Digested peptides must therefore be written with their full sequence, whether
    /// or not anything has read it yet, and must come back as they went in.
    /// </summary>
    [TestFixture]
    public static class PeptideIndexRoundTripTest
    {
        private const BindingFlags PrivateStatic = BindingFlags.NonPublic | BindingFlags.Static;

        [Test]
        public static void DigestedPeptidesComeBackFromThePeptideIndexAsTheyWentIn()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "PeptideIndexRoundTrip-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(folder);
            var digestionParams = new DigestionParams(maxMissedCleavages: 2, minPeptideLength: 5);
            Toml.WriteFile(digestionParams, Path.Combine(folder, "DigestionParameters.toml"), MetaMorpheusTask.tomlConfig);

            Modification oxidation = Known("Common Variable", "Oxidation on M");
            Modification carbamidomethyl = Known("Common Fixed", "Carbamidomethyl on C");
            Modification phospho = Known("Common Biological", "Phosphorylation on S");
            var proteins = new List<Protein>
            {
                new("MSEQSKCMSRTPEKSMLCKDEMSSAKR", "P1", oneBasedModifications: new Dictionary<int, List<Modification>>
                    { { 2, new List<Modification> { phospho } }, { 16, new List<Modification> { phospho } } }),
                new("MKCLMSGRNNMQSAMCRK", "P2", isDecoy: true, oneBasedModifications: new Dictionary<int, List<Modification>>
                    { { 6, new List<Modification> { phospho } } }),
            };

            // nothing reads the peptides' full sequences before they are written
            List<PeptideWithSetModifications> written = proteins
                .SelectMany(p => p.Digest(digestionParams, new List<Modification> { carbamidomethyl }, new List<Modification> { oxidation }))
                .ToList();
            Assert.That(written.Count, Is.GreaterThan(30));
            Assert.That(written.Any(p => p.AllModsOneIsNterminus.Count > 1));

            string path = Path.Combine(folder, MetaMorpheusTask.PeptideIndexFileName);
            typeof(MetaMorpheusTask).GetMethod("WritePeptideIndex", PrivateStatic, new[] { typeof(List<PeptideWithSetModifications>), typeof(string) })
                .Invoke(null, new object[] { written, path });
            MethodInfo readPeptideIndex = typeof(MetaMorpheusTask)
                .GetMethod("ReadPeptideIndex", PrivateStatic, new[] { typeof(string), typeof(List<Protein>) });
            Assert.That(readPeptideIndex, Is.Not.Null, "ReadPeptideIndex(string, List<Protein>) not found; update this test");
            var read = (List<PeptideWithSetModifications>)readPeptideIndex.Invoke(null, new object[] { path, proteins });

            Assert.That(read.Count, Is.EqualTo(written.Count));
            for (int i = 0; i < written.Count; i++)
            {
                Assert.That(read[i].FullSequence, Is.EqualTo(written[i].FullSequence));
                Assert.That(read[i].BaseSequence, Is.EqualTo(written[i].BaseSequence));
                Assert.That(read[i].Protein, Is.SameAs(written[i].Protein));
                Assert.That((read[i].OneBasedStartResidue, read[i].OneBasedEndResidue), Is.EqualTo((written[i].OneBasedStartResidue, written[i].OneBasedEndResidue)));
                // rebuilt from the full sequence, so in sequence order rather than the order digestion added them
                Assert.That(read[i].AllModsOneIsNterminus.OrderBy(kv => kv.Key).Select(kv => (kv.Key, kv.Value.IdWithMotif)).ToList(),
                    Is.EqualTo(written[i].AllModsOneIsNterminus.OrderBy(kv => kv.Key).Select(kv => (kv.Key, kv.Value.IdWithMotif)).ToList()));
                Assert.That(read[i].MonoisotopicMass, Is.EqualTo(written[i].MonoisotopicMass).Within(1e-9));
            }
            Directory.Delete(folder, true);
        }

        private static Modification Known(string type, string idWithMotif) =>
            GlobalVariables.AllModsKnown.First(m => m.ModificationType == type && m.IdWithMotif == idWithMotif);
    }
}
