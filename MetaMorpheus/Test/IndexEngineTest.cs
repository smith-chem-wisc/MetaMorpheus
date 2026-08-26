using EngineLayer;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Digestion;
using Omics.Modifications;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class IndexEngineTest
    {
        /// <summary>
        /// The fragment index is built in parallel over contiguous blocks of peptide ids, and the number
        /// of blocks follows MaxThreadsToUsePerFile. The index that comes out has to be byte-for-byte the
        /// same regardless -- peptide ids must stay ascending within a bin, because the search binary
        /// searches a bin by peptide mass and that only holds while ids ascend.
        /// </summary>
        [Test]
        public static void TestIndexIsIdenticalRegardlessOfThreadCount()
        {
            var proteinList = new List<Protein>
            {
                new Protein("MNNNKQQQMNNNKQQQPEPTIDEKMSSSRTTTKAAAWWWKGGGYYYK", "prot1"),
                new Protein("MKVLINGYGTIGKRVADAVSQQDDMKVIGVSKTRPDFEARMALQKGYDLYVAIPK", "prot2"),
                new Protein("MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAK", "prot3"),
            };

            var fixedModifications = new List<Modification>();
            var variableModifications = new List<Modification>();

            MassBinIndex referenceFragmentIndex = null;
            List<PeptideWithSetModifications> referencePeptideIndex = null;

            foreach (int threadCount in new[] { 1, 2, 3, 8 })
            {
                var commonParameters = new CommonParameters(maxThreadsToUsePerFile: threadCount);

                var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 0,
                    DecoyType.Reverse, commonParameters, null, 30000, false, new List<FileInfo>(),
                    TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

                var results = (IndexingResults)engine.Run();

                if (referenceFragmentIndex == null)
                {
                    referenceFragmentIndex = results.FragmentIndex;
                    referencePeptideIndex = results.PeptideIndex;
                    Assert.That(referenceFragmentIndex.EntryCount, Is.GreaterThan(0), "index came out empty; the test proves nothing");
                    continue;
                }

                Assert.That(results.PeptideIndex.Count, Is.EqualTo(referencePeptideIndex.Count), $"peptide count differs at {threadCount} threads");
                Assert.That(results.FragmentIndex.Length, Is.EqualTo(referenceFragmentIndex.Length), $"bin count differs at {threadCount} threads");
                Assert.That(results.FragmentIndex.EntryCount, Is.EqualTo(referenceFragmentIndex.EntryCount), $"entry count differs at {threadCount} threads");

                for (int bin = 0; bin < referenceFragmentIndex.Length; bin++)
                {
                    if (referenceFragmentIndex.CountInBin(bin) == 0 && results.FragmentIndex.CountInBin(bin) == 0)
                    {
                        continue;
                    }

                    Assert.That(results.FragmentIndex[bin].ToArray(), Is.EqualTo(referenceFragmentIndex[bin].ToArray()),
                        $"bin {bin} differs at {threadCount} threads");
                }
            }
        }

        /// <summary>
        /// Peptide ids ascend within every bin. BinarySearchBinForPrecursorIndex depends on it.
        /// </summary>
        [Test]
        public static void TestFragmentIndexBinsAreSortedByPeptideId()
        {
            var proteinList = new List<Protein>
            {
                new Protein("MNNNKQQQMNNNKQQQPEPTIDEKMSSSRTTTKAAAWWWKGGGYYYK", "prot1"),
                new Protein("MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAK", "prot2"),
            };

            var engine = new IndexingEngine(proteinList, new List<Modification>(), new List<Modification>(), null, null, null, 0,
                DecoyType.Reverse, new CommonParameters(), null, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)engine.Run();

            int binsChecked = 0;
            for (int bin = 0; bin < results.FragmentIndex.Length; bin++)
            {
                ReadOnlySpan<int> ids = results.FragmentIndex[bin];
                for (int i = 1; i < ids.Length; i++)
                {
                    Assert.That(ids[i], Is.GreaterThan(ids[i - 1]), $"bin {bin} is not ascending");
                }
                if (ids.Length > 1)
                {
                    binsChecked++;
                }
            }

            Assert.That(binsChecked, Is.GreaterThan(0), "no multi-entry bins; the test proves nothing");
        }

        [Test]
        public static void TestIndexEngine()
        {
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease p = new Protease("Custom Protease2", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(p.Name, p);
            CommonParameters CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: new DigestionParams(protease: p.Name, minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", CommonParameters));

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.None, CommonParameters,
                fsp, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.That(results.PeptideIndex.Count, Is.EqualTo(5));

            var digestedList = proteinList[0].Digest(CommonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();

            Assert.That(digestedList.Count, Is.EqualTo(5));
            foreach (PeptideWithSetModifications peptide in digestedList)
            {
                Assert.That(results.PeptideIndex.Contains(peptide));

                var fragments = new List<Product>();
                peptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, fragments);

                int positionInPeptideIndex = results.PeptideIndex.IndexOf(peptide);

                foreach (Product fragment in fragments)
                {
                    // mass of the fragment
                    double fragmentMass = fragment.NeutralMass;
                    int integerMassRepresentation = (int)Math.Round(fragmentMass * 1000);

                    // look up the peptides that have fragments with this mass
                    // the result of the lookup is a list of peptide IDs that have this fragment mass
                    ReadOnlySpan<int> fragmentBin = results.FragmentIndex[integerMassRepresentation];

                    // this list should contain this peptide!
                    Assert.That(fragmentBin.Contains(positionInPeptideIndex));
                }
            }
        }

        [Test]
        public static void TestIndexEngineWithWeirdSeq()
        {
            var proteinList = new List<Protein> { new Protein("MQXQ", null) };
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
            {
                modsDictionary.Add(mod, 0);
            }
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease protease = new Protease("Custom Protease", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                scoreCutoff: 1);
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", CommonParameters));
            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, CommonParameters,
                fsp, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.That(results.PeptideIndex.Count, Is.EqualTo(1));

            Assert.That(results.PeptideIndex[0].MonoisotopicMass, Is.NaN);
            Assert.That(results.FragmentIndex.Length, Is.EqualTo(30000000 + 1));
        }

        [Test]
        public static void TestIndexEngineLowRes()
        {
            var proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"indexEngineTestFasta.fasta"), true, DecoyType.Reverse, false, out var dbErrors,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, -1);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            CommonParameters CommonParameters = new CommonParameters(dissociationType: DissociationType.LowCID, maxThreadsToUsePerFile: 1, scoreCutoff: 1, digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", CommonParameters));
            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, CommonParameters,
                fsp, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.That(results.PeptideIndex.Count, Is.EqualTo(10));

            var bubba = results.FragmentIndex;
            var tooBubba = results.PrecursorIndex;


            var digestedList = proteinList[0].Digest(CommonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();
            digestedList.AddRange(proteinList[1].Digest(CommonParameters.DigestionParams, new List<Modification>(), variableModifications));

            Assert.That(digestedList.Count, Is.EqualTo(10));
            foreach (PeptideWithSetModifications peptide in digestedList)
            {
                Assert.That(results.PeptideIndex.Contains(peptide));

                var fragments = new List<Product>();
                peptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, fragments);

                int positionInPeptideIndex = results.PeptideIndex.IndexOf(peptide);

                foreach (Product fragment in fragments.Where(f => f.ProductType == ProductType.b || f.ProductType == ProductType.y))
                {
                    // mass of the fragment
                    double fragmentMass = Math.Round(fragment.NeutralMass / 1.0005079, 0) * 1.0005079;
                    int integerMassRepresentation = (int)Math.Round(fragmentMass * 1000);

                    // look up the peptides that have fragments with this mass
                    // the result of the lookup is a list of peptide IDs that have this fragment mass
                    ReadOnlySpan<int> fragmentBin = results.FragmentIndex[integerMassRepresentation];

                    // this list should contain this peptide!
                    Assert.That(fragmentBin.Contains(positionInPeptideIndex));
                }
            }
            foreach (var fdfd in digestedList)
            {
                Assert.That(results.PeptideIndex.Contains(fdfd));
            }
        }
    }
}