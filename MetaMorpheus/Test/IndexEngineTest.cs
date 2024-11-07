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
                //Assert.Contains(peptide, results.PeptideIndex);

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
                    List<int> fragmentBin = results.FragmentIndex[integerMassRepresentation];

                    // this list should contain this peptide!
                    Assert.That(fragmentBin.Contains(positionInPeptideIndex));
                    //Assert.Contains(positionInPeptideIndex, fragmentBin);
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
                    List<int> fragmentBin = results.FragmentIndex[integerMassRepresentation];

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