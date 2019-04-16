using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using MzLibUtil;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal class VariantSearchTests
    {
        private static CommonParameters CommonParameters = new CommonParameters(
            scoreCutoff: 1,
            digestionParams: new DigestionParams(
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                maxModsForPeptides: 1));

        [Test]
        [TestCase(0, 0, true, "P4V")] // variant is in the detected peptide
        [TestCase(1, 1, false, "PT4KT")] // variant intersects psm, but isn't really identified because PEK/TIDE is the same as PEP/TIDE if you only detect the second peptide
        [TestCase(2, 0, true, "P4PPP")] // intersecting sequence between variant and detected peptide is smaller than the original sequence, so clearly identied
        [TestCase(3, 0, true, "PPP4P")] // peptide is different than original sequence, but the variant site is the same AA
        [TestCase(4, 0, false, "PKPK4PK")]
        [TestCase(5, 1, true, "PTA4KT")] // counterpoint to (1), where the second peptide does distinguish
        [TestCase(6, 0, false, "KKA4K")]
        [TestCase(7, 1, true, "P4V[type:mod on V]")]
        [TestCase(8, 1, true, "P4PP[type:mod on P]P")]
        [TestCase(0, 0, true, "P6V", DecoyType.Reverse)] // variant is in the detected decoy peptide MEDITVEP
        [TestCase(2, 0, true, "P6PPP", DecoyType.Reverse)] // intersecting sequence between variant and detected peptide is smaller than the original sequence, so clearly identied
        [TestCase(3, 0, true, "PPP6P", DecoyType.Reverse)] // peptide is different than original sequence, but the variant site is the same AA
        [TestCase(7, 1, true, "P6V[type:mod on V]", DecoyType.Reverse)]
        [TestCase(8, 1, true, "P6PP[type:mod on P]P", DecoyType.Reverse)]
        [TestCase(9, 0, true, "PTIDEPEPTIDE4PPP")] // intersecting sequence between variant and detected peptide is smaller than the original sequence, so clearly identied
        [TestCase(9, 0, true, "EDITPEPEDITP2PPP", DecoyType.Reverse)] // MPEPPP becomes MPPPEP with the variant beginning at position 2
        public static void SearchTests(int proteinIdx, int peptideIdx, bool containsVariant, string variantPsmShort, DecoyType decoyType = DecoyType.None)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            // Make sure can run the complete search task when multiple compact peptides may correspond to a single PWSM
            SearchTask st = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DecoyType = decoyType,
                    SearchTarget = decoyType == DecoyType.None,
                    ModPeptidesAreDifferent = false
                },
                CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 2), precursorMassTolerance: new PpmTolerance(20)),
            };

            ModificationMotif.TryGetMotif("V", out ModificationMotif motifV);
            Modification mv = new Modification("mod", null, "type", null, motifV, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);

            List<Protein> proteins = new List<Protein>
            {
                new Protein("MPEPTIDE", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "PT", "KT", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPPPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PPP", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPKPKTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 7, "PKPK", "PK", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTAIDE", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "PTA", "KT", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEKKAIDE", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "KKA", "K", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 4, new[] { mv }.ToList() } }) }),
                new Protein("MPEPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 5, new[] { mp }.ToList() } }) }),
                new Protein("MPEPTIDEPEPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "PTIDEPEPTIDE", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
            };
            PeptideWithSetModifications pep = proteins[proteinIdx].GetVariantProteins().SelectMany(p => p.Digest(CommonParameters.DigestionParams, null, null)).ToList()[peptideIdx];

            string xmlName = $"andguiaheov{proteinIdx.ToString()}.xml";
            ProteinDbWriter.WriteXmlDatabase(null, new List<Protein> { proteins[proteinIdx] }, xmlName);

            string mzmlName = $"ajgdiv{proteinIdx.ToString()}.mzML";
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pep });

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"TestSearchWithVariants{proteinIdx.ToString()}");
            Directory.CreateDirectory(outputFolder);

            st.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "");
            var psms = File.ReadAllLines(Path.Combine(outputFolder, "AllPSMs.psmtsv"));
            
            Assert.IsTrue(psms.Any(line => line.Contains($"\t{variantPsmShort}\t" + (containsVariant ? variantPsmShort : "\t"))));

            Directory.Delete(outputFolder, true);
            File.Delete(mzmlName);
            File.Delete(xmlName);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);

            Console.WriteLine($"Analysis time for VariantSearchTests.SearchTests({proteinIdx.ToString()},{peptideIdx.ToString()},{containsVariant.ToString()},{variantPsmShort}): {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        [TestCase("frameshift.xml")]
        [TestCase("frameshift.xml", DecoyType.Reverse)]
        public void MoreTests(string filename, DecoyType decoyType = DecoyType.None)
        {
            string xmlName = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", filename);
            var proteins = ProteinDbLoader.LoadProteinXML(xmlName, decoyType == DecoyType.None, decoyType, null, false, null, out var un);
            var peps = proteins[1].Digest(CommonParameters.DigestionParams, null, null).ToList();
            PeptideWithSetModifications pep = peps[peps.Count - 2];

            string mzmlName = $"ajgdiv{filename}{decoyType.ToString()}.mzML";
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pep });

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"TestSearchWithVariants{filename}{decoyType.ToString()}");
            Directory.CreateDirectory(outputFolder);

            SearchTask st = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DecoyType = decoyType,
                    SearchTarget = decoyType == DecoyType.None,
                    ModPeptidesAreDifferent = false
                },
                CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 2), precursorMassTolerance: new PpmTolerance(20)),
            };

            st.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "");
            var psms = File.ReadAllLines(Path.Combine(outputFolder, "AllPSMs.psmtsv"));

            //Assert.IsTrue(psms.Any(line => line.Contains($"\t{variantPsmShort}\t" + (containsVariant ? variantPsmShort : "\t"))));

            Directory.Delete(outputFolder, true);
            File.Delete(mzmlName);
            //Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }
    }
}