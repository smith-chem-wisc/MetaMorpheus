using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// A glyco PSM is labelled like any other PSM: per candidate, through <see cref="DecoyContaminantTargetLabel"/>, so a
    /// match to an entrapment protein is written <c>ET</c>, not <c>T</c> (MetaMorpheus #2901, the glyco half).
    /// </summary>
    [TestFixture]
    public static class TestOGlycoEntrapmentLabel
    {
        private static string Label(GlycoSpectralMatch gsm)
        {
            int column = Array.IndexOf(GlycoSpectralMatch.GetTabSepHeaderSingle().Split('\t'), SpectrumMatchFromTsvHeader.DecoyContaminantTarget);
            Assert.That(column, Is.GreaterThanOrEqualTo(0), "premise: the glyco header has a Decoy/Contaminant/Target column");
            return gsm.SingleToString().Split('\t')[column];
        }

        private static PeptideWithSetModifications PeptideOn(string accession, bool isDecoy = false, bool isContaminant = false, bool isEntrapment = false)
        {
            var protein = new Protein("PEPTIDESTK", accession, isDecoy: isDecoy, isContaminant: isContaminant, isEntrapment: isEntrapment);
            return new PeptideWithSetModifications("PEPTIDESTK", null, p: protein, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 10);
        }

        private static GlycoSpectralMatch MatchTo(params PeptideWithSetModifications[] candidates)
        {
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 2, 2, true, Polarity.Positive,
                double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=2", double.NaN, null, null, double.NaN, null,
                DissociationType.HCD, 1, null), 1000, 2, null, new CommonParameters());
            var gsm = new GlycoSpectralMatch(candidates[0], 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            foreach (var other in candidates.Skip(1))
                gsm.AddOrReplace(other, 10, 0, true, new List<MatchedFragmentIon>());
            gsm.ResolveAllAmbiguities(); // as the glyco engine does before writing (GlycoSearchEngine.cs:494); it sets IsDecoy/IsContaminant
            return gsm;
        }

        [Test]
        [TestCase(false, false, false, "T")]
        [TestCase(true, false, false, "D")]
        [TestCase(false, true, false, "C")]
        [TestCase(false, false, true, "ET")]
        [TestCase(true, false, true, "ED")]
        public static void GlycoPsm_IsLabelledByItsProtein(bool isDecoy, bool isContaminant, bool isEntrapment, string expected)
        {
            var gsm = MatchTo(PeptideOn("P1", isDecoy, isContaminant, isEntrapment));
            Assert.That(Label(gsm), Is.EqualTo(expected));
        }

        /// <summary>
        /// A PSM shared by a target and an entrapment protein is labelled per candidate, as PsmTsvWriter does, so a reader
        /// can count it as half an entrapment (<see cref="DecoyContaminantTargetLabel.EntrapmentFraction"/>).
        /// </summary>
        [Test]
        public static void GlycoPsm_SharedByATargetAndAnEntrapmentProtein_IsLabelledPerCandidate()
        {
            var gsm = MatchTo(PeptideOn("P1"), PeptideOn("Random_P1", isEntrapment: true));
            Assert.That(gsm.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2), "premise: both candidates are kept");
            Assert.That(Label(gsm), Is.EqualTo("T|ET"));
        }

        /// <summary>
        /// End to end: the O-glyco search on the leukosialin spectrum (TTGSLEPSSGASGPQVSSVK, see
        /// <see cref="GlycoSemiSpecificSearchTests"/>), with leukosialin's accession marking it entrapment or not. The mzLib
        /// loader treats an accession containing "Random" as entrapment.
        /// </summary>
        [Test]
        [TestCase("P16150", "T")]
        [TestCase("Random_P16150", "ET")]
        public static void OGlycoSearch_WritesTheEntrapmentLabel(string accession, string expected)
        {
            string dataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(TestOGlycoEntrapmentLabel)}_{accession}");
            string fastaPath = outputFolder + ".fasta";
            try
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                Directory.CreateDirectory(outputFolder);
                string fasta = File.ReadAllText(Path.Combine(dataFolder, "P16150.fasta"));
                Assert.That(fasta, Does.StartWith(">sp|P16150|"), "premise: the FASTA this test edits");
                File.WriteAllText(fastaPath, fasta.Replace(">sp|P16150|", $">sp|{accession}|"));

                var task = MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(Path.Combine(dataFolder, "GlycoSearchTaskconfigOGlycoTest_Run.toml"));
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) },
                    new List<string> { Path.Combine(dataFolder, "2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf") },
                    new List<DbForTask> { new DbForTask(fastaPath, false) },
                    outputFolder).Run();

                string oGlycoPath = Path.Combine(outputFolder, "Task", "oglyco.psmtsv");
                Assert.That(File.Exists(oGlycoPath), Is.True, "premise: the search identifies the glycopeptide");
                var psms = SpectrumMatchTsvReader.ReadGlycoPsmTsv(oGlycoPath, out var errors);
                Assert.That(errors, Is.Empty);
                var psm = psms.Single(p => p.BaseSeq == "TTGSLEPSSGASGPQVSSVK");
                Assert.That(psm.Accession, Is.EqualTo(accession), "premise: the PSM is on the edited entry");
                Assert.That(psm.DecoyContamTarget, Is.EqualTo(expected));
            }
            finally
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                if (File.Exists(fastaPath))
                    File.Delete(fastaPath);
            }
        }
    }
}
