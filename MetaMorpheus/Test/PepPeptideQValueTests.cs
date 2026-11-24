using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// PEP trains at the peptide level once there are 100 distinct full sequences, and labels its positives by
    /// peptide q-value. When there are more PSMs than <see cref="FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue"/>
    /// but not more peptides, those peptide q-values must still be computed before PEP, or no target is a positive.
    /// Every peptide-level value sits on each peptide's best-scoring match, the one PEP trains on
    /// (<see cref="SpectralMatchGroup.GetBestMatches"/>), never on its lowest-PEP match.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // sets FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue, a process-wide static
    public static class PepPeptideQValueTests
    {
        [Test]
        public static void PepTrainsWhenThereAreMorePsmsThanTheThresholdButFewerPeptides()
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDirectory, "PepPeptideQValueTests");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var task = Nett.Toml.ReadFile<SearchTask>(Path.Combine(testDirectory, @"TestData\Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            task.CommonParameters.QValueCutoffForPepCalculation = 0.01;
            var spectraFiles = new List<string>
            {
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML"),
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML")
            };
            var databases = new List<DbForTask> { new DbForTask(Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta"), false) };

            var threshold = typeof(FdrAnalysisEngine).GetProperty(nameof(FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue))!;
            int original = FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                try
                {
                    threshold.SetValue(null, 500); // 770 PSMs, 308 peptides
                    new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", task) }, spectraFiles, databases, outputFolder).Run();
                }
                finally
                {
                    threshold.SetValue(null, original);
                }
            }

            string results = File.ReadAllText(Path.Combine(outputFolder, "search", "results.txt"));
            Assert.That(results, Does.Not.Contain("Posterior error probability analysis failed"));
            Assert.That(results, Does.Contain("Metrics for Determination of PEP"));
            Directory.Delete(outputFolder, true);
        }

        private static readonly CommonParameters Common = new();
        private const string Residues = "ACDEFGHILNQSTVWY";

        /// <summary>A tryptic peptide of its own protein, one per index.</summary>
        private static PeptideWithSetModifications Peptide(int index, bool decoy)
        {
            var letters = new char[7];
            for (int i = 0, n = index; i < letters.Length; i++, n /= Residues.Length)
                letters[i] = Residues[n % Residues.Length];
            string sequence = "G" + new string(letters) + "K";
            var protein = new Protein(sequence, (decoy ? "DECOY_" : "") + "P" + index, isDecoy: decoy);
            return protein.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .Single(p => p.BaseSequence == sequence);
        }

        private static SpectralMatch Psm(PeptideWithSetModifications peptide, double score, double pep, int scanNumber)
        {
            var spectrum = new MzSpectrum(new double[] { 100, 200 }, new double[] { 1, 1 }, false);
            var scan = new MsDataScan(spectrum, scanNumber, 2, true, Polarity.Positive, scanNumber / 10.0, null, "", MZAnalyzerType.Orbitrap, 2, null, null, null);
            var ms2 = new Ms2ScanWithSpecificMass(scan, peptide.MonoisotopicMass.ToMz(2), 2, "file", Common);
            var ions = new List<MatchedFragmentIon> { new(new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0), 101, 10, 1) };
            var psm = new PeptideSpectralMatch(peptide, 0, score, scanNumber, ms2, Common, ions);
            psm.ResolveAllAmbiguities();
            psm.PsmFdrInfo = new FdrInfo { PEP = pep };
            psm.PeptideFdrInfo = new FdrInfo();
            return psm;
        }

        /// <summary>
        /// Ten targets scored 40 to 31 at PEP 0.01 and five decoys scored 10 to 6, two at PEP 0.1 and three at 0.9. Peptide
        /// X has a best-scoring match (score 50, PEP 0.3) and a weak match with a lower PEP (score 1, PEP 0.001): by score
        /// its best match is above every decoy, by PEP below two of them.
        /// </summary>
        private static (List<SpectralMatch> All, SpectralMatch XBest, SpectralMatch XLowestPep) PeptideX(bool withXLowestPep = true)
        {
            var all = new List<SpectralMatch>();
            int scan = 1;
            for (int t = 0; t < 10; t++)
                all.Add(Psm(Peptide(t, decoy: false), 40 - t, 0.01, scan++));
            for (int d = 0; d < 5; d++)
                all.Add(Psm(Peptide(50 + d, decoy: true), 10 - d, d < 2 ? 0.1 : 0.9, scan++));
            var x = Peptide(100, decoy: false);
            var xBest = Psm(x, 50, 0.3, scan++);
            all.Add(xBest);
            var xLowestPep = withXLowestPep ? Psm(x, 1, 0.001, scan++) : null;
            if (xLowestPep != null) all.Add(xLowestPep);
            return (all, xBest, xLowestPep);
        }

        private static List<SpectralMatch> PeptideRows(List<SpectralMatch> psms, bool includeHighQValuePsms = true) =>
            FilteredPsms.Filter(psms, Common, includeDecoys: true, includeContaminants: true, includeAmbiguous: true,
                includeHighQValuePsms: includeHighQValuePsms, filterAtPeptideLevel: true).FilteredPsmsList;

        /// <summary>
        /// The FDR passes over PEPs already computed (doPep false): X's row, q-value and PEP q-value come from its
        /// best-scoring match, the same as when that match is X's only one, and its lowest-PEP match carries none.
        /// </summary>
        [Test]
        public static void ThePeptideLevelValuesSitOnEachPeptidesBestScoringMatch()
        {
            var (all, xBest, xLowestPep) = PeptideX();
            FdrAnalysisEngine.DoFalseDiscoveryRateAnalysis(all, false, new List<(string, CommonParameters)>(), null, null, Common);
            var (alone, xBestAlone, _) = PeptideX(withXLowestPep: false);
            FdrAnalysisEngine.DoFalseDiscoveryRateAnalysis(alone, false, new List<(string, CommonParameters)>(), null, null, Common);

            Assert.That(PeptideRows(all).Single(p => p.FullSequence == xBest.FullSequence), Is.SameAs(xBest));
            Assert.That(xBest.PeptideFdrInfo.QValue, Is.EqualTo(xBestAlone.PeptideFdrInfo.QValue));
            Assert.That(xBest.PeptideFdrInfo.PEP_QValue, Is.EqualTo(xBestAlone.PeptideFdrInfo.PEP_QValue));
            Assert.That((xLowestPep.PeptideFdrInfo.QValue, xLowestPep.PeptideFdrInfo.PEP_QValue), Is.EqualTo((2.0, 2.0)));
            Assert.That(PeptideRows(all, includeHighQValuePsms: false), Does.Contain(xBest), "X passes at 1% on its best score");
        }

        /// <summary>
        /// The branch with more matches than the threshold but no more peptides computes peptide PEP q-values too, on
        /// the same matches. These sixteen peptides are too few to train on, so the fixture's PEPs stand.
        /// </summary>
        [Test]
        public static void WithFewPeptides_ThePeptidePepQValueSitsOnTheBestScoringMatch()
        {
            var (all, xBest, xLowestPep) = PeptideX();
            int peptides = all.Select(p => p.FullSequence).Distinct().Count();

            var threshold = typeof(FdrAnalysisEngine).GetProperty(nameof(FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue))!;
            int original = FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                try
                {
                    threshold.SetValue(null, peptides);
                    new FdrAnalysisEngine(all, 1, Common, new List<(string, CommonParameters)> { ("file", Common) }, new List<string>()).Run();
                }
                finally
                {
                    threshold.SetValue(null, original);
                }
            }
            Assume.That((xBest.FdrInfo.PEP, xLowestPep.FdrInfo.PEP), Is.EqualTo((0.3, 0.001)), "the fixture's PEPs must stand");

            Assert.That(xBest.PeptideFdrInfo.PEP_QValue, Is.LessThan(2));
            Assert.That(xLowestPep.PeptideFdrInfo.PEP_QValue, Is.EqualTo(2));
        }

        /// <summary>
        /// PEP training's stopping count ranks X by its best-scoring match's PEP (0.3), below two decoys, so X is not
        /// accepted at 0.2; on its lowest-PEP match (0.001) it would be.
        /// </summary>
        [Test]
        public static void PepTrainingCountsAcceptedPeptidesOnTheBestScoringMatches()
        {
            var (all, _, _) = PeptideX();
            Assert.That(PepAnalysisEngine.CountAcceptedPeptides(all, 0.2), Is.EqualTo(10));
        }
    }
}
