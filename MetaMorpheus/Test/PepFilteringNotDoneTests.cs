using EngineLayer;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The PEP q-value filter chosen, but PEP not trained for some matches: those matches are judged by their q-value
    /// notch at the same threshold, the filter says so, and it never drops every match because PEP was not
    /// trained. Also a filter a caller asks for by name.
    /// </summary>
    [TestFixture]
    public static class PepFilteringNotDoneTests
    {
        // The PEP q-value filter, at 0.01
        private static CommonParameters PepChosen => new CommonParameters(identificationFilter: FilterType.PepQValue, filterThreshold: 0.01);

        /// <summary>
        /// <paramref name="count"/> target matches; the first <paramref name="confident"/> have q, notch and PEP
        /// q-value 0.001, the rest 0.5. Without trained PEP, every PEP q-value is FdrInfo's initial 2.
        /// </summary>
        private static List<SpectralMatch> Matches(int count, int confident, bool pepTrained)
        {
            var peptide = new Protein("PEPTIDEK", "target1")
                .Digest(new DigestionParams(protease: "trypsin", minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First();
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());

            var psms = new List<SpectralMatch>();
            for (int i = 0; i < count; i++)
            {
                var psm = new PeptideSpectralMatch(peptide, 0, 10, i, scan, new CommonParameters(), new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                double q = i < confident ? 0.001 : 0.5;
                psm.SetFdrValues(i + 1, 0, q, i + 1, 0, q, pep: pepTrained ? q : 0, pepQValue: pepTrained ? q : 2);
                psm.PsmFdrInfo.PepTrained = pepTrained;
                psms.Add(psm);
            }
            return psms;
        }

        [Test]
        public static void PepNotTrained_FiltersByTheNotchAndSaysSo_InsteadOfDroppingEveryMatch()
        {
            // 100 to 1000 matches: FdrAnalysisEngine trains PEP only above 1000 per digestion agent, so every PEP
            // q-value is still 2. Once, the PEP q-value filter was used anyway and kept none of them.
            var filtered = FilteredPsms.Filter(Matches(300, 200, pepTrained: false), PepChosen);

            Assert.That((filtered.FilterType, filtered.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(200));
            Assert.That(filtered.PepFilteringNotDone, Is.EqualTo(
                "PEP filtering was requested, but PEP was not trained for 300 target PSMs (trypsin: 300), so those were filtered by q-value notch <= 0.01."));
        }

        [Test]
        public static void PepTrainedForPsmsButNotPeptides_SaysPeptideLevel()
        {
            // More than 1000 PSMs but at most 1000 peptides: FdrAnalysisEngine trains PEP on the PSMs, and the
            // peptide-level PEP q-values are not computed. The PSM table filters on PEP; the peptide table falls back to
            // the notch and says so for the peptides.
            var psms = Matches(300, 200, pepTrained: true);
            foreach (var psm in psms)
            {
                var psmFdr = psm.GetFdrInfo(false);
                psm.PeptideFdrInfo = new FdrInfo { QValue = psmFdr.QValue, QValueNotch = psmFdr.QValueNotch, PEP = 0, PEP_QValue = 2, PepTrained = false };
            }

            var psmLevel = FilteredPsms.Filter(psms, PepChosen);
            var peptideLevel = FilteredPsms.Filter(psms, PepChosen, filterAtPeptideLevel: true);

            Assert.That((psmLevel.FilterType, psmLevel.PepFilteringNotDone), Is.EqualTo((FilterType.PepQValue, (string)null)));
            Assert.That(peptideLevel.FilteredPsmsList, Is.Not.Empty, "the peptide is kept on its notch, not dropped on PEP q = 2");
            Assert.That(peptideLevel.PepFilteringNotDone, Is.EqualTo(
                "PEP filtering was requested, but peptide-level PEP q-values were not computed for 1 target peptide (trypsin: 1), so those were filtered by q-value notch <= 0.01."));
        }

        [Test]
        public static void PepTrained_FiltersOnPepAndSaysNothing()
        {
            var filtered = FilteredPsms.Filter(Matches(300, 150, pepTrained: true), PepChosen);

            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(150));
            Assert.That(filtered.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void DefaultSettings_FilterByThePepQValue_AndSayWherePepWasNotTrained()
        {
            // Bottom-up, the default filter is the PEP q-value: untrained, every match is judged by its notch, and it says so.
            var filtered = FilteredPsms.Filter(Matches(300, 200, pepTrained: false), new CommonParameters());

            Assert.That((filtered.FilterType, filtered.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(200));
            Assert.That(filtered.PepFilteringNotDone, Is.EqualTo(
                "PEP filtering was requested, but PEP was not trained for 300 target PSMs (trypsin: 300), so those were filtered by q-value notch <= 0.01."));

            // Top-down, the default is the q-value notch, which says nothing
            var topDown = FilteredPsms.Filter(Matches(300, 200, pepTrained: false), new CommonParameters(digestionParams: new DigestionParams(protease: "top-down")));
            Assert.That((topDown.FilterType, topDown.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));
            Assert.That(topDown.FilteredPsmsList, Has.Count.EqualTo(200));
            Assert.That(topDown.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void TheNotchAlone_Decides_AndThePooledQValueAlone_Decides()
        {
            // A notch q-value below the threshold and a whole-search q-value above it: the notch filter keeps the match and
            // the pooled q-value filter does not.
            var psms = Matches(300, 0, pepTrained: false);
            foreach (var psm in psms.Take(120)) psm.SetFdrValues(1, 0, 0.5, 1, 0, 0.001, 0, 2);

            var byPooled = FilteredPsms.Filter(psms, new CommonParameters(), filterType: FilterType.QValue);
            Assert.That(byPooled.FilteredPsmsList, Is.Empty);
            Assert.That((byPooled.FilterType, byPooled.GetFilterTypeString()), Is.EqualTo((FilterType.QValue, "q-value")));

            var byNotch = FilteredPsms.Filter(psms, new CommonParameters(identificationFilter: FilterType.QValueNotch));

            Assert.That(byNotch.FilteredPsmsList, Has.Count.EqualTo(120));
            Assert.That((byNotch.FilterType, byNotch.FilterThreshold, byNotch.GetFilterTypeString()), Is.EqualTo((FilterType.QValueNotch, 0.01, "q-value notch")));
            Assert.That(byNotch.GetFilterValue(psms[0]), Is.EqualTo(0.001), "the value handed to quantification is the notch q-value");
            Assert.That(byNotch.TargetPsmsAboveThreshold, Is.EqualTo(120));
        }

        [Test]
        public static void AFilterAskedForByName_IsHeldToTheOneThreshold()
        {
            // There is one threshold: a filter asked for by name is held to the search's FilterThreshold
            var filtered = FilteredPsms.Filter(Matches(50, 20, pepTrained: false), PepChosen, filterType: FilterType.QValueNotch);

            Assert.That(filtered.FilterThreshold, Is.EqualTo(0.01), "the notch is held to the search's one threshold");
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(20));
            Assert.That(filtered.PepFilteringNotDone, Is.Null);

            var atOne = FilteredPsms.Filter(Matches(50, 20, pepTrained: false),
                new CommonParameters(identificationFilter: FilterType.PepQValue, filterThreshold: 1.0), filterType: FilterType.QValueNotch);
            Assert.That((atOne.FilterThreshold, atOne.FilteredPsmsList.Count), Is.EqualTo((1.0, 50)));
        }

        /// <summary>
        /// The whole search, with the PEP q-value filter chosen: the two TaGe snips give between 100 and 1000 matches, so
        /// FdrAnalysisEngine never trains PEP. Once every match, peptide and protein group was dropped. Now the matches
        /// are filtered by the notch at 0.01, results.txt and the warning say so, and protein groups are held to the same
        /// 0.01.
        /// </summary>
        [Test]
        public static void SearchTask_PepNotTrained_WritesNotchFilteredResultsAndSaysSo()
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDirectory, "PepFilteringNotDone_SearchTask");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var task = new SearchTask { CommonParameters = PepChosen };
            var spectraFiles = new List<string>
            {
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML"),
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML")
            };
            var databases = new List<DbForTask> { new DbForTask(Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta"), false) };
            var warnings = new List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);

            // QvalueThresholdOverride is process-wide and would make the FDR engine train PEP on any count.
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                MetaMorpheusTask.WarnHandler += onWarn;
                try
                {
                    typeof(FdrAnalysisEngine).GetProperty("QvalueThresholdOverride")!.SetValue(null, false);
                    new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", task) }, spectraFiles, databases, outputFolder).Run();
                }
                finally
                {
                    MetaMorpheusTask.WarnHandler -= onWarn;
                }
            }

            string searchFolder = Path.Combine(outputFolder, "search");
            string[] results = File.ReadAllLines(Path.Combine(searchFolder, "results.txt"));
            const string said = "PEP filtering was requested, but PEP was not trained for 258 target PSMs (trypsin: 258), so those were filtered by q-value notch <= 0.01.";
            const string saidForPeptides = "PEP filtering was requested, but peptide-level PEP q-values were not computed for 88 target peptides (trypsin: 88), so those were filtered by q-value notch <= 0.01.";
            int first = Array.IndexOf(results, said);
            Assert.That(first, Is.GreaterThan(-1));
            Assert.That(results.Skip(first).Take(5), Is.EqualTo(new[]
            {
                said, // just above the PSM count it qualifies
                "All target PSMs with pep q-value <= 0.01: 258",
                saidForPeptides, // just above the peptide count
                "All target peptides with pep q-value <= 0.01: 88",
                "All target protein groups with pep q-value <= 0.01 (1% FDR): 84",
            }));
            Assert.That(warnings.Where(w => w.StartsWith("PEP filtering was requested")), Is.EqualTo(new[] { said }), "one warning per search");

            // Protein groups: held to the filter's 0.01
            string[] groups = File.ReadAllLines(Directory.GetFiles(searchFolder, "All*ProteinGroups.tsv").Single());
            int qColumn = groups[0].Split('\t').ToList().IndexOf("Protein QValue");
            Assert.That(qColumn, Is.Not.EqualTo(-1));
            var groupQValues = groups.Skip(1).Select(l => double.Parse(l.Split('\t')[qColumn], System.Globalization.CultureInfo.InvariantCulture)).ToList();
            Assert.That(groupQValues, Has.Count.EqualTo(84));
            Assert.That(groupQValues, Has.All.LessThanOrEqualTo(0.01));

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void GetFilterValue_IsTheValueTheMatchWasJudgedBy()
        {
            var psm = Matches(1, 0, pepTrained: false)[0];
            psm.SetFdrValues(1, 0, 0.3, 1, 0, 0.2, 0.05, 0.1);

            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.QValueNotch).GetFilterValue(psm), Is.EqualTo(0.2));
            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.PepQValue).GetFilterValue(psm), Is.EqualTo(0.2),
                "PEP not trained: the notch");
            psm.PsmFdrInfo.PepTrained = true;
            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.PepQValue).GetFilterValue(psm), Is.EqualTo(0.1));
        }

        [Test]
        public static void AFilterTypeWithNoCase_Throws()
        {
            // A FilterType added later without its case here must fail loudly, not be read as another value: the first match
            // it judges throws, and so does naming it.
            var psm = Matches(1, 0, pepTrained: false)[0];

            Assert.Throws<ArgumentOutOfRangeException>(() => FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: (FilterType)99));
            Assert.Throws<ArgumentOutOfRangeException>(() => FilteredPsms.Passes(psm, false, (FilterType)99, 0.01));
            Assert.Throws<ArgumentOutOfRangeException>(() => FilteredPsms.Describe((FilterType)99, 0.01));
        }
    }
}
