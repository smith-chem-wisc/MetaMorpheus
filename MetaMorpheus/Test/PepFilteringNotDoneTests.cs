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
    /// PEP filtering requested (a PEP q-value threshold below the q-value threshold) but not possible: the filter says
    /// what it did instead, and never drops every match because PEP was not computed. Also the q-value notch filter a
    /// caller can ask for by name.
    /// </summary>
    [TestFixture]
    public static class PepFilteringNotDoneTests
    {
        // What the GUI sends when its PEP box is ticked.
        private static CommonParameters GuiPepBox => new CommonParameters(qValueThreshold: 1.0, pepQValueThreshold: 0.01);

        /// <summary>
        /// <paramref name="count"/> target matches; the first <paramref name="confident"/> have q, notch and PEP
        /// q-value 0.001, the rest 0.5. A PEP q-value of 2 is FdrInfo's "never computed".
        /// </summary>
        private static List<SpectralMatch> Matches(int count, int confident, bool pepComputed)
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
                psm.SetFdrValues(i + 1, 0, q, i + 1, 0, q, pep: pepComputed ? q : 0, pepQValue: pepComputed ? q : 2);
                psms.Add(psm);
            }
            return psms;
        }

        [Test]
        public static void PepNotComputed_FiltersByQValueAndSaysSo_InsteadOfDroppingEveryMatch()
        {
            // 100 to 1000 matches: FdrAnalysisEngine trains PEP only above 1000 per digestion agent, so every PEP
            // q-value is still 2. Before, the PEP q-value filter was chosen anyway and kept none of them.
            var filtered = FilteredPsms.Filter(Matches(300, 200, pepComputed: false), GuiPepBox);

            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.QValue));
            Assert.That(filtered.FilterThreshold, Is.EqualTo(0.01));
            Assert.That(filtered.FilteringNotPerformed, Is.False);
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(200));
            Assert.That(filtered.PepFilteringNotDone,
                Is.EqualTo("PEP filtering was requested, but PEP was not computed for these PSMs, so results were filtered by q-value."));
        }

        [Test]
        public static void PepComputedForPsmsButNotPeptides_SaysPeptideLevel()
        {
            // More than 1000 PSMs but at most 1000 peptides: FdrAnalysisEngine trains PEP on the PSMs, and the
            // peptide-level PEP q-values stay 2. The PSM table filters on PEP; the peptide table falls back to the
            // q-value and must not say PEP was not computed for the PSMs.
            var psms = Matches(300, 200, pepComputed: true);
            foreach (var psm in psms)
            {
                var psmFdr = psm.GetFdrInfo(false);
                psm.PeptideFdrInfo = new FdrInfo { QValue = psmFdr.QValue, QValueNotch = psmFdr.QValueNotch, PEP = 0, PEP_QValue = 2 };
            }

            var psmLevel = FilteredPsms.Filter(psms, GuiPepBox);
            var peptideLevel = FilteredPsms.Filter(psms, GuiPepBox, filterAtPeptideLevel: true);

            Assert.That((psmLevel.FilterType, psmLevel.PepFilteringNotDone), Is.EqualTo((FilterType.PepQValue, (string)null)));
            Assert.That((peptideLevel.FilterType, peptideLevel.FilterThreshold, peptideLevel.FilteringNotPerformed), Is.EqualTo((FilterType.QValue, 0.01, false)));
            Assert.That(peptideLevel.FilteredPsmsList, Is.Not.Empty, "the peptide is kept on its q-value, not dropped on PEP q = 2");
            Assert.That(peptideLevel.PepFilteringNotDone,
                Is.EqualTo("PEP filtering was requested, but PEP q-values were not computed at the peptide level, so results were filtered by q-value."));
        }

        [Test]
        public static void FewerThan100_IsNotFilteredAndSaysSo()
        {
            var filtered = FilteredPsms.Filter(Matches(50, 20, pepComputed: false), GuiPepBox);

            Assert.That(filtered.FilteringNotPerformed, Is.True);
            Assert.That(filtered.FilterThreshold, Is.EqualTo(1));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(50), "behaviour unchanged: nothing is filtered");
            Assert.That(filtered.PepFilteringNotDone,
                Is.EqualTo("PEP filtering was requested, but there were fewer than 100 PSMs (50, counting decoys), so results were not filtered."));
        }

        [Test]
        public static void PepComputed_FiltersOnPepAndSaysNothing()
        {
            var filtered = FilteredPsms.Filter(Matches(300, 150, pepComputed: true), GuiPepBox);

            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(150));
            Assert.That(filtered.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void DefaultSettings_KeepTheQValueFilterAndSayNothing()
        {
            var filtered = FilteredPsms.Filter(Matches(300, 200, pepComputed: false), new CommonParameters());

            Assert.That((filtered.FilterType, filtered.FilterThreshold, filtered.FilteringNotPerformed), Is.EqualTo((FilterType.QValue, 0.01, false)));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(200));
            Assert.That(filtered.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void QValueNotch_ByName_FiltersOnTheNotchAlone()
        {
            // A notch q-value below the threshold and a whole-search q-value above it: the q-value filter (q AND notch)
            // drops the match, the notch filter keeps it.
            var psms = Matches(300, 0, pepComputed: false);
            foreach (var psm in psms.Take(120)) psm.SetFdrValues(1, 0, 0.5, 1, 0, 0.001, 0, 2);

            var byQValue = FilteredPsms.Filter(psms, new CommonParameters());
            var byNotch = FilteredPsms.Filter(psms, new CommonParameters(), filterType: FilterType.QValueNotch);

            Assert.That(byQValue.FilteredPsmsList, Is.Empty);
            Assert.That(byNotch.FilteredPsmsList, Has.Count.EqualTo(120));
            Assert.That((byNotch.FilterType, byNotch.FilterThreshold, byNotch.GetFilterTypeString()), Is.EqualTo((FilterType.QValueNotch, 0.01, "q-value notch")));
            Assert.That(byNotch.GetFilterValue(psms[0]), Is.EqualTo(0.001), "the value handed to quantification is the notch q-value");
            Assert.That(byNotch.TargetPsmsAboveThreshold, Is.EqualTo(120));
        }

        [Test]
        public static void QValueNotch_ByName_HasNoCountRule()
        {
            // The fewer-than-100 rule belongs to the inferred PEP filter; a filter asked for by name always filters.
            var filtered = FilteredPsms.Filter(Matches(50, 20, pepComputed: false), GuiPepBox, filterType: FilterType.QValueNotch);

            Assert.That(filtered.FilteringNotPerformed, Is.False);
            Assert.That(filtered.FilterThreshold, Is.EqualTo(1.0), "a filter asked for by name is held to its own threshold");
            Assert.That(filtered.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void PepQValue_ByName_WithoutPepQValues_FiltersByQValueAndSaysSo()
        {
            // Not a throw: a search would run to the end and then fail before writing anything.
            var fellBack = FilteredPsms.Filter(Matches(300, 200, pepComputed: false), GuiPepBox, filterType: FilterType.PepQValue);
            Assert.That((fellBack.FilterType, fellBack.FilterThreshold, fellBack.FilteredPsmsList.Count), Is.EqualTo((FilterType.QValue, 0.01, 200)));
            Assert.That(fellBack.PepFilteringNotDone, Is.EqualTo("PEP filtering was requested, but PEP was not computed for these PSMs, so results were filtered by q-value."));

            // With PEP q-values it filters on them, at the PEP q-value threshold.
            var filtered = FilteredPsms.Filter(Matches(300, 150, pepComputed: true), GuiPepBox, filterType: FilterType.PepQValue);
            Assert.That((filtered.FilterType, filtered.FilterThreshold, filtered.FilteredPsmsList.Count), Is.EqualTo((FilterType.PepQValue, 0.01, 150)));
        }

        /// <summary>
        /// The whole search, with the GUI's PEP box ticked: the two TaGe snips give between 100 and 1000 matches, so
        /// FdrAnalysisEngine never trains PEP. On master every match, peptide and protein group was dropped. Now the
        /// matches are filtered by q-value at 0.01, results.txt and the warning say so, and protein groups are held
        /// to the same 0.01, not to the q-value threshold of 1.
        /// </summary>
        [Test]
        public static void SearchTask_PepNotComputed_WritesQValueFilteredResultsAndSaysSo()
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDirectory, "PepFilteringNotDone_SearchTask");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var task = new SearchTask { CommonParameters = GuiPepBox };
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
            const string said = "PEP filtering was requested, but PEP was not computed for these PSMs, so results were filtered by q-value.";
            const string saidForPeptides = "PEP filtering was requested, but PEP q-values were not computed at the peptide level, so results were filtered by q-value.";
            // On master the three counts were "pep q-value <= 0.01: 0" for PSMs, peptides and protein groups.
            Assert.That(results.Skip(5).Take(5), Is.EqualTo(new[]
            {
                said, // just above the PSM count it qualifies
                "All target PSMs with q-value <= 0.01: 258",
                saidForPeptides, // just above the peptide count
                "All target peptides with q-value <= 0.01: 88",
                "All target protein groups with q-value <= 0.01 (1% FDR): 84",
            }));
            Assert.That(warnings.Where(w => w.StartsWith("PEP filtering was requested")), Is.EqualTo(new[] { said }), "one warning per search");

            // Protein groups: held to the matches' 0.01, not to the q-value threshold of 1.
            string[] groups = File.ReadAllLines(Directory.GetFiles(searchFolder, "All*ProteinGroups.tsv").Single());
            int qColumn = groups[0].Split('\t').ToList().IndexOf("Protein QValue");
            Assert.That(qColumn, Is.Not.EqualTo(-1));
            var groupQValues = groups.Skip(1).Select(l => double.Parse(l.Split('\t')[qColumn], System.Globalization.CultureInfo.InvariantCulture)).ToList();
            Assert.That(groupQValues, Has.Count.EqualTo(84));
            Assert.That(groupQValues, Has.All.LessThanOrEqualTo(0.01));

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void GetFilterValue_FollowsTheFilterType()
        {
            var psm = Matches(1, 0, pepComputed: false)[0];
            psm.SetFdrValues(1, 0, 0.3, 1, 0, 0.2, 0.05, 0.1);

            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.QValue).GetFilterValue(psm), Is.EqualTo(0.3));
            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.QValueNotch).GetFilterValue(psm), Is.EqualTo(0.2));
            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.PepQValue).GetFilterValue(psm), Is.EqualTo(0.1));
        }

        [Test]
        public static void AFilterTypeWithNoCase_Throws()
        {
            // A FilterType added later without its case here must fail loudly, not be read as the q-value.
            var psm = Matches(1, 0, pepComputed: false)[0];
            var unknown = FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: (FilterType)99);

            Assert.Throws<ArgumentOutOfRangeException>(() => unknown.GetFilterValue(psm));
            Assert.Throws<ArgumentOutOfRangeException>(() => unknown.GetFilterTypeString());
        }
    }
}
