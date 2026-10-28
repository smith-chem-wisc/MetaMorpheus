using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.ModificationAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// <see cref="CommonParameters.IdentificationFilter"/>: a search names the FDR value that decides identification,
    /// instead of MetaMorpheus inferring it from the thresholds. The FDR engine records, on each match and at each
    /// level, whether PEP was trained for the match's digestion agent (<see cref="FdrInfo.PepTrained"/>), and every
    /// consumer judges each match by that record, strictly below the threshold.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // the whole-search tests set FdrAnalysisEngine.QvalueThresholdOverride, a process-wide static
    public static class IdentificationFilterSettingTests
    {
        private static CommonParameters Named(IdentificationFilterMode mode, double qValueThreshold = 1.0, double pepQValueThreshold = 0.01) =>
            new CommonParameters(qValueThreshold: qValueThreshold, pepQValueThreshold: pepQValueThreshold, identificationFilter: mode);

        /// <summary>
        /// <paramref name="count"/> target matches to <paramref name="sequence"/>, digested with <paramref name="protease"/>;
        /// the first <paramref name="confident"/> have q, notch and PEP q-value 0.001, the rest 0.5. When PEP was not
        /// trained, every PEP q-value is FdrInfo's initial 2, as FdrAnalysisEngine leaves it.
        /// </summary>
        private static List<SpectralMatch> Matches(string sequence, string protease, int count, int confident, bool pepTrained, double pep = 0)
        {
            var peptide = new Protein(sequence, "acc_" + sequence)
                .Digest(new DigestionParams(protease: protease, minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First();
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());

            // a match's DigestionParams, and so its digestion agent, come from the CommonParameters it is made with
            var searchParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease, minPeptideLength: 1));
            var psms = new List<SpectralMatch>();
            for (int i = 0; i < count; i++)
            {
                var psm = new PeptideSpectralMatch(peptide, 0, 10, i, scan, searchParameters, new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                double q = i < confident ? 0.001 : 0.5;
                psm.SetFdrValues(i + 1, 0, q, i + 1, 0, q, pep: pep, pepQValue: pepTrained ? q : 2);
                psm.PsmFdrInfo.PepTrained = pepTrained;
                psms.Add(psm);
            }
            return psms;
        }

        // Trypsin trained PEP; Arg-C, with too few matches, did not.
        private static List<SpectralMatch> TwoProteases() =>
            Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true)
                .Concat(Matches("PEPTIDER", "Arg-C", 200, 100, pepTrained: false))
                .ToList();

        [Test]
        public static void DefaultIsThresholds_AndTheTomlKeepsANamedFilter()
        {
            Assert.That(new CommonParameters().IdentificationFilter, Is.EqualTo(IdentificationFilterMode.Thresholds));

            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, "IdentificationFilterSetting.toml");
            var task = new SearchTask { CommonParameters = Named(IdentificationFilterMode.PepQValueElseQValueNotch) };
            Nett.Toml.WriteFile(task, toml, MetaMorpheusTask.tomlConfig);
            var read = Nett.Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig);
            File.Delete(toml);

            Assert.That(read.CommonParameters.IdentificationFilter, Is.EqualTo(IdentificationFilterMode.PepQValueElseQValueNotch));
            Assert.That(read.CommonParameters.CloneWithNewTerminus().IdentificationFilter, Is.EqualTo(IdentificationFilterMode.PepQValueElseQValueNotch));
        }

        [Test]
        public static void PepElseNotch_AProteaseWithoutPep_IsJudgedOnTheNotch_NotOnPepQValue2()
        {
            var filtered = FilteredPsms.Filter(TwoProteases(), Named(IdentificationFilterMode.PepQValueElseQValueNotch));

            // On a PEP q-value of 2, every Arg-C match would be dropped.
            Assert.That(filtered.FilteredPsmsList.Count(p => p.DigestionParams.DigestionAgent.Name == "trypsin"), Is.EqualTo(150));
            Assert.That(filtered.FilteredPsmsList.Count(p => p.DigestionParams.DigestionAgent.Name == "Arg-C"), Is.EqualTo(100));
            Assert.That(filtered.GetFilterTypeString(), Is.EqualTo("pep q-value or q-value notch"));
            Assert.That(filtered.GetThresholdString(), Is.EqualTo("< 0.01"));
            Assert.That(IdentificationRule.Describe(TwoProteases(), peptideLevel: false, Named(IdentificationFilterMode.PepQValueElseQValueNotch)),
                Is.EqualTo("q-value notch < 0.01 (PEP not trained: 200 PSMs) for Arg-C; pep q-value < 0.01 for trypsin"));
            Assert.That(IdentificationRule.FellBack(TwoProteases(), peptideLevel: false, Named(IdentificationFilterMode.PepQValueElseQValueNotch)), Is.True);
        }

        [Test]
        public static void PepQValue_AProteaseWithoutPep_FailsNamingIt()
        {
            var ex = Assert.Throws<MetaMorpheusException>(() =>
                FilteredPsms.Filter(TwoProteases(), Named(IdentificationFilterMode.PepQValue)));
            Assert.That(ex!.Message, Does.Contain("PEP was not trained at the PSM level for Arg-C (200 PSMs)"));
            Assert.That(ex.Message, Does.Not.Contain("trypsin"));

            // Every protease trained: it filters on the PEP q-value.
            var filtered = FilteredPsms.Filter(Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true), Named(IdentificationFilterMode.PepQValue));
            Assert.That((filtered.FilterType, filtered.FilterThreshold, filtered.FilteredPsmsList.Count), Is.EqualTo((FilterType.PepQValue, 0.01, 150)));
        }

        [Test]
        public static void FilteringTheFilteredMatchesAgain_GivesTheSameAnswer()
        {
            // pcruzparri's case: every confident target has PEP exactly 0. Judged from the values alone, the second pass
            // would read that as a failed training and switch filters. The record on each match does not change.
            var once = FilteredPsms.Filter(Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true, pep: 0),
                Named(IdentificationFilterMode.PepQValueElseQValueNotch));
            var twice = FilteredPsms.Filter(once, Named(IdentificationFilterMode.PepQValueElseQValueNotch));

            Assert.That((twice.FilterType, twice.FilteredPsmsList.Count), Is.EqualTo((once.FilterType, once.FilteredPsmsList.Count)));
            Assert.That(once.FilterType, Is.EqualTo(FilterType.PepQValue));
        }

        [Test]
        public static void CountPsm_AgreesWithTheFilter()
        {
            // The Spectrum Match Count column comes from CountPsm, per protease, inside the FDR engine.
            var matches = TwoProteases();
            var commonParameters = Named(IdentificationFilterMode.PepQValueElseQValueNotch);
            foreach (var protease in matches.GroupBy(m => m.DigestionParams.DigestionAgent))
                FdrAnalysisEngine.CountPsm(protease.ToList(), commonParameters);

            var kept = FilteredPsms.Filter(matches, commonParameters).FilteredPsmsList;
            Assert.That(matches.First(m => m.BaseSequence == "PEPTIDER").PsmCount, Is.EqualTo(kept.Count(m => m.BaseSequence == "PEPTIDER")).And.EqualTo(100));
            Assert.That(matches.First(m => m.BaseSequence == "PEPTIDEK").PsmCount, Is.EqualTo(kept.Count(m => m.BaseSequence == "PEPTIDEK")).And.EqualTo(150));
        }

        [Test]
        public static void ANamedFilter_IsStrict_TheThresholdsFilterIsNot()
        {
            var atThreshold = Matches("PEPTIDEK", "trypsin", 1, 0, pepTrained: false);
            atThreshold[0].SetFdrValues(1, 0, 0.01, 1, 0, 0.01, 0, 2);

            Assert.That(FilteredPsms.Filter(atThreshold, new CommonParameters()).FilteredPsmsList, Has.Count.EqualTo(1), "q-value and notch <= 0.01, unchanged");
            Assert.That(FilteredPsms.Filter(atThreshold, Named(IdentificationFilterMode.QValueNotch, qValueThreshold: 0.01)).FilteredPsmsList, Is.Empty,
                "the named filter is mzLib's QuantifiedPsmRule: strictly below");
        }

        [Test]
        public static void PsmsWithin1PercentFdr_FollowsTheNamedFilter()
        {
            // 120 trypsin matches fail the q-value but pass the notch; Arg-C has no PEP and 100 confident matches.
            var matches = TwoProteases();
            foreach (var psm in matches.Where(m => m.BaseSequence == "PEPTIDEK").Take(120))
            {
                psm.SetFdrValues(1, 0, 0.5, 1, 0, 0.001, 0, 0.5);
                psm.PsmFdrInfo.PepTrained = true; // SetFdrValues makes a new FdrInfo; trypsin's PEP was trained
            }

            Assert.That(FdrAnalysisEngine.CountWithin1PercentFdr(matches, new CommonParameters()), Is.EqualTo(30 + 100),
                "q-value <= 0.01, as before: the 120 with a q-value of 0.5 are out");
            Assert.That(FdrAnalysisEngine.CountWithin1PercentFdr(matches, Named(IdentificationFilterMode.QValueNotch, qValueThreshold: 0.01)), Is.EqualTo(150 + 100),
                "the notch: the 120 are in");
            Assert.That(FdrAnalysisEngine.CountWithin1PercentFdr(matches, Named(IdentificationFilterMode.PepQValueElseQValueNotch)), Is.EqualTo(30 + 100),
                "trypsin on its PEP q-value (the 120 have 0.5), Arg-C on the notch");
        }

        [Test]
        public static void ProteinGroupAtExactlyTheThreshold_IsWrittenOnlyUnderTheThresholds()
        {
            // AllProteinGroups.tsv and quantification share PostSearchAnalysisTask.ProteinGroupIsWritten.
            var proteinGroup = new ProteinGroup(new HashSet<Omics.IBioPolymer> { new Protein("MEDEEK", "prot1") },
                new HashSet<Omics.IBioPolymerWithSetMods>(), new HashSet<Omics.IBioPolymerWithSetMods>())
            {
                QValue = 0.01
            };
            var isWritten = typeof(PostSearchAnalysisTask).GetMethod("ProteinGroupIsWritten",
                System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)!;

            bool Written(CommonParameters commonParameters)
            {
                var task = new PostSearchAnalysisTask
                {
                    CommonParameters = commonParameters,
                    Parameters = new PostSearchAnalysisParameters { SearchParameters = new SearchParameters { WriteHighQValuePsms = false } }
                };
                return (bool)isWritten.Invoke(task, new object[] { proteinGroup, 0.01 })!;
            }

            Assert.That(Written(new CommonParameters()), Is.True, "q-value <= 0.01, as before");
            Assert.That(Written(Named(IdentificationFilterMode.QValueNotch, qValueThreshold: 0.01)), Is.False, "a named filter is strict");
        }

        [Test]
        public static void ModificationCounts_FollowTheNamedFilter()
        {
            // A match that fails the q-value but passes the notch: counted under QValueNotch, not under the thresholds.
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif);
            var mod = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            var protein = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: new Dictionary<int, List<Modification>> { { 2, new List<Modification> { mod } } });
            var peptide = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0,
                new Dictionary<int, Modification> { { 2, mod } }, 0);
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, 0.5, 1, 0, 0.001, 0, 2);
            var fsp = new List<(string, CommonParameters)> { ("", new CommonParameters()) };

            int Counted(CommonParameters commonParameters)
            {
                var results = (ModificationAnalysisResults)new ModificationAnalysisEngine(new List<SpectralMatch> { psm }, commonParameters, fsp, new List<string>()).Run();
                return results.CountOfEachModSeenOnProteins.Values.Sum() + results.CountOfModsSeenAndLocalized.Values.Sum();
            }

            Assert.That(Counted(new CommonParameters()), Is.EqualTo(0));
            Assert.That(Counted(Named(IdentificationFilterMode.QValueNotch, qValueThreshold: 0.01)), Is.GreaterThan(0));
        }

        /// <summary>
        /// A whole SearchTask on the two TaGe snips: 258 target PSMs and no decoys. PEP is not trained, and with the
        /// override (<paramref name="trainPep"/>) training is attempted and fails.
        /// </summary>
        private static (string Folder, string[] Results, List<string> Warnings) Search(string name, CommonParameters commonParameters, bool trainPep,
            bool writeIdentificationFiles = false)
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDirectory, "IdentificationFilterSetting_" + name);
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var task = new SearchTask
            {
                CommonParameters = commonParameters,
                SearchParameters = new SearchParameters { WriteMzId = writeIdentificationFiles, WritePepXml = writeIdentificationFiles }
            };
            var spectraFiles = new List<string>
            {
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML"),
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML")
            };
            var databases = new List<DbForTask> { new DbForTask(Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta"), false) };
            var warnings = new List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);

            var qValueOverride = typeof(FdrAnalysisEngine).GetProperty("QvalueThresholdOverride")!;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                MetaMorpheusTask.WarnHandler += onWarn;
                try
                {
                    qValueOverride.SetValue(null, trainPep);
                    Directory.CreateDirectory(outputFolder);
                    task.RunTask(outputFolder, databases, spectraFiles, name);
                }
                finally
                {
                    qValueOverride.SetValue(null, false);
                    MetaMorpheusTask.WarnHandler -= onWarn;
                }
            }
            return (outputFolder, File.ReadAllLines(Path.Combine(outputFolder, "results.txt")), warnings);
        }

        [Test]
        public static void SearchTask_PepElseNotch_PepNotTrained_FallsBackToTheNotchAndSaysSo()
        {
            var (folder, results, warnings) = Search("fallback", Named(IdentificationFilterMode.PepQValueElseQValueNotch), trainPep: false);

            string psmLine = results.Single(l => l.StartsWith("Identification filter for PSMs: "));
            Assert.That(psmLine, Does.Match(@"^Identification filter for PSMs: q-value notch < 0\.01 \(PEP not trained: \d+ PSMs\)$"));
            Assert.That(results, Has.One.Match(@"^Identification filter for peptides: q-value notch < 0\.01 \(PEP not trained: \d+ peptides\)$"));
            Assert.That(warnings, Has.Some.EqualTo(psmLine));
            Assert.That(int.Parse(results.Single(l => l.StartsWith("All target PSMs with q-value notch < 0.01: ")).Split(": ")[1]), Is.GreaterThan(0));
            Assert.That(results, Has.One.StartsWith("All target protein groups with q-value notch < 0.01: "));
            Assert.That(File.ReadAllText(Path.Combine(folder, "AutoGeneratedManuscriptProse.txt")),
                Does.Contain("Identifications were filtered on one value throughout (PepQValueElseQValueNotch): PSMs at q-value notch < 0.01"));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The HeLa subset searched as PepIterativeTrainingTests searches it: with the override, PEP trains on it.
        /// </summary>
        [Test]
        public static void FdrEngine_RecordsWhetherItTrainedPep_AndALaterPassKeepsTheRecord()
        {
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams());
            var fsp = new List<(string, CommonParameters)> { ("TaGe_SA_HeLa_04_subset_longestSeq.mzML", commonParameters) };
            string dataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            var msDataFile = new MyFileManager(true).LoadFile(dataFile, commonParameters);
            List<Protein> proteins = UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinFasta(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, UsefulProteomicsDatabases.DecoyType.Reverse, false, out _,
                UsefulProteomicsDatabases.ProteinDbLoader.UniprotAccessionRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotFullNameRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotFullNameRegex,
                UsefulProteomicsDatabases.ProteinDbLoader.UniprotGeneNameRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotOrganismRegex, -1);
            var scans = MetaMorpheusTask.GetMs2Scans(msDataFile, dataFile, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] psmArray = new PeptideSpectralMatch[scans.Length];
            new EngineLayer.ClassicSearch.ClassicSearchEngine(psmArray, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteins, new SinglePpmAroundZeroSearchMode(5), commonParameters, fsp, null, new List<string>(), false).Run();
            var psms = psmArray.Where(p => p != null).ToList();

            var qValueOverride = typeof(FdrAnalysisEngine).GetProperty("QvalueThresholdOverride")!;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                try
                {
                    qValueOverride.SetValue(null, false);
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();
                    Assert.That(psms.Count, Is.LessThanOrEqualTo(1000));
                    Assert.That(psms.Select(p => (p.PsmFdrInfo.PepTrained, p.PeptideFdrInfo.PepTrained)).Distinct(), Is.EqualTo(new[] { (false, false) }),
                        "1000 matches or fewer: PEP is not trained");

                    // The override's (decoys + 1) / targets q-values leave no positive example below the default 0.005.
                    qValueOverride.SetValue(null, true);
                    commonParameters.QValueCutoffForPepCalculation = 0.1;
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();
                    Assert.That(psms.Select(p => p.PsmFdrInfo.PEP).Distinct().Count(), Is.GreaterThan(1), "a trained model, not PEPs left at 0");
                    Assert.That(psms.Select(p => (p.PsmFdrInfo.PepTrained, p.PeptideFdrInfo.PepTrained)).Distinct(), Is.EqualTo(new[] { (true, true) }));

                    // A later pass that reuses the PEPs (disambiguation, the variant tables) trains nothing and keeps the record.
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>(), doPEP: false).Run();
                    Assert.That(psms.All(p => p.PsmFdrInfo.PepTrained && p.PeptideFdrInfo.PepTrained));
                }
                finally
                {
                    qValueOverride.SetValue(null, false);
                }
            }

            var named = Named(IdentificationFilterMode.PepQValueElseQValueNotch, pepQValueThreshold: 0.5);
            Assert.That(FilteredPsms.Filter(psms, named).FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(IdentificationRule.Describe(psms, peptideLevel: false, named), Is.EqualTo("pep q-value < 0.5"));
            Assert.That(IdentificationRule.FellBack(psms, peptideLevel: false, named), Is.False);
        }

        [Test]
        public static void SearchTask_PepElseNotch_PepTrainingFailed_IsNotTreatedAsTrained()
        {
            // The TaGe snips have no decoy to train on: PepAnalysisEngine gives up, every PEP stays 0, and the FDR engine
            // still ranks them into PEP q-values. Those are not a trained PEP's, and the record says so.
            var (folder, results, _) = Search("trainingFailed", Named(IdentificationFilterMode.PepQValueElseQValueNotch), trainPep: true);

            Assert.That(results, Has.One.Match(@"^Identification filter for PSMs: q-value notch < 0\.01 \(PEP not trained: 258 PSMs\)$"));

            Directory.Delete(folder, true);
        }

        [Test]
        public static void SearchTask_PepQValue_PepNotTrained_FailsTheTask()
        {
            var ex = Assert.Throws<MetaMorpheusException>(() => Search("fails", Named(IdentificationFilterMode.PepQValue), trainPep: false));
            Assert.That(ex!.Message, Does.StartWith("The identification filter is the PEP q-value, but PEP was not trained at the PSM level for trypsin ("));

            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "IdentificationFilterSetting_fails"), true);
        }

        [Test]
        public static void SearchTask_NamedFilter_ReachesMzIdentMlAndPepXml()
        {
            // At 0.05, so a file that still said the hard-coded 0.01 would fail.
            var (folder, _, _) = Search("writers", Named(IdentificationFilterMode.QValueNotch, qValueThreshold: 0.05, pepQValueThreshold: 1.0),
                trainPep: false, writeIdentificationFiles: true);
            string individual = Path.Combine(folder, "Individual File Results");

            XDocument mzId = XDocument.Load(Path.Combine(individual, "TaGe_SA_A549_3_snip.mzID"));
            var items = mzId.Descendants().Where(e => e.Name.LocalName == "SpectrumIdentificationItem").ToList();
            Assert.That(items, Is.Not.Empty);
            Assert.That(items.Select(e => (string)e.Attribute("passThreshold")), Has.All.EqualTo("true"),
                "every written match passed the named filter, including those with a q-value notch between 0.01 and 0.05");
            var thresholds = mzId.Descendants().Where(e => e.Name.LocalName == "cvParam" && ((string)e.Attribute("name"))?.EndsWith("FDR threshold") == true).ToList();
            Assert.That(thresholds.Select(e => (string)e.Attribute("value")), Is.EqualTo(new[] { "0.05", "0.05" }));

            string pepXml = File.ReadAllText(Path.Combine(individual, "TaGe_SA_A549_3_snip.pep.XML"));
            Assert.That(pepXml, Does.Contain("name=\"Identification Filter\" value=\"q-value notch &lt; 0.05\""));

            Directory.Delete(folder, true);
        }
    }
}
