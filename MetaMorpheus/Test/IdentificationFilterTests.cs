using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using GuiFunctions.Util;
using EngineLayer.ModificationAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
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
    /// <see cref="CommonParameters.IdentificationFilter"/>: the user picks the value that decides identification, the PEP
    /// q-value or the q-value notch, and every step of a search uses it at that filter's own threshold. Where PEP was not
    /// trained for a match, that match is judged by the q-value notch instead, and the search says so.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // the whole-search tests set FdrAnalysisEngine.QvalueThresholdOverride, a process-wide static
    public static class IdentificationFilterTests
    {
        private static CommonParameters Chosen(FilterType filter, double qValueThreshold = 0.01, double pepQValueThreshold = 0.01) =>
            new CommonParameters(qValueThreshold: qValueThreshold, pepQValueThreshold: pepQValueThreshold, identificationFilter: filter);

        /// <summary>
        /// <paramref name="count"/> target matches to <paramref name="sequence"/>, digested with <paramref name="protease"/>.
        /// The first <paramref name="confident"/> have q-value, notch and PEP q-value 0.001, the rest 0.5. When PEP was not
        /// trained, every PEP q-value is FdrInfo's initial 2, as FdrAnalysisEngine leaves it.
        /// </summary>
        private static List<SpectralMatch> Matches(string sequence, string protease, int count, int confident, bool pepTrained,
            CleavageSpecificity specificity = CleavageSpecificity.Full)
        {
            var peptide = new Protein(sequence, "acc_" + sequence)
                .Digest(new DigestionParams(protease: protease, minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First();
            peptide.CleavageSpecificityForFdrCategory = specificity;
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());

            // a match's digestion agent comes from the CommonParameters it is made with
            var searchParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease, minPeptideLength: 1));
            var psms = new List<SpectralMatch>();
            for (int i = 0; i < count; i++)
            {
                var psm = new PeptideSpectralMatch(peptide, 0, 10, i, scan, searchParameters, new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                double q = i < confident ? 0.001 : 0.5;
                psm.SetFdrValues(i + 1, 0, q, i + 1, 0, q, pep: 0, pepQValue: pepTrained ? q : 2);
                psm.PsmFdrInfo.PepTrained = pepTrained;
                psms.Add(psm);
            }
            return psms;
        }

        private static void Set(SpectralMatch psm, double q, double notch, double pepQ, bool pepTrained = true, double pep = 0)
        {
            psm.SetFdrValues(1, 0, q, 1, 0, notch, pep, pepQ);
            psm.PsmFdrInfo.PepTrained = pepTrained;
        }

        private static int CountFor(FilteredPsms filtered, string protease) =>
            filtered.FilteredPsmsList.Count(p => p.DigestionParams.DigestionAgent.Name == protease);

        // Spec 1: the user picks the filter; the thresholds never choose it.
        [Test]
        public static void TheSearchNamesItsFilter_TheThresholdsDoNotChooseIt()
        {
            var trained = Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true);

            // Equal thresholds meant the q-value on master.
            var pep = FilteredPsms.Filter(trained, Chosen(FilterType.PepQValue, qValueThreshold: 0.01, pepQValueThreshold: 0.01));
            Assert.That(pep.FilterType, Is.EqualTo(FilterType.PepQValue));

            // A PEP threshold below the q-value threshold meant PEP on master.
            var notch = FilteredPsms.Filter(trained, Chosen(FilterType.QValueNotch, qValueThreshold: 1.0, pepQValueThreshold: 0.01));
            Assert.That((notch.FilterType, notch.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 1.0)));
        }

        // Spec 1b: the default is the q-value notch alone, not the q-value and the notch together.
        [Test]
        public static void TheDefaultFilter_IsTheQValueNotchAlone()
        {
            Assert.That(new CommonParameters().IdentificationFilter, Is.EqualTo(FilterType.QValueNotch));

            var psms = Matches("PEPTIDEK", "trypsin", 300, 0, pepTrained: false);
            foreach (var psm in psms.Take(120)) Set(psm, q: 0.5, notch: 0.001, pepQ: 2, pepTrained: false);

            var filtered = FilteredPsms.Filter(psms, new CommonParameters());
            Assert.That((filtered.FilterType, filtered.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(120), "the notch alone decides");
        }

        // Spec 1: saved in the TOML; a task that names no filter gets the default; the pooled q-value is not a choice.
        [Test]
        public static void TheChoiceIsSavedInTheToml_ATaskThatNamesNoneGetsTheNotch_AndTheQValueIsNotAChoice()
        {
            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, "IdentificationFilterTests.toml");
            try
            {
                Nett.Toml.WriteFile(new SearchTask { CommonParameters = Chosen(FilterType.PepQValue) }, toml, MetaMorpheusTask.tomlConfig);
                string[] lines = File.ReadAllLines(toml);
                Assert.That(lines, Has.One.Match(@"^\s*IdentificationFilter = ""PepQValue""$"));
                Assert.That(Nett.Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig).CommonParameters.IdentificationFilter,
                    Is.EqualTo(FilterType.PepQValue));

                File.WriteAllLines(toml, lines.Where(l => !l.TrimStart().StartsWith("IdentificationFilter")));
                Assert.That(Nett.Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig).CommonParameters.IdentificationFilter,
                    Is.EqualTo(FilterType.QValueNotch), "a task saved before the filter was named");
            }
            finally
            {
                File.Delete(toml);
            }

            Assert.Throws<ArgumentException>(() => new CommonParameters(identificationFilter: FilterType.QValue));
        }

        // Spec 2: each filter is held to its own threshold, and a PEP fallback to the notch at the PEP threshold.
        [Test]
        public static void EachFilterIsHeldToItsOwnThreshold_AndAPepFallbackUsesThePepThreshold()
        {
            var psm = Matches("PEPTIDEK", "trypsin", 1, 0, pepTrained: true)[0];

            Set(psm, q: 0.5, notch: 0.03, pepQ: 0.5);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.QValueNotch, qValueThreshold: 0.05, pepQValueThreshold: 0.01)).FilteredPsmsList,
                Has.Count.EqualTo(1), "the notch at the q-value threshold");

            Set(psm, q: 0.5, notch: 0.5, pepQ: 0.03);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.PepQValue, qValueThreshold: 0.01, pepQValueThreshold: 0.05)).FilteredPsmsList,
                Has.Count.EqualTo(1), "the PEP q-value at the PEP threshold");

            Set(psm, q: 0.5, notch: 0.03, pepQ: 2, pepTrained: false);
            var fellBack = FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.PepQValue, qValueThreshold: 1.0, pepQValueThreshold: 0.05));
            Assert.That(fellBack.FilteredPsmsList, Has.Count.EqualTo(1), "untrained: the notch, at the PEP threshold");
            Assert.That(fellBack.GetFilterValue(psm), Is.EqualTo(0.03), "the value handed to quantification is the one it was judged by");

            Set(psm, q: 0.001, notch: 0.07, pepQ: 2, pepTrained: false);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.PepQValue, qValueThreshold: 1.0, pepQValueThreshold: 0.05)).FilteredPsmsList,
                Is.Empty, "the fallback is the notch alone, not the q-value");
        }

        // Spec 5 (README C): one protease trained PEP, the other did not.
        [Test]
        public static void PepNotTrainedForOneProtease_ItsMatchesAreJudgedByTheNotch_AndTheSentenceNamesIt()
        {
            var both = Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true)
                .Concat(Matches("PEPTIDER", "Arg-C", 200, 100, pepTrained: false)).ToList();

            var filtered = FilteredPsms.Filter(both, Chosen(FilterType.PepQValue));

            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(CountFor(filtered, "trypsin"), Is.EqualTo(150));
            Assert.That(CountFor(filtered, "Arg-C"), Is.EqualTo(100), "judged by the notch, not dropped on a PEP q-value of 2");
            Assert.That(filtered.PepFilteringNotDone, Is.EqualTo(
                "PEP filtering was requested, but PEP was not trained for 200 target PSMs (Arg-C: 200), so those were filtered by q-value notch <= 0.01."));
        }

        // Spec 5: a semi- or non-specific search runs FDR, and trains PEP, once per cleavage-specificity category.
        [Test]
        public static void PepNotTrainedForOneCategory_IsReportedForThatCategory()
        {
            var both = Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true, specificity: CleavageSpecificity.Full)
                .Concat(Matches("PEPTIDER", "trypsin", 200, 100, pepTrained: false, specificity: CleavageSpecificity.Semi)).ToList();

            var filtered = FilteredPsms.Filter(both, Chosen(FilterType.PepQValue));

            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(250));
            Assert.That(filtered.PepFilteringNotDone, Does.Contain("(trypsin (Semi): 200)"));
        }

        // Spec 5 (README E): fewer than 100 matches are filtered, not kept whole.
        [Test]
        public static void FewerThan100Matches_AreFiltered()
        {
            var filtered = FilteredPsms.Filter(Matches("PEPTIDER", "Arg-C", 50, 20, pepTrained: false), Chosen(FilterType.PepQValue, qValueThreshold: 1.0));

            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(20));
            Assert.That(filtered.PepFilteringNotDone, Does.StartWith("PEP filtering was requested, but PEP was not trained for 50 target PSMs"));
        }

        // Spec 3 (README F): the Spectrum Match Count column follows the chosen filter.
        [Test]
        public static void SpectrumMatchCount_FollowsTheChosenFilter()
        {
            var psms = Matches("PEPTIDEK", "trypsin", 300, 0, pepTrained: true);
            for (int i = 0; i < 150; i++) Set(psms[i], q: i < 100 ? 0.001 : 0.5, notch: i < 100 ? 0.001 : 0.5, pepQ: 0.001);

            FdrAnalysisEngine.CountPsm(psms, Chosen(FilterType.PepQValue));

            Assert.That(psms.Select(p => p.PsmCount).Distinct(), Is.EqualTo(new[] { 150 }), "the 150 that pass the PEP q-value");
        }

        /// <summary>One protein, one peptide, one PSM at the given score, notch q-value and PEP, PEP trained.</summary>
        private static (ProteinGroup Group, SpectralMatch Psm) OneGroup(string accession, string sequence, bool decoy, double score, double notch, double pep)
        {
            var protein = new Protein(sequence, accession, isDecoy: decoy);
            var peptide = protein.Digest(new DigestionParams(protease: "trypsin", minPeptideLength: 1), new List<Modification>(), new List<Modification>()).First();
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());
            var psm = new PeptideSpectralMatch(peptide, 0, score, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            Set(psm, q: 0.001, notch: notch, pepQ: 0.001, pep: pep);
            var group = new ProteinGroup(new HashSet<IBioPolymer> { protein }, new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });
            return (group, psm);
        }

        private static Dictionary<string, double> ProteinQValues(CommonParameters commonParameters)
        {
            // X: best score 30, notch 0.0001, PEP 0.5. Y: best score 50, notch 0.005, PEP 0.0001. A decoy between them on
            // both. Ranked by the notch: X, decoy, Y. Ranked by PEP: Y, decoy, X.
            var all = new[]
            {
                OneGroup("X", "ACDEFK", decoy: false, score: 30, notch: 0.0001, pep: 0.5),
                OneGroup("Y", "GHILMK", decoy: false, score: 50, notch: 0.005, pep: 0.0001),
                OneGroup("DECOY_Z", "NQSTVK", decoy: true, score: 40, notch: 0.003, pep: 0.001),
            };
            var results = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(all.Select(a => a.Group).ToList(), all.Select(a => a.Psm).ToList(),
                noOneHitWonders: false, treatModPeptidesAsDifferentPeptides: true, mergeIndistinguishableProteinGroups: false, commonParameters,
                new List<(string, CommonParameters)>(), new List<string>()).Run();
            return results.SortedAndScoredProteinGroups.ToDictionary(g => g.ProteinGroupName, g => g.QValue);
        }

        // Spec 3: protein scoring and protein FDR use the chosen filter and stop inferring their own.
        [Test]
        public static void ProteinFdr_UsesTheChosenFilter_NotTheThresholds()
        {
            var byPep = ProteinQValues(Chosen(FilterType.PepQValue, qValueThreshold: 0.01, pepQValueThreshold: 0.01));
            Assert.That((byPep["Y"], byPep["X"]), Is.EqualTo((0.0, 0.5)), "ranked by PEP, though the thresholds are equal");

            var byNotch = ProteinQValues(Chosen(FilterType.QValueNotch, qValueThreshold: 1.0, pepQValueThreshold: 0.01));
            Assert.That((byNotch["X"], byNotch["Y"]), Is.EqualTo((0.0, 0.5)), "ranked by the notch, though the PEP threshold is lower");
        }

        /// <summary>A match to the first peptide of <paramref name="protein"/>, on its own scan.</summary>
        private static SpectralMatch Match(Protein protein, int notch, double score, int scanNumber)
        {
            var peptide = protein.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).First();
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                scanNumber, 2, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=" + scanNumber, double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());
            var psm = new PeptideSpectralMatch(peptide, notch, score, scanNumber, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            return psm;
        }

        // Spec 5b: the engines' counts in results.txt follow the chosen filter at its threshold and name it.
        [Test]
        public static void TheEngineCounts_FollowTheChosenFilter_AndNameIt()
        {
            // Through FdrAnalysisEngine.Run, which computes the q-values itself: 10 decoys in notch 0 outscore 100 targets
            // in notch 1, so every target fails the pooled q-value and passes the notch.
            var target = new Protein("PEPTIDEK", "target", databaseFilePath: "temp");
            var decoy = new Protein("EDITPEPK", "DECOY_target", isDecoy: true, databaseFilePath: "temp");
            var psms = Enumerable.Range(1, 10).Select(i => Match(decoy, notch: 0, score: 100, scanNumber: i))
                .Concat(Enumerable.Range(11, 100).Select(i => Match(target, notch: 1, score: 50, scanNumber: i)))
                .ToList();
            var commonParameters = new CommonParameters();
            var fdr = (FdrAnalysisResults)new FdrAnalysisEngine(psms, 2, commonParameters,
                new List<(string, CommonParameters)> { ("File", commonParameters) }, new List<string>(), doPEP: false).Run();
            Assert.That(psms.Where(p => !p.IsDecoy).All(p => p.FdrInfo.QValue > 0.01 && p.FdrInfo.QValueNotch <= 0.01));
            Assert.That(fdr.PsmsWithin1PercentFdr, Is.EqualTo(100));
            Assert.That(fdr.ToString(), Does.Contain("PSMs with q-value notch <= 0.01: 100"));

            // The localized-modification counts.
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif);
            var mod = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            var protein = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: new Dictionary<int, List<Modification>> { { 2, new List<Modification> { mod } } });
            var peptide = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0,
                new Dictionary<int, Modification> { { 2, mod } }, 0);
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());
            var modPsm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            modPsm.ResolveAllAmbiguities();
            Set(modPsm, q: 0.5, notch: 0.001, pepQ: 2, pepTrained: false);
            var mods = (ModificationAnalysisResults)new ModificationAnalysisEngine(new List<SpectralMatch> { modPsm }, commonParameters,
                new List<(string, CommonParameters)> { ("", commonParameters) }, new List<string>()).Run();
            Assert.That(mods.CountOfEachModSeenOnProteins.Values.Sum() + mods.CountOfModsSeenAndLocalized.Values.Sum(), Is.GreaterThan(0));
            Assert.That(mods.ToString(), Does.Contain("Localized mods seen with q-value notch <= 0.01:"));

            // The protein groups, at the chosen threshold.
            var groups = ProteinQValues(Chosen(FilterType.QValueNotch, qValueThreshold: 0.05));
            Assert.That(groups["X"], Is.EqualTo(0));
            var proteinResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(new List<ProteinGroup>(), new List<SpectralMatch>(),
                false, true, false, Chosen(FilterType.QValueNotch, qValueThreshold: 0.05), new List<(string, CommonParameters)>(), new List<string>()).Run();
            Assert.That(proteinResults.ToString(), Does.Contain("Number of proteins with q-value <= 0.05: 0"));
        }

        // Spec 5b: mzIdentML's pass flags follow the chosen filter, for matches and protein groups.
        [Test]
        public static void MzIdentMl_PassThreshold_FollowsTheChosenFilter()
        {
            var protein = new Protein("PEPTIDEK", "ACC1", databaseFilePath: "temp");
            var notchOnly = Match(protein, notch: 0, score: 10, scanNumber: 1);
            Set(notchOnly, q: 0.5, notch: 0.001, pepQ: 2, pepTrained: false);
            var qValueOnly = Match(protein, notch: 0, score: 10, scanNumber: 2);
            Set(qValueOnly, q: 0.001, notch: 0.5, pepQ: 2, pepTrained: false);
            var psms = new List<SpectralMatch> { notchOnly, qValueOnly };
            var group = new ProteinGroup(new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods>(psms.SelectMany(p => p.BestMatchingBioPolymersWithSetMods).Select(p => p.SpecificBioPolymer)),
                new HashSet<IBioPolymerWithSetMods>())
            {
                QValue = 0.03
            };
            group.SequenceCoverageFraction.Add(0.5);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "IdentificationFilterTests_passThreshold.mzID");
            try
            {
                var filter = FilteredPsms.Filter(psms, Chosen(FilterType.QValueNotch, qValueThreshold: 0.05), includeHighQValuePsms: true);
                MzIdentMLWriter.WriteMzIdentMl(psms, new List<ProteinGroup> { group }, new List<Modification>(), new List<Modification>(),
                    null, new List<DigestionAgent>(), new MzLibUtil.PpmTolerance(20), new MzLibUtil.PpmTolerance(20), 0, path, true,
                    identificationFilter: filter);
                var elements = XDocument.Load(path).Descendants().ToList();
                string[] Of(string name) => elements.Where(e => e.Name.LocalName == name).Select(e => (string)e.Attribute("passThreshold")).ToArray();
                Assert.That(Of("SpectrumIdentificationItem"), Is.EqualTo(new[] { "true", "false" }), "the notch decides each match");
                Assert.That(Of("ProteinDetectionHypothesis"), Is.EqualTo(new[] { "true" }), "a group at 0.03 passes a 0.05 filter");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        /// <summary>
        /// The FDR engine records, on each match and at each level, whether PEP was trained for it; a later pass that
        /// trains nothing keeps the record. The HeLa subset, searched as the PEP tests search it.
        /// </summary>
        [Test]
        public static void TheFdrEngineRecordsWhetherItTrainedPep_AndALaterPassKeepsTheRecord()
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

                    // A later pass that reuses the PEPs (the individual-file and variant tables) trains nothing and keeps the record.
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>(), doPEP: false).Run();
                    Assert.That(psms.All(p => p.PsmFdrInfo.PepTrained && p.PeptideFdrInfo.PepTrained));
                }
                finally
                {
                    qValueOverride.SetValue(null, false);
                }
            }
        }

        /// <summary>
        /// A whole SearchTask on the two TaGe snips: 258 target PSMs and no decoys. PEP is not trained; with
        /// <paramref name="trainPep"/> training is attempted and fails.
        /// </summary>
        private static (string Folder, string[] Results, List<string> Warnings) Search(string name, CommonParameters commonParameters, bool trainPep,
            bool writeIndividualFiles = false, bool writeIdentificationFiles = false)
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDirectory, "IdentificationFilterTests_" + name);
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var task = new SearchTask
            {
                CommonParameters = commonParameters,
                SearchParameters = new SearchParameters
                {
                    WriteIndividualFiles = writeIndividualFiles,
                    WriteMzId = writeIdentificationFiles,
                    WritePepXml = writeIdentificationFiles
                }
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

        // Spec 5 (README H, J): a failed training is not a trained PEP.
        [Test]
        public static void SearchTask_PepTrainingFailed_FallsBackToTheNotch_AndSaysSo()
        {
            var (folder, results, warnings) = Search("trainingFailed", Chosen(FilterType.PepQValue, qValueThreshold: 1.0), trainPep: true);
            foreach (string line in results.Where(l => l.StartsWith("PEP filtering") || l.StartsWith("All target") || l.StartsWith("Posterior")))
                TestContext.WriteLine("results.txt: " + (line.Length > 200 ? line.Substring(0, 200) + "..." : line));

            Assert.That(results, Has.Some.StartsWith("Posterior error probability analysis failed"));
            int psmLine = Array.FindIndex(results, l => l.StartsWith("All target PSMs with "));
            Assert.That(results[psmLine], Does.StartWith("All target PSMs with pep q-value <= 0.01: "));
            Assert.That(results[psmLine - 1], Is.EqualTo(
                "PEP filtering was requested, but PEP was not trained for 258 target PSMs (trypsin: 258), so those were filtered by q-value notch <= 0.01."));
            Assert.That(warnings, Has.One.EqualTo(results[psmLine - 1]), "one warning per search");

            Directory.Delete(folder, true);
        }

        // Spec 4 (README I): the per-file lines follow the same rule as the all-files lines.
        [Test]
        public static void SearchTask_PerFileLines_FollowTheSameRule_AndCarryTheSentence()
        {
            var (folder, results, warnings) = Search("perFile", Chosen(FilterType.PepQValue, qValueThreshold: 1.0), trainPep: false, writeIndividualFiles: true);
            foreach (string line in results.Where(l => l.StartsWith("PEP filtering") || l.Contains(" - Target ") || l.StartsWith("All target")))
                TestContext.WriteLine("results.txt: " + line);

            var countLines = Enumerable.Range(0, results.Length)
                .Where(i => results[i].Contains(" - Target PSMs with ") || results[i].Contains(" - Target peptides with ")
                    || results[i].StartsWith("All target PSMs with ") || results[i].StartsWith("All target peptides with "))
                .ToList();
            Assert.That(countLines, Has.Count.EqualTo(6), "all files and two files, PSMs and peptides");
            foreach (int i in countLines)
            {
                Assert.That(results[i], Does.Contain("pep q-value <= 0.01: "), "the same filter and threshold on every line");
                Assert.That(results[i - 1], Does.StartWith("PEP filtering was requested, but "), "the fallback sentence above every count it qualifies");
            }
            Assert.That(warnings.Count(w => w.StartsWith("PEP filtering was requested")), Is.EqualTo(1), "one warning per search");

            Directory.Delete(folder, true);
        }

        // Spec 5b (README I): mzIdentML and pepXML state the chosen filter and threshold.
        [Test]
        public static void SearchTask_TheChosenFilter_ReachesMzIdentMlAndPepXml()
        {
            // At 0.05, so a file that still said the hard-coded 0.01 or the q-value setting would fail.
            var (folder, _, _) = Search("writers", Chosen(FilterType.QValueNotch, qValueThreshold: 0.05, pepQValueThreshold: 1.0),
                trainPep: false, writeIdentificationFiles: true);
            string individual = Path.Combine(folder, "Individual File Results");

            XDocument mzId = XDocument.Load(Path.Combine(individual, "TaGe_SA_A549_3_snip.mzID"));
            var thresholds = mzId.Descendants().Where(e => e.Name.LocalName == "cvParam" && ((string)e.Attribute("name"))?.EndsWith("FDR threshold") == true).ToList();
            Assert.That(thresholds.Select(e => (string)e.Attribute("value")), Is.EqualTo(new[] { "0.05", "0.05" }));

            string pepXml = File.ReadAllText(Path.Combine(individual, "TaGe_SA_A549_3_snip.pep.XML"));
            Assert.That(pepXml, Does.Contain("name=\"Identification Filter\" value=\"q-value notch &lt;= 0.05\""));

            Directory.Delete(folder, true);
        }

        // Spec 8: a chosen filter at 1 is allowed, with a warning at the start and a line in results.txt.
        [Test]
        public static void SearchTask_AFilterAtOne_RunsAndWarnsThatItKeepsEveryMatch()
        {
            var (folder, results, warnings) = Search("atOne", Chosen(FilterType.QValueNotch, qValueThreshold: 1.0), trainPep: false);

            const string sentence = "This search filters by q-value notch <= 1, which keeps every match.";
            Assert.That(warnings, Has.One.EqualTo(sentence));
            Assert.That(results, Has.One.EqualTo(sentence));

            Directory.Delete(folder, true);
        }

        // Review of #2982 (nbollis): the setter's refusal must say what is refused and what is allowed.
        [Test]
        public static void ThePooledQValue_IsRefused_AndTheMessageNamesTheSettingAndBothChoices()
        {
            var ex = Assert.Throws<ArgumentException>(() => new CommonParameters(identificationFilter: FilterType.QValue));
            Assert.That(ex.Message, Does.Contain("IdentificationFilter").And.Contain("PepQValue").And.Contain("QValueNotch")
                .And.Contain("pooled q-value"));

            // A task TOML that names it does not load, and the message a user sees names the setting.
            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, "IdentificationFilterTests_qvalue.toml");
            try
            {
                Nett.Toml.WriteFile(new SearchTask { CommonParameters = Chosen(FilterType.QValueNotch) }, toml, MetaMorpheusTask.tomlConfig);
                File.WriteAllLines(toml, File.ReadAllLines(toml).Select(l => l.TrimStart().StartsWith("IdentificationFilter")
                    ? "IdentificationFilter = \"QValue\"" : l));
                var loadError = Assert.Catch(() => Nett.Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig));
                var messages = new List<string>();
                for (Exception e = loadError; e != null; e = e.InnerException) messages.Add(e.Message);
                Assert.That(messages, Has.Some.Contain("IdentificationFilter").And.Some.Contain("PepQValue"));
            }
            finally
            {
                File.Delete(toml);
            }
        }

        private static string TestData(string file) => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", file);

        private static CommonParameters LoadCommon(string tomlText)
        {
            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, $"IdentificationFilterTests_{Guid.NewGuid():N}.toml");
            try
            {
                File.WriteAllText(toml, tomlText);
                return Nett.Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig).CommonParameters;
            }
            finally
            {
                File.Delete(toml);
            }
        }

        private static string TomlOf(CommonParameters commonParameters, bool withoutTheFilter)
        {
            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, $"IdentificationFilterTests_{Guid.NewGuid():N}.toml");
            try
            {
                Nett.Toml.WriteFile(new SearchTask { CommonParameters = commonParameters }, toml, MetaMorpheusTask.tomlConfig);
                var lines = File.ReadAllLines(toml).Where(l => !withoutTheFilter || !l.TrimStart().StartsWith("IdentificationFilter"));
                return string.Join(Environment.NewLine, lines);
            }
            finally
            {
                File.Delete(toml);
            }
        }

        // Review of #2982 (nbollis), loading old TOMLs. A task saved with the old GUI's PEP box ticked names no filter and holds
        // QValueThreshold 1.0 and PepQValueThreshold 0.01 (TestData/Task2-SearchTaskconfig.toml, as on master). It loads as the PEP
        // q-value filter at 0.01, says so once, and is not reported as a filter at 1.
        [Test]
        public static void AnOldGuiPepTask_LoadsAsThePepFilter_AndSaysSo()
        {
            var common = Nett.Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig).CommonParameters;

            Assert.That((common.IdentificationFilter, common.PepQValueThreshold, common.QValueThreshold), Is.EqualTo((FilterType.PepQValue, 0.01, 1.0)));
            Assert.That(common.IdentificationFilterFromLegacyThresholds, Is.True);
            Assert.That(SearchTask.LegacyFilterSentence(common), Is.EqualTo(
                "IdentificationFilter read from this task's legacy thresholds (PepQValueThreshold 0.01 < QValueThreshold 1): PEP q-value at 0.01"));
            Assert.That(SearchTask.FilterKeepsEveryMatch(common), Is.Null, "a migrated PEP task filters at 0.01, not at 1");
            Assert.That(common.Clone().IdentificationFilterFromLegacyThresholds, Is.True, "a copy remembers where its filter came from");
        }

        // Review of #2982 (Alexander-Sol), his reproduction: an unmodified old PEP TOML (master's Task2-SearchTaskconfig.toml,
        // byte for byte) admitted ten trained matches with PEP q-value and notch 0.5, as a notch filter at its saved
        // QValueThreshold of 1.0. Read as the PEP filter at its saved 0.01, it rejects all ten.
        [Test]
        public static void AnUnmodifiedOldPepToml_RejectsTenTrainedMatchesAtOneHalf()
        {
            var common = Nett.Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig).CommonParameters;
            var psms = Matches("PEPTIDEK", "trypsin", 10, 0, pepTrained: true);
            Assume.That(psms.Select(p => (p.FdrInfo.PEP_QValue, p.FdrInfo.QValueNotch, p.FdrInfo.PepTrained)).Distinct(), Is.EqualTo(new[] { (0.5, 0.5, true) }));

            var filtered = FilteredPsms.Filter(psms, common);

            Assert.That((filtered.FilterType, filtered.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            Assert.That(filtered.FilteredPsmsList, Is.Empty, "none of the ten passes the saved 0.01");
            Assert.That(psms.Where(p => FilteredPsms.Passes(p, false, common.IdentificationFilter, FilteredPsms.Threshold(common, common.IdentificationFilter))), Is.Empty);
        }

        // Review of #2982 (Alexander-Sol): opening and saving an old PEP task in the Search task window keeps it a PEP task.
        // The window's load and save rule is TaskWindowIdentificationFilter.
        [Test]
        public static void TheSearchWindow_ShowsAnOldPepTaskTicked_AndSavesItsFilterByName()
        {
            var old = Nett.Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig).CommonParameters;

            bool ticked = TaskWindowIdentificationFilter.PepBoxTicked(old);
            string qText = TaskWindowIdentificationFilter.QValueThresholdText(old), pepText = TaskWindowIdentificationFilter.PepQValueThresholdText(old);
            Assert.That((ticked, qText, pepText), Is.EqualTo((true, "0.01", "0.01")), "the PEP box ticked, at the saved 0.01");

            var (filter, qValueThreshold, pepQValueThreshold) = TaskWindowIdentificationFilter.ForSave(ticked, qText, pepText);
            Assert.That((filter, qValueThreshold, pepQValueThreshold), Is.EqualTo((FilterType.PepQValue, 1.0, 0.01)));
            string saved = TomlOf(new CommonParameters(qValueThreshold: qValueThreshold, pepQValueThreshold: pepQValueThreshold, identificationFilter: filter),
                withoutTheFilter: false);
            Assert.That(saved.Split(Environment.NewLine), Has.One.Match(@"^\s*IdentificationFilter = ""PepQValue""$"));
            var reopened = LoadCommon(saved);
            Assert.That((TaskWindowIdentificationFilter.PepBoxTicked(reopened), reopened.IdentificationFilterFromLegacyThresholds), Is.EqualTo((true, false)));

            // A new task keeps the new default: the box unticked, saved as the notch.
            var fresh = new SearchTask().CommonParameters;
            Assert.That(TaskWindowIdentificationFilter.PepBoxTicked(fresh), Is.False);
            Assert.That(TaskWindowIdentificationFilter.ForSave(false, TaskWindowIdentificationFilter.QValueThresholdText(fresh), TaskWindowIdentificationFilter.PepQValueThresholdText(fresh)),
                Is.EqualTo((FilterType.QValueNotch, 0.01, 1.0)));
        }

        [Test]
        public static void AnOldDefaultTask_LoadsAsTheNotch()
        {
            var common = LoadCommon(TomlOf(new CommonParameters(), withoutTheFilter: true));

            Assert.That((common.IdentificationFilter, common.QValueThreshold, common.PepQValueThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01, 1.0)));
            Assert.That(common.IdentificationFilterFromLegacyThresholds, Is.False);
            Assert.That(SearchTask.LegacyFilterSentence(common), Is.Null);
        }

        [Test]
        public static void ATaskThatNamesTheNotch_KeepsIt_WhateverItsThresholds()
        {
            var named = new CommonParameters(qValueThreshold: 1.0, pepQValueThreshold: 0.01, identificationFilter: FilterType.QValueNotch);
            var common = LoadCommon(TomlOf(named, withoutTheFilter: false));

            Assert.That(common.IdentificationFilter, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(common.IdentificationFilterFromLegacyThresholds, Is.False);
        }

        [Test]
        public static void AnOldGuiPepTask_SavedAgain_NamesItsFilter()
        {
            var old = Nett.Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            string saved = TomlOf(old.CommonParameters, withoutTheFilter: false);

            Assert.That(saved.Split(Environment.NewLine), Has.One.Match(@"^\s*IdentificationFilter = ""PepQValue""$"));
            var reloaded = LoadCommon(saved);
            Assert.That(reloaded.IdentificationFilter, Is.EqualTo(FilterType.PepQValue));
            Assert.That(reloaded.IdentificationFilterFromLegacyThresholds, Is.False, "the rule applies once per old file");
        }

        [Test]
        public static void SearchTask_AnOldGuiPepTask_WarnsOnce_AndResultsTxtSaysWhereTheFilterCameFrom()
        {
            var common = LoadCommon(TomlOf(Chosen(FilterType.PepQValue, qValueThreshold: 1.0, pepQValueThreshold: 0.01), withoutTheFilter: true));
            Assert.That(common.IdentificationFilterFromLegacyThresholds, Is.True);

            var (folder, results, warnings) = Search("legacyToml", common, trainPep: false);
            const string sentence = "IdentificationFilter read from this task's legacy thresholds (PepQValueThreshold 0.01 < QValueThreshold 1): PEP q-value at 0.01";
            Assert.That(warnings, Has.One.EqualTo(sentence));
            Assert.That(results, Has.One.EqualTo(sentence));
            Assert.That(results, Has.Some.StartWith("All target PSMs with pep q-value <= 0.01: "));
            Assert.That(warnings, Has.None.Contain("keeps every match"));

            Directory.Delete(folder, true);
        }

        // Review of #2982 (nbollis), the peptide/PSM PepTrained mismatch. With more than 1000 matches but no more than 1000
        // peptides, the FDR engine trained PEP but computed no peptide-level PEP q-values, so the PSM level was stamped trained
        // and the peptide level not: under the PEP filter the PSM tables used PEP and the peptide tables the notch. Both levels
        // now follow the training.
        [Test]
        public static void PepTrainedWithFewPeptides_BothLevelsFollowTheTraining()
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
            psms.ForEach(p => p.ResolveAllAmbiguities());
            int peptides = psms.Select(p => p.FullSequence).Distinct().Count();
            Assume.That(psms.Count, Is.GreaterThan(peptides), "the fixture needs more matches than peptides");

            var threshold = typeof(FdrAnalysisEngine).GetProperty(nameof(FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue))!;
            int original = FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                try
                {
                    // More matches than the threshold, but no more peptides: the engine's second branch.
                    threshold.SetValue(null, peptides);
                    commonParameters.QValueCutoffForPepCalculation = 0.1;
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();
                }
                finally
                {
                    threshold.SetValue(null, original);
                }
            }

            Assert.That(psms.Select(p => p.PsmFdrInfo.PepTrained).Distinct(), Is.EqualTo(new[] { true }), "the fixture trains PEP");
            Assert.That(psms.Select(p => p.PeptideFdrInfo.PepTrained).Distinct(), Is.EqualTo(new[] { true }), "the peptide level follows the training");
            var peptideStandIns = FilteredPsms.Filter(psms, commonParameters, includeDecoys: true, includeContaminants: true, includeAmbiguous: true,
                includeHighQValuePsms: true, filterAtPeptideLevel: true).FilteredPsmsList;
            Assert.That(peptideStandIns.Select(p => p.PeptideFdrInfo.PEP_QValue), Has.All.LessThan(2), "peptide-level PEP q-values were computed");
        }

        // Review of #2982 (nbollis), SetFdrValues: replacing a match's q-values keeps the record of whether its PEP was trained.
        [Test]
        public static void SetFdrValues_KeepsTheRecordOfWhetherPepWasTrained()
        {
            var psm = Matches("PEPTIDEK", "trypsin", 1, 1, pepTrained: true)[0];
            Assume.That(psm.PsmFdrInfo.PepTrained, Is.True);

            psm.SetFdrValues(1, 0, 0.002, 1, 0, 0.003, 0.01, 0.004);

            Assert.That(psm.PsmFdrInfo.PepTrained, Is.True);
            Assert.That((psm.FdrInfo.QValue, psm.FdrInfo.QValueNotch, psm.FdrInfo.PEP_QValue), Is.EqualTo((0.002, 0.003, 0.004)));
        }

        // Review of #2982 (nbollis), rounding: results.txt states the threshold exactly; 0.005 is not 0.01.
        [Test]
        public static void SearchTask_ResultsTxt_StatesTheThresholdExactly()
        {
            var (folder, results, _) = Search("exactThreshold", Chosen(FilterType.QValueNotch, qValueThreshold: 0.005, pepQValueThreshold: 1.0), trainPep: false);

            var countLines = results.Where(l => l.StartsWith("All target ") && l.Contains(" with q-value notch <= ")).ToList();
            Assert.That(countLines, Has.Count.EqualTo(3), "PSMs, peptides and protein groups");
            Assert.That(countLines, Has.All.Contain("q-value notch <= 0.005"));

            Directory.Delete(folder, true);
        }

        // Review of #2982 (nbollis): the pooled q-value is no filter, so a caller cannot ask for it by name either.
        [Test]
        public static void ThePooledQValue_CannotBeAskedForByName()
        {
            var psms = Matches("PEPTIDEK", "trypsin", 10, 10, pepTrained: true);

            Assert.Throws<ArgumentException>(() => FilteredPsms.Filter(psms, new CommonParameters(), filterType: FilterType.QValue));
            Assert.Throws<ArgumentException>(() => FilteredPsms.Passes(psms[0], false, FilterType.QValue, 0.01));
        }
    }
}
