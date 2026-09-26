using Chemistry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using PredictionClients.Koina.AbstractClasses;
using Omics.Digestion;
using EngineLayer.DatabaseLoading;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The opt-in spectral-angle path in PostSearchAnalysisTask.ComputeSpectrumSimilarity.
    ///
    /// Every PSM here has its ambiguities resolved, because GetSpectralMatchesWithoutComputedSpectralAngle
    /// filters out any PSM whose FullSequence is null - an unresolved PSM never reaches the code under test.
    ///
    /// Nothing here touches the network except the one [Category("ExternalService")] test. The remote
    /// call is replaced through PostSearchAnalysisTask.SpectrumPredictor, and the real Prosit wrapper is
    /// exercised with an input the model rejects before it builds a request.
    /// </summary>
    [TestFixture]
    public class SupplementalSpectralSimilarityTests
    {
        private static readonly CommonParameters CommonParams = new(
            digestionParams: new DigestionParams(protease: "trypsin"),
            scoreCutoff: 1,
            productMassTolerance: new PpmTolerance(20),
            precursorMassTolerance: new PpmTolerance(5));

        private readonly List<(SpectralLibrary Library, string Path)> _libraries = new();
        private readonly List<string> _warnings = new();

        private void CaptureWarning(object sender, StringEventArgs e) => _warnings.Add(e.S);

        [SetUp]
        public void SetUp()
        {
            _warnings.Clear();
            MetaMorpheusTask.WarnHandler += CaptureWarning;
        }

        [TearDown]
        public void TearDown()
        {
            MetaMorpheusTask.WarnHandler -= CaptureWarning;
            foreach (var (library, path) in _libraries)
            {
                library.CloseConnections();
                File.Delete(path);
            }
            _libraries.Clear();
        }

        /// <summary>
        /// A PSM whose ambiguities are resolved, so FullSequence is populated and the PSM actually
        /// survives the filter. hcdEnergy null means the scan carries no recorded energy.
        /// </summary>
        private static SpectralMatch ResolvedPsm(string sequence, int charge, string hcdEnergy = null)
        {
            // digested from a real protein rather than constructed directly: ResolveAllAmbiguities
            // reaches through to the parent to decide decoy status, and a parentless peptide NREs
            var peptide = new Protein(sequence, "accession")
                .Digest(CommonParams.DigestionParams, new List<Omics.Modifications.Modification>(),
                        new List<Omics.Modifications.Modification>())
                .First(p => p.BaseSequence == sequence);

            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            var matched = products.Take(5)
                .Select(p => new MatchedFragmentIon(p, p.ToMz(1), 100, 1))
                .ToList();

            var scan = new MsDataScan(
                massSpectrum: new MzSpectrum(matched.Select(m => m.Mz).ToArray(), matched.Select(_ => 100.0).ToArray(), false),
                oneBasedScanNumber: 1, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.0, scanWindowRange: new MzRange(50, 2000), scanFilter: "FTMS + p NSI",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1000, injectionTime: null,
                noiseData: null, nativeId: "scan=1", hcdEnergy: hcdEnergy);

            var psm = new PeptideSpectralMatch(peptide, 0, 50, 0,
                new Ms2ScanWithSpecificMass(scan, peptide.MonoisotopicMass.ToMz(charge), charge, "test.raw", CommonParams),
                CommonParams, matched);

            psm.ResolveAllAmbiguities();
            psm.SpectralAngle = -1;
            return psm;
        }

        /// <summary>
        /// A task with predictions switched on for an HCD peptide search, and the remote call stubbed to
        /// return nothing.
        /// </summary>
        private static PostSearchAnalysisTask TaskWith(params SpectralMatch[] psms)
        {
            var task = new PostSearchAnalysisTask
            {
                CommonParameters = CommonParams,
                SpectrumPredictor = _ => new List<LibrarySpectrum>()
            };
            task.Parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { UsePredictedSpectraForSpectralAngle = true },
                AllSpectralMatches = psms.ToList(),
                OutputFolder = Path.GetTempPath(),
                SearchTaskId = "test",
                SearchTaskResults = new MyTaskResults(task)
            };
            return task;
        }

        /// <summary>A PSM whose matched ions are replaced, keeping its annotations.</summary>
        private static void SetMatchedIons(SpectralMatch psm, IEnumerable<MatchedFragmentIon> ions)
        {
            var list = ions.ToList();
            psm.MatchedFragmentIons.Clear();
            psm.MatchedFragmentIons.AddRange(list);
        }

        private static MatchedFragmentIon WithMzAndIntensity(MatchedFragmentIon ion, double mz, double intensity) =>
            new(ion.NeutralTheoreticalProduct, mz, intensity, ion.Charge);

        private static LibrarySpectrum LibrarySpectrumFor(SpectralMatch psm) =>
            new(psm.FullSequence, psm.ScanPrecursorMonoisotopicPeakMz, psm.ScanPrecursorCharge,
                psm.MatchedFragmentIons, 1.0);

        /// <summary>
        /// A real, file-backed library. SpectralLibrary.GetAllLibrarySpectra reads spectra from disk
        /// by byte offset, so setting Results on a default-constructed instance does not produce a
        /// readable library - it throws on enumeration. The file is deleted in TearDown.
        /// </summary>
        private SpectralLibrary LibraryOf(params SpectralMatch[] psms)
        {
            string path = Path.Combine(Path.GetTempPath(), "spectralAngle_" + Guid.NewGuid().ToString("N") + ".msp");
            File.WriteAllLines(path, psms.Select(psm => LibrarySpectrumFor(psm).ToString()));
            var library = new SpectralLibrary(new List<string> { path });
            _libraries.Add((library, path));
            return library;
        }

        // ---------- the opt-in gate ----------

        /// <summary>
        /// Prediction is a call to Koina, a third-party web service, and turning it on changes the PEP
        /// features. A search must not start doing either unless the user asked, so the default stays off.
        /// </summary>
        [Test]
        public void PredictedSpectraAreOffByDefault()
        {
            Assert.That(new SearchParameters().UsePredictedSpectraForSpectralAngle, Is.False);
        }

        /// <summary>
        /// Off means the method does nothing: no prediction is requested, and not even a PSM the library
        /// covers is rescored, so a default search is identical to one without the feature. The stub
        /// throws if asked, so a regression that re-enables the call fails here.
        /// </summary>
        [Test]
        public void WithPredictionsOffNothingIsRequestedOrScored()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);
            task.Parameters.SearchParameters.UsePredictedSpectraForSpectralAngle = false;
            bool predictorCalled = false;
            task.SpectrumPredictor = _ => { predictorCalled = true; throw new InvalidOperationException("must not be called"); };

            task.ComputeSpectrumSimilarity(LibraryOf(covered));

            Assert.That(predictorCalled, Is.False, "no prediction may be requested with the flag off");
            Assert.That(covered.SpectralAngle, Is.EqualTo(-1));
            Assert.That(uncovered.SpectralAngle, Is.EqualTo(-1));
            Assert.That(_warnings, Is.Empty);
        }

        // ---------- the scoring loop ----------

        /// <summary>A PSM the library covers comes out with a real angle.</summary>
        [Test]
        public void PsmCoveredByTheLibraryGetsAComputedAngle()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);

            task.ComputeSpectrumSimilarity(LibraryOf(psm));

            Assert.That(psm.SpectralAngle, Is.InRange(0.0, 1.0),
                "a PSM present in the library must be scored, not left at the sentinel");
        }

        /// <summary>
        /// A PSM the library misses is scored against its predicted spectrum, and only that PSM is sent
        /// for prediction.
        /// </summary>
        [Test]
        public void PsmTheLibraryMissedIsScoredAgainstItsPrediction()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);
            List<FragmentIntensityPredictionInput> requested = null;
            task.SpectrumPredictor = inputs =>
            {
                requested = inputs;
                return new List<LibrarySpectrum> { LibrarySpectrumFor(uncovered) };
            };

            task.ComputeSpectrumSimilarity(LibraryOf(covered));

            Assert.That(requested.Select(i => (i.FullSequence, i.PrecursorCharge)),
                Is.EquivalentTo(new[] { (uncovered.FullSequence, 3) }));
            Assert.That(covered.SpectralAngle, Is.InRange(0.0, 1.0));
            Assert.That(uncovered.SpectralAngle, Is.InRange(0.0, 1.0));
        }

        /// <summary>
        /// A PSM with neither a library spectrum nor a prediction keeps the sentinel rather than being
        /// scored against somebody else's spectrum.
        /// </summary>
        [Test]
        public void PsmMissingFromANonEmptyLookupKeepsTheSentinel()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);

            task.ComputeSpectrumSimilarity(LibraryOf(covered));

            Assert.That(covered.SpectralAngle, Is.GreaterThan(-1));
            Assert.That(uncovered.SpectralAngle, Is.EqualTo(-1),
                "a PSM the lookup does not contain must not borrow another peptide's spectrum");
        }

        /// <summary>
        /// SpectralAngle 0 is a legitimate (terrible) score, and only a negative value means
        /// "not computed". A PSM that already has an angle must not be recomputed.
        /// </summary>
        [Test]
        public void AlreadyScoredPsmsAreLeftAlone()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            psm.SpectralAngle = 0;
            var task = TaskWith(psm);

            task.ComputeSpectrumSimilarity(LibraryOf(psm));

            Assert.That(psm.SpectralAngle, Is.EqualTo(0), "zero is a score, not a sentinel");
        }

        // ---------- the angle itself ----------

        /// <summary>
        /// Ions are paired by annotation, so an m/z offset between the matched ions and the prediction does
        /// not move the angle. Identical intensities at 50 ppm (wider than a fixed 20 ppm window) score 1.
        /// Pairing by m/z scored exactly this case 0.
        /// </summary>
        [Test]
        public void MassErrorDoesNotChangeTheAngle()
        {
            var psm = ResolvedPsm("ELVISLIVESK", 2);
            var ideal = psm.MatchedFragmentIons.Select((m, i) => WithMzAndIntensity(m, m.Mz, 100 + 50 * i)).ToList();
            SetMatchedIons(psm, ideal.Select(m => WithMzAndIntensity(m, m.Mz * (1 + 50e-6), m.Intensity)));
            var task = TaskWith(psm);
            task.SpectrumPredictor = _ => new List<LibrarySpectrum>
            {
                new(psm.FullSequence, psm.ScanPrecursorMonoisotopicPeakMz, 2, ideal, 1.0)
            };

            task.ComputeSpectrumSimilarity(null);

            Assert.That(psm.SpectralAngle, Is.EqualTo(1).Within(1e-7));
        }

        /// <summary>
        /// A value, not a range. Two library ions of equal intensity, one observed: the square-root
        /// vectors are (1, 0) and (1, 1), the cosine is 1/sqrt(2), so the normalized angle is exactly 0.5.
        /// An experimental ion the library lacks is ignored, as in the spectral-library search.
        /// </summary>
        [Test]
        public void HalfTheLibraryObservedScoresOneHalf()
        {
            var ions = ResolvedPsm("ELVISLIVESK", 2).MatchedFragmentIons;
            var library = new List<MatchedFragmentIon> { WithMzAndIntensity(ions[0], ions[0].Mz, 4), WithMzAndIntensity(ions[1], ions[1].Mz, 4) };
            var experimental = new List<MatchedFragmentIon> { WithMzAndIntensity(ions[0], ions[0].Mz, 9), WithMzAndIntensity(ions[2], ions[2].Mz, 1000) };

            Assert.That(PostSearchAnalysisTask.SpectralAngleByAnnotation(experimental, library), Is.EqualTo(0.5).Within(1e-12));
        }

        /// <summary>The same fragment at another charge is a different peak and does not pair.</summary>
        [Test]
        public void ChargeIsPartOfTheAnnotation()
        {
            var ion = ResolvedPsm("ELVISLIVESK", 2).MatchedFragmentIons[0];
            var library = new List<MatchedFragmentIon> { ion };
            var experimental = new List<MatchedFragmentIon> { new(ion.NeutralTheoreticalProduct, ion.Mz, ion.Intensity, 2) };

            Assert.That(PostSearchAnalysisTask.SpectralAngleByAnnotation(experimental, library), Is.EqualTo(0));
        }

        /// <summary>
        /// A library or predicted spectrum with no ions cannot be compared: the PSM keeps the sentinel rather
        /// than throwing between the search and FDR.
        /// </summary>
        [Test]
        public void AnEmptyLibrarySpectrumLeavesTheSentinel()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);
            task.SpectrumPredictor = _ => new List<LibrarySpectrum>
            {
                new(psm.FullSequence, psm.ScanPrecursorMonoisotopicPeakMz, 2, new List<MatchedFragmentIon>(), 1.0)
            };

            Assert.DoesNotThrow(() => task.ComputeSpectrumSimilarity(null));
            Assert.That(psm.SpectralAngle, Is.EqualTo(-1));
        }

        // ---------- what Prosit HCD can score ----------

        /// <summary>
        /// Prosit 2020 HCD predicts b and y ions. On any other fragmentation the matched ions never pair, and
        /// every PSM would carry a 0 that PEP reads as real, so nothing is requested or scored and the user is
        /// told, in results.txt as well as the log.
        /// </summary>
        [TestCase(DissociationType.ETD)]
        [TestCase(DissociationType.EThcD)]
        [TestCase(DissociationType.LowCID)]
        [TestCase(DissociationType.Autodetect)]
        public void DissociationTypesPrositCannotPredictAreSkipped(DissociationType dissociationType)
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);
            task.CommonParameters = new CommonParameters(dissociationType: dissociationType);
            task.SpectrumPredictor = _ => throw new InvalidOperationException("must not be called");

            task.ComputeSpectrumSimilarity(LibraryOf(psm));

            Assert.That(psm.SpectralAngle, Is.EqualTo(-1));
            Assert.That(_warnings, Has.Exactly(1).Contains("not computed").And.Contains(dissociationType.ToString()));
            Assert.That(task.Parameters.SearchTaskResults.ToString(), Does.Contain("not computed"));
        }

        /// <summary>A file-specific override to ETD disqualifies the search just as the task setting would.</summary>
        [Test]
        public void AFileSpecificDissociationTypeIsCheckedToo()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);
            task.FileSpecificParameters = new List<(string, CommonParameters)> { ("etd.raw", new CommonParameters(dissociationType: DissociationType.ETD)) };
            task.SpectrumPredictor = _ => throw new InvalidOperationException("must not be called");

            task.ComputeSpectrumSimilarity(null);

            Assert.That(_warnings, Has.Exactly(1).Contains("ETD"));
        }

        /// <summary>CID is close enough to be useful, but the user is told the model is HCD.</summary>
        [Test]
        public void CidIsScoredWithAWarning()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);
            task.CommonParameters = new CommonParameters(dissociationType: DissociationType.CID);

            task.ComputeSpectrumSimilarity(LibraryOf(psm));

            Assert.That(psm.SpectralAngle, Is.EqualTo(1).Within(1e-7));
            Assert.That(_warnings, Has.Exactly(1).Contains("approximation"));
        }

        /// <summary>An oligo search must not send its sequences to Koina as if they were peptides.</summary>
        [Test]
        public void OligoSearchesAreSkipped()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);
            task.SpectrumPredictor = _ => throw new InvalidOperationException("must not be called");
            var analyteType = GlobalVariables.AnalyteType;
            try
            {
                GlobalVariables.AnalyteType = AnalyteType.Oligo;
                task.ComputeSpectrumSimilarity(null);
            }
            finally
            {
                GlobalVariables.AnalyteType = analyteType;
            }

            Assert.That(_warnings, Has.Exactly(1).Contains("Oligo"));
        }

        // ---------- the record in results.txt ----------

        /// <summary>
        /// With the flag on, q-values depend on whether Koina answered. The counts go into results.txt so
        /// two runs that disagree can be explained afterwards.
        /// </summary>
        [Test]
        public void TheOutcomeIsRecordedInTheTaskSummary()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);
            task.SpectrumPredictor = _ => new List<LibrarySpectrum> { LibrarySpectrumFor(uncovered) };

            task.ComputeSpectrumSimilarity(LibraryOf(covered));

            Assert.That(task.Parameters.SearchTaskResults.ToString(), Does.Contain(
                "2 spectral matches without an angle, 1 spectral library entries available; 1 peptide/charge pairs requested, 1 predicted, 0 rejected by the model; 2 angles assigned."));
        }

        [Test]
        public void AFailedRequestIsRecordedInTheTaskSummary()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);
            task.SpectrumPredictor = _ => throw new HttpRequestException("Koina is down");

            task.ComputeSpectrumSimilarity(null);

            Assert.That(task.Parameters.SearchTaskResults.ToString(), Does.Contain("request failed (HttpRequestException: Koina is down)"));
        }

        /// <summary>
        /// The prediction-to-library conversion for a MODIFIED peptide, offline. The FullSequence is the one
        /// MetaMorpheus writes (fixed carbamidomethyl C, variable oxidised M); under MapToInputFullSequence
        /// mzLib rebuilds the peptide from that string, and the spectrum must come back under the same key
        /// the PSM looks up with. A throw here would land in the catch and silently drop every prediction.
        /// </summary>
        [Test]
        public void AModifiedPredictionComesBackUnderThePsmKey()
        {
            var carbamidomethyl = GlobalVariables.AllModsKnown.First(m => m.IdWithMotif == "Carbamidomethyl on C" && m.ModificationType == "Common Fixed");
            var oxidation = GlobalVariables.AllModsKnown.First(m => m.IdWithMotif == "Oxidation on M" && m.ModificationType == "Common Variable");
            var peptide = new Protein("PEPCTIDEMK", "accession")
                .Digest(CommonParams.DigestionParams, new List<Omics.Modifications.Modification> { carbamidomethyl },
                        new List<Omics.Modifications.Modification> { oxidation })
                .First(p => p.AllModsOneIsNterminus.Count == 2);
            var model = new PredictedProsit(new PeptideFragmentIntensityPrediction(
                peptide.FullSequence, peptide.FullSequence, 2,
                new List<string> { "b2+1", "y3+1", "y4+1" }, new List<double> { 0, 0, 0 }, new List<double> { 0.2, 1.0, 0.5 }));

            var spectra = TaskWith().LibrarySpectraFrom(model);

            Assert.That(spectra, Has.Count.EqualTo(1));
            Assert.That(spectra[0].Sequence, Is.EqualTo(peptide.FullSequence));
            Assert.That(spectra[0].MatchedFragmentIons.Select(i => i.Annotation), Is.EquivalentTo(new[] { "b2+1", "y3+1", "y4+1" }));
            var y3 = spectra[0].MatchedFragmentIons.Single(i => i.Annotation == "y3+1");
            Assert.That(y3.NeutralTheoreticalProduct.NeutralMass,
                Is.EqualTo(new PeptideWithSetModifications("EM[Common Variable:Oxidation on M]K", GlobalVariables.AllModsKnownDictionary).MonoisotopicMass).Within(1e-6),
                "the oxidised M must be in the fragment masses");
        }

        /// <summary>
        /// Whole searches with the flag on that the feature cannot serve, from the call site to results.txt.
        /// Semi-specific searches run FDR before post-search analysis, and an ETD search has no b/y ions to
        /// pair; both must say so in results.txt rather than silently match a run with the flag off. Neither
        /// reaches the network, so this runs offline.
        /// </summary>
        [TestCase("semi")]
        [TestCase("etd")]
        public void ASearchTheFeatureCannotServeSaysSoInResults(string kind)
        {
            var searchTask = kind == "semi"
                ? new SearchTask
                {
                    SearchParameters = new SearchParameters
                    {
                        SearchType = SearchType.NonSpecific,
                        LocalFdrCategories = new List<FdrCategory> { FdrCategory.FullySpecific, FdrCategory.SemiSpecific },
                        UsePredictedSpectraForSpectralAngle = true
                    },
                    CommonParameters = new CommonParameters(scoreCutoff: 11,
                        digestionParams: new DigestionParams(minPeptideLength: 7, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.N))
                }
                : new SearchTask
                {
                    SearchParameters = new SearchParameters { DoLabelFreeQuantification = false, UsePredictedSpectraForSpectralAngle = true },
                    CommonParameters = new CommonParameters(dissociationType: DissociationType.ETD)
                };
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", kind == "semi" ? "tinySemi.mgf" : "SmallCalibratible_Yeast.mzML");
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", kind == "semi" ? "semiTest.fasta" : "smalldb.fasta");
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "PredictedAngleSkipped_" + kind);
            Directory.CreateDirectory(folder);
            try
            {
                searchTask.RunTask(folder, new List<DbForTask> { new(database, false) }, new List<string> { spectra }, "");

                string results = File.ReadAllText(Path.Combine(folder, "results.txt"));
                Assert.That(results, Does.Contain("Predicted spectral angles were requested but not computed: "
                    + (kind == "semi" ? "semi- and non-specific searches" : "Prosit 2020 predicts HCD spectra, and this search uses ETD")));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>A Prosit model whose predictions are set directly, standing in for a completed request.</summary>
        private class PredictedProsit : PredictionClients.Koina.SupportedModels.FragmentIntensityModels.Prosit2020IntensityHCD
        {
            public PredictedProsit(params PeptideFragmentIntensityPrediction[] predictions)
                : base(fragmentIonMappingMode: PredictionClients.Koina.Util.FragmentIonMappingMode.MapToInputFullSequence)
            {
                Predictions = predictions.ToList();
                ValidInputsMask = predictions.Select(p => p.FragmentIntensities != null).ToArray();
            }
        }

        // ---------- failure is not fatal ----------

        /// <summary>
        /// Koina being down must cost the predicted angles, not the search: this runs before FDR, so an
        /// exception escaping here would lose every output file of a search that had already finished.
        /// The library half still scores, and the user is told why the rest did not.
        /// </summary>
        [Test]
        public void AServiceFailureWarnsAndLeavesTheSentinel()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);
            task.SpectrumPredictor = _ => throw new HttpRequestException("Koina is down");

            Assert.DoesNotThrow(() => task.ComputeSpectrumSimilarity(LibraryOf(covered)));

            Assert.That(covered.SpectralAngle, Is.GreaterThan(-1), "the library half must survive the outage");
            Assert.That(uncovered.SpectralAngle, Is.EqualTo(-1));
            Assert.That(_warnings, Has.Exactly(1).Contains("HttpRequestException").And.Contains("Koina is down"));
        }

        /// <summary>
        /// The real Prosit wrapper, offline: charge 7 is outside what Prosit accepts, so the model drops
        /// the input before building a request. The PSM keeps the sentinel and the drop is reported
        /// rather than silently lost.
        /// </summary>
        [Test]
        public void InputsPrositRejectsAreReportedAndNeverSent()
        {
            var psm = ResolvedPsm("PEPTIDEK", 7);
            var task = TaskWith(psm);
            task.SpectrumPredictor = null;

            Assert.DoesNotThrow(() => task.ComputeSpectrumSimilarity(null));

            Assert.That(psm.SpectralAngle, Is.EqualTo(-1));
            Assert.That(_warnings, Has.Exactly(1).Contains("could not predict 1 of 1"));
        }

        /// <summary>
        /// SearchTask closes the spectral library before post-search analysis runs, unless the run is
        /// updating it. Reading it here therefore hits a closed file on an ordinary search with a
        /// supplied library. That must cost the angles, not the whole completed search.
        /// </summary>
        [Test]
        public void ClosedLibraryDegradesInsteadOfLosingTheSearch()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);

            var library = LibraryOf(psm);
            library.CloseConnections();

            Assert.DoesNotThrow(() => task.ComputeSpectrumSimilarity(library));
            Assert.That(psm.SpectralAngle, Is.EqualTo(-1),
                "a closed library yields no spectra, so the PSM keeps the sentinel");
        }

        // ---------- the combined lookup ----------

        /// <summary>
        /// The lookup holds the library's spectra and the predictions for what the library missed,
        /// under the key PSMs query with, and a prediction never displaces a measured spectrum.
        /// </summary>
        [Test]
        public void CombinedLookupHoldsLibraryAndPredictedSpectra()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);
            var predictedForCovered = new LibrarySpectrum(covered.FullSequence, 500, 2, covered.MatchedFragmentIons, 10);
            task.SpectrumPredictor = _ => new List<LibrarySpectrum> { LibrarySpectrumFor(uncovered), predictedForCovered };

            var lookup = task.BuildCombinedSpectrumLookup(new List<SpectralMatch> { covered, uncovered }, LibraryOf(covered));

            Assert.That(lookup.Keys, Is.EquivalentTo(new[] { (covered.FullSequence, 2), (uncovered.FullSequence, 3) }));
            Assert.That(lookup[(covered.FullSequence, 2)], Is.Not.SameAs(predictedForCovered),
                "a measured spectrum beats a predicted one");
        }

        /// <summary>A measured spectrum beats a predicted one for the same peptide and charge.</summary>
        [Test]
        public void RealLibrarySpectraAreNotOverwrittenByPredictions()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var real = LibrarySpectrumFor(psm);
            var lookup = new Dictionary<(string, int), LibrarySpectrum> { [("PEPTIDEK", 2)] = real };
            var predicted = new LibrarySpectrum("PEPTIDEK", 500, 2, psm.MatchedFragmentIons, 10);

            PostSearchAnalysisTask.MergePredictedSpectra(new[] { predicted }, lookup);

            Assert.That(lookup[("PEPTIDEK", 2)], Is.SameAs(real));
        }

        // ---------- what would be sent for prediction ----------

        /// <summary>
        /// Only the PSMs the library did not cover are requested. Predicting one we already have a
        /// real spectrum for would be a wasted call and a worse answer.
        /// </summary>
        [Test]
        public void OnlyPsmsTheLibraryMissedAreQueuedForPrediction()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var lookup = new Dictionary<(string, int), LibrarySpectrum>
            {
                [(covered.FullSequence, covered.ScanPrecursorCharge)] = LibrarySpectrumFor(covered)
            };

            var inputs = PostSearchAnalysisTask.BuildPredictionInputs(
                new List<SpectralMatch> { covered, uncovered }, lookup);

            Assert.That(inputs.Keys, Is.EquivalentTo(new[] { (uncovered.FullSequence, 3) }));
        }

        /// <summary>
        /// The same peptide and charge seen twice is one request, not two. Duplicates are the normal
        /// case: the same peptide is routinely matched in many scans.
        /// </summary>
        [Test]
        public void DuplicateSequenceAndChargeIsRequestedOnce()
        {
            var inputs = PostSearchAnalysisTask.BuildPredictionInputs(
                new List<SpectralMatch> { ResolvedPsm("PEPTIDEK", 2), ResolvedPsm("PEPTIDEK", 2) },
                new Dictionary<(string, int), LibrarySpectrum>());

            Assert.That(inputs, Has.Count.EqualTo(1));
        }

        /// <summary>
        /// The same peptide at a different charge is a different spectrum and must be requested
        /// separately - the dedup key has to include the charge.
        /// </summary>
        [Test]
        public void SameSequenceAtADifferentChargeIsRequestedSeparately()
        {
            var inputs = PostSearchAnalysisTask.BuildPredictionInputs(
                new List<SpectralMatch> { ResolvedPsm("PEPTIDEK", 2), ResolvedPsm("PEPTIDEK", 3) },
                new Dictionary<(string, int), LibrarySpectrum>());

            Assert.That(inputs, Has.Count.EqualTo(2));
        }

        /// <summary>
        /// The energy comes from the scan through SpectralMatch.CollisionalEnergy, which every search
        /// PSM carries - not from Ms2Scan, which a search never sets.
        /// </summary>
        [Test]
        public void CollisionEnergyComesFromTheScanWhenItRecordedOne()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2, hcdEnergy: "27.0");

            Assert.That(psm.Ms2Scan, Is.Null, "a search PSM has no Ms2Scan; the energy must not depend on it");
            Assert.That(PostSearchAnalysisTask.ResolveCollisionEnergy(psm), Is.EqualTo(27));
        }

        /// <summary>
        /// 30 is a fallback, not a measurement. Pinned so a change to it is a deliberate decision
        /// about what energy unlabelled spectra are predicted at, rather than an accident.
        /// </summary>
        [Test]
        public void CollisionEnergyFallsBackToThirtyWithoutARecordedEnergy()
        {
            Assert.That(PostSearchAnalysisTask.ResolveCollisionEnergy(ResolvedPsm("PEPTIDEK", 2, hcdEnergy: null)),
                Is.EqualTo(30));
        }

        [Test]
        public void CollisionEnergyFallsBackWhenTheScanEnergyIsNotANumber()
        {
            Assert.That(PostSearchAnalysisTask.ResolveCollisionEnergy(ResolvedPsm("PEPTIDEK", 2, hcdEnergy: "not a number")),
                Is.EqualTo(30));
        }

        // ---------- live ----------

        /// <summary>
        /// End to end against Koina: an unmodified tryptic peptide comes back with a real angle. An
        /// outage is reported by the production code as a warning, which this test turns into a skip.
        /// </summary>
        [Test]
        [Category("ExternalService")]
        public void PrositPredictionProducesARealAngle()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2, hcdEnergy: "28");
            var task = TaskWith(psm);
            task.SpectrumPredictor = null;

            task.ComputeSpectrumSimilarity(null);

            if (_warnings.Any(w => w.StartsWith("Predicted spectra were unavailable")))
            {
                Assert.Ignore("Skipping external-service test: Koina unavailable. " + string.Join(" | ", _warnings));
            }
            Assert.That(psm.SpectralAngle, Is.InRange(0.0, 1.0));
        }
    }
}
