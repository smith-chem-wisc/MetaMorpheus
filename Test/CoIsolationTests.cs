using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class CoIsolationTests
    {
        #region Public Methods

        [Test]
        public static void TestCoIsolation()
        {
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNNKNDNK", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            Proteomics.Peptide pep1 = new Proteomics.Peptide("NNNK");
            Proteomics.Peptide pep2 = new Proteomics.Peptide("NDNK");

            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);

            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            IMzmlScan[] Scans = new IMzmlScan[2];
            double[] ms1intensities = new double[] { 0.8, 0.8, 0.2, 0.02, 0.2, 0.02 };
            double[] ms1mzs = dist1.Masses.Concat(dist2.Masses).OrderBy(b => b).Select(b => b.ToMz(1)).ToArray();

            double selectedIonMz = ms1mzs[1];

            MzmlMzSpectrum MS1 = new MzmlMzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MzmlScan(1, MS1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, null);

            double[] ms2intensities = new double[] { 1, 1, 1, 1, 1 };
            double[] ms2mzs = new double[] { 146.106.ToMz(1), 228.086.ToMz(1), 229.07.ToMz(1), 260.148.ToMz(1), 342.129.ToMz(1) };
            MzmlMzSpectrum MS2 = new MzmlMzSpectrum(ms2mzs, ms2intensities, false);
            double isolationMZ = selectedIonMz;
            Scans[1] = new MzmlScanWithPrecursor(2, MS2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, selectedIonMz, null, null, isolationMZ, 2.5, DissociationType.HCD, 1, null, null);

            var myMsDataFile = new FakeMsDataFile(Scans);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 50;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
            int maximumMissedCleavages = 2;
            int? minPeptideLength = null;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false, InitiatorMethionineBehavior.Variable, false,1);
            var searchResults = (SearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.Psms.Length);

            // Two matches for this single scan! Corresponding to two co-isolated masses
            Assert.AreEqual(2, searchResults.Psms[0].Length);

            Assert.IsTrue(searchResults.Psms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.Psms[0][0].ScanNumber);

            var ojdfkj = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { searchResults.Psms[0][0], searchResults.Psms[0][1] } }, proteinList, searchModes, protease, maximumMissedCleavages, null, null, InitiatorMethionineBehavior.Variable, fixedModifications, variableModifications, 4096, new List<string>(), TerminusType.None).Run();

            foreach (var huh in searchResults.Psms[0])
            {
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(ojdfkj.CompactPeptideToProteinPeptideMatching);
            }

            Assert.AreEqual("NNNK", searchResults.Psms[0][0].BaseSequence);
            Assert.AreEqual("NDNK", searchResults.Psms[0][1].BaseSequence);
        }

        #endregion Public Methods
    }
}