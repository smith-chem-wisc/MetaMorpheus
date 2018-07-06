using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal static class AddComplementaryIonsTest
    {
        [Test]
        public static void TestAddComplementaryIonsClassic()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("QXQ", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new OpenSearchMode();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Protease protease = new Protease("Custom Protease3", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1), scoreCutoff: 1, addCompIons: false);
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            CommonParameters CommonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1), scoreCutoff: 1, addCompIons: true);
            PeptideSpectralMatch[] allPsmsArray2 = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray2, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters2, new List<string>()).Run();

            double scoreT = allPsmsArray2[0].Score;
            double scoreF = allPsmsArray[0].Score;

            // Single search mode
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            // Single ms2 scan
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            Assert.IsTrue(scoreT > 1);

            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

            Assert.IsTrue(scoreT == scoreF * 3 && scoreT > scoreF + 2);
        }

        [Test]
        public static void TestComplementaryIons_ModernSearch()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
            {
                modsDictionary.Add(mod, 0);
            }
            int ii = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)ii);
                ii++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)ii);
                ii++;
            }

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                SearchTarget = true,
            };
            Protease protease = new Protease("singleN4", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1), scoreCutoff: 1);
            CommonParameters withCompIons = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1), scoreCutoff: 1, addCompIons: true);

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications,
                new List<ProductType> { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            // without complementary ions
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            // with complementary ions
            PeptideSpectralMatch[] allPsmsArray2 = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray2, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, withCompIons, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            // Single ms2 scan
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);
            Assert.That(allPsmsArray[0] != null);
            Assert.That(allPsmsArray2[0] != null);

            Assert.IsTrue(allPsmsArray2[0].Score > 1);

            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

            Assert.IsTrue(allPsmsArray2[0].Score <= allPsmsArray[0].Score * 2 && allPsmsArray2[0].Score > allPsmsArray[0].Score + 3);
        }

        [Test]
        public static void TestComplementaryIons_MatchIonsScore()
        {
            TestDataFile t = new TestDataFile();
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 300;
            //The below theoretical does not accurately represent B-Y ions
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { precursorMass + (2 * Constants.ProtonMass) - 275.1350, precursorMass + (2 * Constants.ProtonMass) - 258.127, precursorMass + (2 * Constants.ProtonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants.ProtonMass) - 147.0764, precursorMass + (2 * Constants.ProtonMass) - 70, precursorMass + (2 * Constants.ProtonMass) - 60, precursorMass + (2 * Constants.ProtonMass) - 50, 257.1244, 258.127, 275.1350 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            double scoreT = MetaMorpheusEngine.CalculatePeptideScoreOld(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, precursorMass, new List<DissociationType> { DissociationType.HCD }, true, 0);
            double scoreF = MetaMorpheusEngine.CalculatePeptideScoreOld(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, precursorMass, new List<DissociationType> { DissociationType.HCD }, false, 0);
            Assert.IsTrue(scoreT == scoreF * 2 && scoreT > scoreF + 1);
        }

        [Test]
        public static void TestComplementaryIons_MatchIons()
        {
            TestDataFile t = new TestDataFile(0.0001);
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 402.18629720155;
            //The below theoretical does not accurately represent B-Y ions
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { 50, 60, 70, 147.0764 - Constants.ProtonMass, 200, 215, 230, 245, precursorMass + Constants.ProtonMass - 147.0764, 258.127, 275.1350, precursorMass + (2 * Constants.ProtonMass) - 70, precursorMass + (2 * Constants.ProtonMass) - 60, precursorMass + (2 * Constants.ProtonMass) - 50 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<double> matchedIonMassesT = new List<double>();
            List<double> matchedDaErrorT = new List<double>();
            List<double> matchedPpmErrorT = new List<double>();
            List<double> matchedIonIntensityT = new List<double>();
            List<double> matchedIonMassesF = new List<double>();
            List<double> matchedDaErrorF = new List<double>();
            List<double> matchedPpmErrorF = new List<double>();
            List<double> matchedIonIntensityF = new List<double>();

            List<int> matchedIonSeriesT = new List<int>();
            List<int> matchedIonSeriesF = new List<int>();

            MetaMorpheusEngine.MatchIonsOld(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonSeriesT, matchedIonMassesT, matchedDaErrorT, matchedPpmErrorT, matchedIonIntensityT, precursorMass, ProductType.BnoB1ions, true);
            MetaMorpheusEngine.MatchIonsOld(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonSeriesF, matchedIonMassesF, matchedDaErrorF, matchedPpmErrorF, matchedIonIntensityF, precursorMass, ProductType.BnoB1ions, false);

            //Test the number of series is doubled
            Assert.IsTrue(matchedIonSeriesT.Count == matchedIonSeriesF.Count * 2);
            //Test the number of ions is doubled
            Assert.IsTrue(matchedIonMassesT.Count == matchedIonMassesF.Count * 2);
            //Test the number of da errors is doubled
            Assert.IsTrue(matchedDaErrorT.Count == matchedDaErrorF.Count * 2);
            //test the number of ppm errors is doubled
            Assert.IsTrue(matchedPpmErrorT.Count == matchedPpmErrorF.Count * 2);
            //test the number of the intensity values is doubled
            Assert.IsTrue(matchedIonIntensityT.Count == matchedIonIntensityF.Count * 2);
            foreach (double d in matchedDaErrorF)
            {
                Assert.IsTrue(d <= 0.01);
            }
            foreach (double d in matchedDaErrorT)
            {
                Assert.IsTrue(d <= 0.01);
            }
        }
    }
}