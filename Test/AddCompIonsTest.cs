using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal static class AddCompIonsTest
    {
        [Test]
        public static void TestAddCompIonsClassic()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("QXQ", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new OpenSearchMode();

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease protease = new Protease("Custom Protease3", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1),
                scoreCutoff: 1,
                addCompIons: false);
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            CommonParameters CommonParameters2 = new CommonParameters(
                digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1),
                scoreCutoff: 1,
                addCompIons: true);

            PeptideSpectralMatch[] allPsmsArray2 = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray2, listOfSortedms2Scans, variableModifications, fixedModifications, null, proteinList, searchModes, CommonParameters2, new List<string>()).Run();

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
        public static void TestCompIons_ModernSearch()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
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
            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease protease = new Protease("singleN4", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1), scoreCutoff: 1);
            CommonParameters withCompIons = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1), scoreCutoff: 1, addCompIons: true);

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null,
                 1, DecoyType.Reverse, CommonParameters, SearchParameters.MaxFragmentSize, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            // without complementary ions
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            // with complementary ions
            PeptideSpectralMatch[] allPsmsArray2 = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray2, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, withCompIons, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

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

        public static void TestCompIons_MatchIonsScore()
        {
            TestDataFile t = new TestDataFile();
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 300;
            //The below theoretical does not accurately represent B-Y ions
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { precursorMass + (2 * Constants.ProtonMass) - 275.1350, precursorMass + (2 * Constants.ProtonMass) - 258.127, precursorMass + (2 * Constants.ProtonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants.ProtonMass) - 147.0764, precursorMass + (2 * Constants.ProtonMass) - 70, precursorMass + (2 * Constants.ProtonMass) - 60, precursorMass + (2 * Constants.ProtonMass) - 50, 257.1244, 258.127, 275.1350 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<Product> productsWithLocalizedMassDiff = new List<Product>();
            foreach (double d in sorted_theoretical_product_masses_for_this_peptide)
            {
                NeutralTerminusFragment frag = new NeutralTerminusFragment(FragmentationTerminus.Both, d, 1, 1);
                productsWithLocalizedMassDiff.Add(new Product(ProductType.b, frag, 0));
            }
            CommonParameters commonParametersNoComp = new CommonParameters { ProductMassTolerance = new AbsoluteTolerance(0.01) };
            CommonParameters commonParametersWithComp = new CommonParameters(productMassTolerance: new AbsoluteTolerance(0.01), addCompIons: true);

            MsDataScan scan = t.GetOneBasedScan(2);
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", new CommonParameters());
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersNoComp);

            List<MatchedFragmentIon> matchedCompIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersWithComp);
            matchedCompIons.AddRange(matchedIons);

            // score when the mass-diff is on this residue
            double localizedScore = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedIons);
            double scoreNormal = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedIons);
            double scoreComp = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedCompIons);
            Assert.IsTrue(scoreNormal * 2 == scoreComp && scoreComp > scoreNormal + 1);
        }

        [Test]
        public static void TestCompIons_MatchIons()
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

            //MetaMorpheusEngine.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonSeriesT, matchedIonMassesT, matchedDaErrorT, matchedPpmErrorT, matchedIonIntensityT, precursorMass, ProductType.B, true);
            //MetaMorpheusEngine.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonSeriesF, matchedIonMassesF, matchedDaErrorF, matchedPpmErrorF, matchedIonIntensityF, precursorMass, ProductType.B, false);

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

        [Test]
        public static void AddCompIonsCommonParams()
        {
            CommonParameters cp = new CommonParameters(null, DissociationType.HCD, DissociationType.Unknown, true, true, 3, 12, true, true, 1, 
                5, 200, 0.01, null, null, false, false, true, false, null, null, null, -1, null, null, null, 1, true, 4, 1);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
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
            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease protease = new Protease("Test", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<SilacLabel>(),
                 1, DecoyType.Reverse, cp, SearchParameters.MaxFragmentSize, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(cp.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            // without complementary ions
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            var mse = new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, cp, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

        }
    }
}