using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using IO.MzML;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using UsefulProteomicsDatabases;
using Chemistry;
using System;
using MassSpectrometry;

namespace Test
{
    [TestFixture]
    public class MatchIonsOfAllCharges
    {
        [Test]
        public static void TestMatchIonsOfAllChargesBottomUp()
        {
            CommonParameters CommonParameters = new CommonParameters();

            MetaMorpheusTask.DetermineAnalyteType(CommonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();

            var proteinList = new List<Protein>
            {
                new Protein("AAAHSSLK", ""),new Protein("RQPAQPR", ""),new Protein("EKAEAEAEK", "")
            };
            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"));


            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //search by new method of looking for all charges 
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchMode, CommonParameters, null, null, new List<string>(), true).Run();
            var psm = allPsmsArray.Where(p => p != null).ToList();
            Assert.That(psm[1].MatchedFragmentIons.Count == 14);
            //there are ions with same product type and same fragment number but different charges 
            Assert.That(psm[1].MatchedFragmentIons[8].NeutralTheoreticalProduct.ProductType == psm[1].MatchedFragmentIons[9].NeutralTheoreticalProduct.ProductType &&
                psm[1].MatchedFragmentIons[8].NeutralTheoreticalProduct.FragmentNumber == psm[1].MatchedFragmentIons[9].NeutralTheoreticalProduct.FragmentNumber &&
                psm[1].MatchedFragmentIons[8].Charge != psm[1].MatchedFragmentIons[9].Charge);
            Assert.That(psm[2].MatchedFragmentIons.Count == 14);
            Assert.That(psm[4].MatchedFragmentIons.Count == 16);

            //search by old method of looking for only one charge 
            PeptideSpectralMatch[] allPsmsArray_oneCharge = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray_oneCharge, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchMode, CommonParameters, null, null, new List<string>(), false).Run();
            var psm_oneCharge = allPsmsArray_oneCharge.Where(p => p != null).ToList();

            //compare 2 scores , they should have same integer part but new search has a little higher score than old search
            Assert.That(psm[1].Score > psm_oneCharge[1].Score);
            Assert.AreEqual(Math.Truncate(psm[1].Score), 12);
            Assert.AreEqual(Math.Truncate(psm_oneCharge[1].Score), 12);

            //compare 2 results and evaluate the different matched ions
            var peptideTheorProducts = new List<Product>();
            Assert.That(psm_oneCharge[1].MatchedFragmentIons.Count == 12);
            var differences = psm[1].MatchedFragmentIons.Except(psm_oneCharge[1].MatchedFragmentIons);
            psm[1].BestMatchingPeptides.First().Peptide.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
            foreach(var ion in differences)
            {
                foreach(var product in peptideTheorProducts)
                {
                    if (product.Annotation.ToString().Equals(ion.NeutralTheoreticalProduct.Annotation.ToString()))
                    {
                        //to see if the different matched ions are qualified
                        Assert.That(CommonParameters.ProductMassTolerance.Within(ion.Mz.ToMass(ion.Charge), product.NeutralMass));
                    }
                }
            }


        }

        [Test]
        public static void TestMatchIonsOfAllChargesTopDown()
        {
            CommonParameters CommonParameters = new CommonParameters(
               digestionParams: new DigestionParams(protease: "top-down"),
               scoreCutoff: 1,
               assumeOrphanPeaksAreZ1Fragments: false);

            MetaMorpheusTask.DetermineAnalyteType(CommonParameters);

            // test output file name (should be proteoform and not peptide)
            Assert.That(GlobalVariables.AnalyteType == "Proteoform");

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein>
            {
                new Protein("MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", "P40312")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //search by new method of looking for all charges 
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchMode, CommonParameters, null, null, new List<string>(), true).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 62);
           

            //search by old method of looking for only one charge 
            PeptideSpectralMatch[] allPsmsArray_oneCharge = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray_oneCharge, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchMode, CommonParameters, null, null, new List<string>(), false).Run();

            var psm_oneCharge = allPsmsArray_oneCharge.Where(p => p != null).FirstOrDefault();
            Assert.That(psm_oneCharge.MatchedFragmentIons.Count == 47);

            //compare 2 scores , they should have same integer but new search has a little higher score than old search
            Assert.That(psm.Score > psm_oneCharge.Score);
            Assert.AreEqual(Math.Truncate(psm.Score), 47);
            Assert.AreEqual(Math.Truncate(psm_oneCharge.Score), 47);

            //compare 2 results and evaluate the different matched ions
            var peptideTheorProducts = new List<Product>();
            var differences = psm.MatchedFragmentIons.Except(psm_oneCharge.MatchedFragmentIons);
            psm.BestMatchingPeptides.First().Peptide.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
            foreach (var ion in differences)
            {
                foreach (var product in peptideTheorProducts)
                {
                    if (product.Annotation.ToString().Equals(ion.NeutralTheoreticalProduct.Annotation.ToString()))
                    {
                        //to see if the different matched ions are qualified
                        Assert.That(CommonParameters.ProductMassTolerance.Within(ion.Mz.ToMass(ion.Charge), product.NeutralMass));
                    }
                }
            }
        }

        [Test]
        public static void TestLookingForAcceptableIsotopicEnvelopes()
        {
            CommonParameters CommonParameters = new CommonParameters();

            MetaMorpheusTask.DetermineAnalyteType(CommonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();

            var proteinList = new List<Protein>
            {
                new Protein("AAAHSSLK", ""),new Protein("RQPAQPR", ""),new Protein("EKAEAEAEK", "")
            };
            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"));


            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            var ms2ScanTest = listOfSortedms2Scans[0];

            //test when all the masses are smaller than min mass
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMinimum(100),0);
            //test when all the mass are bigger than min mass
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMinimum(582), ms2ScanTest.ExperimentalFragments.Length);

            //test normal conditions, find the closest index which is bigger than the given min mass 
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMinimum(100.1),2);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMinimum(111.1), 9);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMinimum(152.05972487589443),35);


            //test when all the masses are bigger than the given max mass
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMaximum(100), -1);
            //test when all the masses are smaller than the given max mass
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMaximum(582), ms2ScanTest.ExperimentalFragments.Length - 1);

            //test normal conditions, find the closest index which is smaller than the given max mass 
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMaximum(100.2), 1);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMaximum(111.4), 8);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMaximum(152.05972487589443), 35);

            //test when all the masses are not in the given range
            //test when all the masses are too small
            var test1 = ms2ScanTest.GetClosestExperimentalIsotopicEnvelopeList(50, 95);
            Assert.AreEqual(test1, null);
            //test when all the masses are too big
            var test2 = ms2ScanTest.GetClosestExperimentalIsotopicEnvelopeList(582, 682);
            Assert.AreEqual(test2, null);
            //test when the mass which is bigger than given min mass is bigger than the mass which is smaller than the given max mass
            //for example: the mass array is [1,2,3,4,5], the given min mass is 2.2, the given max mass is 2.8
            var test3 = ms2ScanTest.GetClosestExperimentalIsotopicEnvelopeList(110, 111);
            Assert.AreEqual(test3, null);


            //test normal conditions:look for IsotopicEnvelopes which are in the range of acceptable mass 
            var test4 = ms2ScanTest.GetClosestExperimentalIsotopicEnvelopeList(120, 130);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMinimum(120), 15);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMaximum(130), 23);
            IsotopicEnvelope[] expected4 = ms2ScanTest.ExperimentalFragments.Skip(15).Take(9).ToArray();
            Assert.AreEqual(test4, expected4);

            var test5 = ms2ScanTest.GetClosestExperimentalIsotopicEnvelopeList(400, 500);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMinimum(400), 150);
            Assert.AreEqual(ms2ScanTest.GetClosestFragmentMassMaximum(500), 156);
            IsotopicEnvelope[] expected5 = ms2ScanTest.ExperimentalFragments.Skip(150).Take(7).ToArray();
            Assert.AreEqual(test5, expected5);
        }

    }
}
