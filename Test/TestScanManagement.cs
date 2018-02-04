using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;

using NUnit.Framework;
using TaskLayer;


namespace Test
{
    [TestFixture]
    class TestScanManagement
    {
        #region Public Method

        [Test]
        public static void TestGetCombinedMs2Scans()
        {
            var myMsDataFile = new TestDataFile(5);
            
            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Dictionary<int, double> listOfScanPrecusor = new Dictionary<int, double>();

            List<Ms2ScanWithSpecificMass> test = new List<Ms2ScanWithSpecificMass>();

            foreach (var ms2scan in myMsDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
            {
                if (listOfScanPrecusor.Contains(new KeyValuePair<int, double>(ms2scan.OneBasedScanNumber, ms2scan.SelectedIonMZ)))
                {
                    break;
                }
                else
                {
                    listOfScanPrecusor.Add(ms2scan.OneBasedScanNumber, ms2scan.SelectedIonMZ);
                    List<int> currentScanMS2OneBasedScanNumber = new List<int>();
                    if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                    {
                        for (int i = 1; i < 7; i++)
                        {
                            if (ms2scan.OneBasedScanNumber + i <= myMsDataFile.NumSpectra)
                            {
                                
                                var x = myMsDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().ElementAt(i);
                                currentScanMS2OneBasedScanNumber.Add(x.OneBasedScanNumber);
                                if (x.MsnOrder == 2 && x.SelectedIonMZ == ms2scan.SelectedIonMZ)
                                {
                                    ms2scan.MassSpectrum.XArray.Concat(x.MassSpectrum.XArray);
                                    ms2scan.MassSpectrum.YArray.Concat(x.MassSpectrum.YArray);
                                }
                                if (x.MsnOrder == 3 && currentScanMS2OneBasedScanNumber.Contains(x.OneBasedPrecursorScanNumber.Value))
                                {
                                    ms2scan.MassSpectrum.XArray.Concat(x.MassSpectrum.XArray);
                                    ms2scan.MassSpectrum.YArray.Concat(x.MassSpectrum.YArray);
                                }
                            }

                        }
                    }
                    test.Add(new Ms2ScanWithSpecificMass(ms2scan, ms2scan.SelectedIonMonoisotopicGuessMz.Value, ms2scan.SelectedIonChargeStateGuess.Value, ""));

                }
                
            }

            Assert.AreEqual(5, myMsDataFile.NumSpectra);
        }

        #endregion 
    }
}
