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

using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using MzLibUtil;
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

            Dictionary<int, double> listOfScanPrecusor;

            foreach (var ms2scan in myMsDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
            {
                if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                {                    
                    for (int i = 1; i < 10; i++)
                    {
                        var x = myMsDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().ElementAt( ms2scan.OneBasedScanNumber + i);
                        if (x.MsnOrder != 1)
                        {

                        }
                    }
                }
            }
            Assert.AreEqual(5, myMsDataFile.NumSpectra);
        }

        #endregion 
    }
}
