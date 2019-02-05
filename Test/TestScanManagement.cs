using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    internal class TestScanManagement
    {
        //[Test]
        public static void TestGetCombinedMs2Scans()
        {
            var myMsDataFile = new TestDataFile(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //Write prime code to combine MS2MS3
            Dictionary<int, double> listOfScanPrecusor = new Dictionary<int, double>();

            List<MsDataScan> ListOfSortedMsScans = new List<MsDataScan>();
            List<Ms2ScanWithSpecificMass> test = new List<Ms2ScanWithSpecificMass>();

            foreach (var ms2scan in myMsDataFile.GetAllScansList().Where(x => x.MsnOrder != 1))
            {
                if (ms2scan.MsnOrder == 2 && !listOfScanPrecusor.Contains(new KeyValuePair<int, double>(ms2scan.OneBasedPrecursorScanNumber.Value, ms2scan.SelectedIonMZ.Value)))
                {
                    if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                    {
                        listOfScanPrecusor.Add(ms2scan.OneBasedPrecursorScanNumber.Value, ms2scan.SelectedIonMZ.Value);
                        List<int> currentScanMS2OneBasedScanNumber = new List<int>
                        {
                            ms2scan.OneBasedScanNumber
                        };
                        var mz2 = ms2scan.MassSpectrum.XArray.ToList();
                        var intensities2 = ms2scan.MassSpectrum.YArray.ToList();
                        for (int i = 1; i < 7; i++)
                        {
                            if (ms2scan.OneBasedScanNumber + i <= myMsDataFile.NumSpectra)
                            {
                                var x = myMsDataFile.GetOneBasedScan(ms2scan.OneBasedScanNumber + i);
                                //var x = myMsDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().ElementAt(i);

                                if (x.MsnOrder == 2 && x.SelectedIonMZ == ms2scan.SelectedIonMZ)
                                {
                                    currentScanMS2OneBasedScanNumber.Add(x.OneBasedScanNumber);
                                    mz2.AddRange(x.MassSpectrum.XArray.ToList());
                                    intensities2.AddRange(x.MassSpectrum.YArray.ToList());
                                }
                                if (x.MsnOrder == 3 && currentScanMS2OneBasedScanNumber.Contains(x.OneBasedPrecursorScanNumber.Value))
                                {
                                    mz2.AddRange(x.MassSpectrum.XArray.ToList());
                                    intensities2.AddRange(x.MassSpectrum.YArray.ToList());
                                }
                            }
                        }
                        var MassSpectrum2 = new MzSpectrum(mz2.ToArray(), intensities2.ToArray(), false);
                        ListOfSortedMsScans.Add(new MsDataScan(MassSpectrum2, ms2scan.OneBasedScanNumber, ms2scan.MsnOrder, ms2scan.IsCentroid, Polarity.Positive, ms2scan.RetentionTime,
                            ms2scan.ScanWindowRange, ms2scan.ScanFilter, ms2scan.MzAnalyzer, ms2scan.TotalIonCurrent, ms2scan.InjectionTime, null, "", ms2scan.SelectedIonMZ, ms2scan.SelectedIonChargeStateGuess, ms2scan.SelectedIonIntensity, ms2scan.IsolationMz, null, ms2scan.DissociationType, ms2scan.OneBasedPrecursorScanNumber, ms2scan.SelectedIonMonoisotopicGuessMz));
                    }
                }
            }
            foreach (var ms2scan in ListOfSortedMsScans.Where(x => x.MsnOrder != 1))
            {
                test.Add(new Ms2ScanWithSpecificMass(ms2scan, ms2scan.SelectedIonMonoisotopicGuessMz.Value, ms2scan.SelectedIonChargeStateGuess.Value, "", new CommonParameters()));
            }
            var testToArray = test.OrderBy(b => b.PrecursorMass).ToArray();

        }
    }
}