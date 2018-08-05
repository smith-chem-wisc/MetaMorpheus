using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkAnalysis;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using OxyPlot;
using OxyPlot.Pdf;

namespace Test
{
    [TestFixture]
    public static class GlyTest
    {
        [Test]
        public static void GlyTest_FragmentIons()
        {
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, cIons: true, zDotIons: true, scoreCutoff: 2, digestionParams: new DigestionParams(minPeptideLength: 5));
            var xlSearchParameters = new XlSearchParameters { XlCharge_2_3_PrimeFragment = false };

            //Create databases contain protein.
            var proteinList = new List<Protein> { new Protein("DANNTQFQFTSR", "25170") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Carbamidomethyl of C", "Common Fixed", motif2, TerminusLocalization.Any, 57.02146372068994);
            var variableModifications = new List<ModificationWithMass>() { mod1 };
            var fixedModifications = new List<ModificationWithMass>() { mod2 };
            var localizeableModifications = new List<ModificationWithMass>();

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y};
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            //Generate digested peptide lists.
            List<PeptideWithSetModifications> digestedList = new List<PeptideWithSetModifications>();
            foreach (var item in proteinList)
            {
                var digested = item.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                digestedList.AddRange(digested);
            }

            foreach (var fdfd in digestedList)
            {
                fdfd.CompactPeptide(TerminusType.None);
            }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, lp, 1, DecoyType.Reverse, new List<DigestionParams>
            { commonParameters.DigestionParams }, commonParameters, 30000, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var fragmentIndexCount = indexResults.FragmentIndex.Count(p => p != null);
            var fragmentIndexAll = indexResults.FragmentIndex.Select((s, j) => new { j, s }).Where(p => p.s != null).Select(t => t.j).ToList();
            Assert.IsTrue(fragmentIndexAll.Count() > 0);

            //Get MS2 scans.
            var myMsDataFile = new GlyTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).ToArray();
            
            //TwoPassCrosslinkSearchEngine.Run().
            List<PsmCross> newPsms = new List<PsmCross>();
            new TwoPassCrosslinkSearchEngine(newPsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, lp, 0, commonParameters, false, true, xlSearchParameters.XlPrecusorMsTl, null, xlSearchParameters.CrosslinkSearchTop, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, xlSearchParameters.XlCharge_2_3, xlSearchParameters.XlCharge_2_3_PrimeFragment, new List<string> { }).Run();

            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            new CrosslinkAnalysisEngine(newPsms, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, lp, null, null, TerminusType.None, commonParameters, new List<string> { }).Run();


        }

        [Test]
        public static void GlyTest_BinarySearch()
        {
            double[] array = new double[] { 3.44, 3.45, 4.55, 4.55, 4.55, 4.55, 5.66, 5.66, 6.77, 6.77 };
            double x = 3.55;
            double y = 4.55;
            double z = 5.44;
            var xid = Array.BinarySearch(array, x);
            if (xid<0){ xid = ~xid;}          
            var yid = Array.BinarySearch(array, y);
            if (yid < 0) { yid = ~yid; }
            else{ while (yid < array.Length && array[yid+1]==array[yid]){yid++;} }
            //else{ while (yid > 0 && array[yid - 1] == array[yid]) { yid--; } }
            var zid = Array.BinarySearch(array, z);
            if (zid < 0) { zid = ~zid -1; }
        }
    }

    internal class GlyTestDataFile : MsDataFile
    {
        public GlyTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/25170.mgf");
            double[] mz;
            double[] intensity;
            GetScan(raw, out mz, out intensity);

            var mz1 = new double[] { 2644.06992.ToMz(2) };
            var intensities1 = new double[] { 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Orbitrap, 1000, 1, null, "scan=1") };

            var MassSpectrum2 = new MzSpectrum(mz, intensity, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Orbitrap, intensity.Sum(), 1.0, null, "scan=2", 2644.06992.ToMz(2),
                2, 1, 2644.06992.ToMz(2), 2, DissociationType.HCD, 1, 2644.06992.ToMz(2)));

            Scans = ScansHere.ToArray();
        }

        public string FilePath
        {
            get
            {
                return "GlyTestDataFile";
            }
        }

        public string Name
        {
            get
            {
                return "GlyTestDataFile";
            }
        }

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, null, Scans[0].NativeId);
        }

        private void GetScan(string filePath, out double[] mz, out double[]intensity)
        {
            List<double> mzList = new List<double>();
            List<double> intensityList = new List<double>();
            
            using (StreamReader glycans = new StreamReader(filePath))
            {            
                while (glycans.Peek() != -1)
                {
                    string line = glycans.ReadLine();
                    if (!line.StartsWith("B") && !line.StartsWith("T") && !line.StartsWith("R") && !line.StartsWith("P") && !line.StartsWith("T") && !line.StartsWith("C") && !line.StartsWith("E"))
                    {       
                        mzList.Add(Convert.ToDouble(line.Split(null)[0]));
                        intensityList.Add(Convert.ToDouble(line.Split(null)[1]));
                    }
                }
            }
            mz = mzList.ToArray();
            intensity = intensityList.ToArray();
        }
    }
}
