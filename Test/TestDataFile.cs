using Chemistry;
using EngineLayer;
using IO.MzML;
using MassSpectrometry;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    internal class TestDataFile : IMsDataFile<IMzmlScan>
    {

        #region Private Fields

        private readonly List<IMzmlScan> Scans;

        #endregion Private Fields

        #region Public Constructors

        public TestDataFile()
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 402.18629720155.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            Scans = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1) };
            var mz2 = new double[] { 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            Scans.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType.HCD, 1, 402.18629720155.ToMz(2), 1));
        }

        public TestDataFile(List<PeptideWithSetModifications> pepWithSetModss)
        {
            Scans = new List<IMzmlScan>();
            for (int i = 0; i < pepWithSetModss.Count; i++)
            {
                var pepWithSetMods = pepWithSetModss[i];
                var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) };
                var intensities1 = new double[] { 1, 1, 1 };
                var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

                Scans.Add(new MzmlScan(2 * i + 1, MassSpectrum1, 1, true, Polarity.Positive, 2 * i, new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Orbitrap, 1000, 1));

                List<double> mz2 = new List<double>();
                List<double> intensities2 = new List<double>();
                foreach (var aok in pepWithSetMods.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                {
                    mz2.Add(aok.ToMz(1));
                    mz2.Add((aok + 1.003).ToMz(1));
                    intensities2.Add(1);
                    intensities2.Add(1);
                }
                var MassSpectrum2 = new MzmlMzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
                Scans.Add(new MzmlScanWithPrecursor(2 * i + 2, MassSpectrum2, 2, true, Polarity.Positive, 2 * i + 1, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 2 * i + 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 1));
            }
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods)
        {
            var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            Scans = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1) };

            List<double> mz2 = new List<double>();
            List<double> intensities2 = new List<double>();
            foreach (var aok in pepWithSetMods.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
            {
                mz2.Add(aok.ToMz(1));
                mz2.Add((aok + 1.003).ToMz(1));
                intensities2.Add(1);
                intensities2.Add(1);
            }
            var MassSpectrum2 = new MzmlMzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
            var scan2 = new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 1);
            scan2.ComputeSelectedPeakIntensity(MassSpectrum1);
            scan2.ComputeMonoisotopicPeakIntensity(MassSpectrum1);
            Scans.Add(scan2);
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods, string v)
        {
            if (v.Equals("quadratic"))
            {
                // Add three ms1 peaks with charge 2, exact
                var MassSpectrum1 = new MzmlMzSpectrum(new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) }, new double[] { 1, 1, 1 }, false);

                List<double> mz2 = new List<double>();
                List<double> intensities2 = new List<double>();
                foreach (var aok in pepWithSetMods.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                {
                    var t1 = aok.ToMz(1);
                    var c = 0.0000001;
                    mz2.Add(t1 + c * Math.Pow(t1, 2));

                    Console.WriteLine("orig: " + aok.ToMz(1) + " new: " + (t1 + c * Math.Pow(t1, 2)));

                    var t2 = (aok + 1.003).ToMz(1);
                    mz2.Add(t2 + c * Math.Pow(t2, 2));
                    intensities2.Add(1);
                    intensities2.Add(1);
                }
                var MassSpectrum2 = new MzmlMzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);

                var scan2 = new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 1);
                scan2.ComputeSelectedPeakIntensity(MassSpectrum1);
                scan2.ComputeMonoisotopicPeakIntensity(MassSpectrum1);
                Scans = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000,1),
               scan2 };
            }
        }

        #endregion Public Constructors

        #region Public Properties

        public string FilePath
        {
            get
            {
                return "TestDataFile";
            }
        }

        public string Name
        {
            get
            {
                return "TestDataFile";
            }
        }

        public int NumSpectra
        {
            get
            {
                return Scans.Count;
            }
        }

        #endregion Public Properties

        //public IEnumerator<IMsDataScan<IMzSpectrum<IMzPeak>, IMzPeak>> GetEnumerator()
        //{
        //    return Scans.GetEnumerator();
        //}

        //public IMsDataScan<IMzSpectrum<IMzPeak>, IMzPeak> GetOneBasedScan(int oneBasedScanNumber)
        //{
        //    return Scans[oneBasedScanNumber - 1];
        //}

        //IEnumerator IEnumerable.GetEnumerator()
        //{
        //    return Scans.GetEnumerator();
        //}

        //IMsDataScan<DefaultMzSpectrum> IMsDataFile<DefaultMzSpectrum>.GetOneBasedScan(int oneBasedScanNumber)
        //{
        //    return Scans[oneBasedScanNumber - 1];
        //}

        //IEnumerator<IMsDataScan<DefaultMzSpectrum>> IEnumerable<IMsDataScan<DefaultMzSpectrum>>.GetEnumerator()
        //{
        //    return Scans.GetEnumerator();
        //}

        //public override int GetClosestOneBasedSpectrumNumber(double retentionTime)
        //{
        //    throw new NotImplementedException();
        //}

        //public override void Open()
        //{
        //    throw new NotImplementedException();
        //}

        //public override void Close()
        //{
        //    throw new NotImplementedException();
        //}

        //protected override MzmlScan GetMsDataOneBasedScanFromFile(int oneBasedSpectrumNumber)
        //{
        //    throw new NotImplementedException();
        //}

        //protected override int GetNumSpectra()
        //{
        //    throw new NotImplementedException();
        //}

        //protected override TestScan GetMsDataOneBasedScanFromFile(int oneBasedSpectrumNumber)
        //{
        //    throw new NotImplementedException();
        //}

        #region Public Methods

        public IMzmlScan GetOneBasedScan(int oneBasedScanNumber)
        {
            return Scans[oneBasedScanNumber - 1];
        }

        public IEnumerable<TestScan> GetMsScansInTimeRange(double firstRT, double lastRT)
        {
            throw new NotImplementedException();
        }

        public void LoadAllScansInMemory()
        {
        }

        public IEnumerator<IMzmlScan> GetEnumerator()
        {
            return Scans.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return Scans.GetEnumerator();
        }

        public int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            throw new NotImplementedException();
        }

        public void Close()
        {
            throw new NotImplementedException();
        }

        public void Open()
        {
            throw new NotImplementedException();
        }

        IEnumerable<IMzmlScan> IMsDataFile<IMzmlScan>.GetMsScansInTimeRange(double firstRT, double lastRT)
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods

    }
}