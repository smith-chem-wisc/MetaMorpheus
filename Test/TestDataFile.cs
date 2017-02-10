using Chemistry;
using EngineLayer;
using IO.MzML;
using MassSpectrometry;

using Spectra;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    internal class TestDataFile : IMsDataFile<TestScan>
    {

        #region Private Fields

        private readonly List<TestScan> Scans;

        #endregion Private Fields

        #region Public Constructors

        public TestDataFile()
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 402.18629720155.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            Scans = new List<TestScan> { new TestScan(1, 1, MassSpectrum1, 1) };
            var mz2 = new double[] { 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            Scans.Add(new TestScanWithPrecursor(2, 2, MassSpectrum2, 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, 1));
        }

        public TestDataFile(List<PeptideWithSetModifications> pepWithSetModss)
        {
            Scans = new List<TestScan>();
            for (int i = 0; i < pepWithSetModss.Count; i++)
            {
                var pepWithSetMods = pepWithSetModss[i];
                var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) };
                var intensities1 = new double[] { 1, 1, 1 };
                var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

                Scans.Add(new TestScan(2 * i + 1, 2 * i, MassSpectrum1, 2 * i));

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
                Scans.Add(new TestScanWithPrecursor(2 * i + 2, 2 * i + 1, MassSpectrum2, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2 * i + 1, 2 * i + 1));
            }
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods)
        {
            var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            Scans = new List<TestScan> { new TestScan(1, 1, MassSpectrum1, 1) };

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
            Scans.Add(new TestScanWithPrecursor(2, 2, MassSpectrum2, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1));
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

        #region Public Methods


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

        public TestScan GetOneBasedScan(int oneBasedScanNumber)
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

        public IEnumerator<TestScan> GetEnumerator()
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

        #endregion Public Methods

    }
}