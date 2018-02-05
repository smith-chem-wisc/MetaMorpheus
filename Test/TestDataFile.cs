using Chemistry;
using EngineLayer;
using IO.MzML;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    internal class TestDataFile : MsDataFile<IMzmlScan>
    {
        #region Public Constructors

        public TestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 402.18629720155.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };
            var mz2 = new double[] { 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType.HCD, 1, 402.18629720155.ToMz(2), 1, "scan=2"));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(double closeMassDifference) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 402.18629720155.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };
            var mz2 = new double[] { 50, 60, 70, 147.0764, 258.132 - closeMassDifference - Constants.protonMass, 258.132 - Constants.protonMass, 275.1350 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType.HCD, 1, 402.18629720155.ToMz(2), 1, "scan=2"));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(bool emptyScan) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50 };
            var intensities1 = new double[] { 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };
            var mz2 = new double[] { 1 };
            var intensities2 = new double[] { 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType.HCD, 1, 402.18629720155.ToMz(2), 1, "scan=2"));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(string slightlyLargerDataFile) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 630.27216.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };
            var mz2 = new double[] { 50, 60, 70, 76.0393, 133.0608, 147.0764, 190.0822, 247.1037, 257.1244, 258.127, 275.1350, 385.1830, 442.2045, 630.27216 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 630.27216.ToMz(2), 2, 1, 630.27216.ToMz(2), 2, DissociationType.HCD, 1, 630.27216.ToMz(2), 1, "scan=2"));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(List<PeptideWithSetModifications> pepWithSetModss, bool additionalMasses = false) : base(pepWithSetModss.Count * 2, new SourceFile(@"no nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
        {
            var ScansHere = new List<IMzmlScan>();
            for (int i = 0; i < pepWithSetModss.Count; i++)
            {
                var pepWithSetMods = pepWithSetModss[i];
                var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(3), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(3), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(3), pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) };
                var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
                var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

                ScansHere.Add(new MzmlScan(2 * i + 1, MassSpectrum1, 1, true, Polarity.Positive, 2 * i, new MzLibUtil.MzRange(0, 10000), "gg", MZAnalyzerType.Orbitrap, 1000, 1, "scan=1"));

                List<double> mz2 = new List<double>();
                List<double> intensities2 = new List<double>();
                IEnumerable<double> additionalMassesArray;
                if (additionalMasses)
                    additionalMassesArray = new List<double> { 260.08307817722, 397.14199003569, 498.18966850487, 612.23259594625, 683.2697097314, 146.10552769922, 217.14264148437 };
                else
                    additionalMassesArray = new List<double>();
                foreach (var aok in pepWithSetMods.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }).Concat(additionalMassesArray))
                {
                    mz2.Add(aok.ToMz(1));
                    mz2.Add((aok + 1.003).ToMz(1));
                    intensities2.Add(1);
                    intensities2.Add(1);
                }
                var MassSpectrum2 = new MzmlMzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
                ScansHere.Add(new MzmlScanWithPrecursor(2 * i + 2, MassSpectrum2, 2, true, Polarity.Positive, 2 * i + 1, new MzLibUtil.MzRange(0, 10000), "gg", MZAnalyzerType.Orbitrap, 234734, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 2 * i + 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 1, "scan=2"));
            }
            Scans = ScansHere.ToArray();
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods) : base(2, new SourceFile(@"no nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
        {
            var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };

            List<double> mz2 = new List<double>();
            List<double> intensities2 = new List<double>();
            foreach (var aok in pepWithSetMods.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
            {
                mz2.Add(aok.ToMz(1));
                mz2.Add((aok + 1.003).ToMz(1));
                intensities2.Add(1);
                intensities2.Add(1);
            }
            var MassSpectrum2 = new MzmlMzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
            var scan2 = new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 1, "scan=2");
            scan2.ComputeSelectedPeakIntensity(MassSpectrum1);
            scan2.ComputeMonoisotopicPeakIntensity(MassSpectrum1);
            ScansHere.Add(scan2);
            Scans = ScansHere.ToArray();
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods, string v) : base(2, new SourceFile(null, null, null, null, null))
        {
            if (v.Equals("quadratic"))
            {
                // Add three ms1 peaks with charge 2, exact
                var MassSpectrum1 = new MzmlMzSpectrum(new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(2) }, new double[] { 1, 1, 1 }, false);

                List<double> mz2 = new List<double>();
                List<double> intensities2 = new List<double>();
                foreach (var aok in pepWithSetMods.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
                {
                    var t1 = aok.ToMz(1);
                    var c = 0.0000001;
                    mz2.Add(t1 + c * Math.Pow(t1, 2));
                    var t2 = (aok + 1.003).ToMz(1);
                    mz2.Add(t2 + c * Math.Pow(t2, 2));
                    intensities2.Add(1);
                    intensities2.Add(1);
                }
                var MassSpectrum2 = new MzmlMzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);

                var scan2 = new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 2, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 1, "scan=2");
                scan2.ComputeSelectedPeakIntensity(MassSpectrum1);
                scan2.ComputeMonoisotopicPeakIntensity(MassSpectrum1);
                var ScansHere = new List<IMzmlScan>  { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000,1, "scan=1"),
               scan2 };
                Scans = ScansHere.ToArray();
            }
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods, int charge, double intensity, double rt) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(charge), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(charge), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(charge) };
            var intensities1 = new double[] { intensity, intensity * 10, intensity / 10 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, rt, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };

            List<double> mz2 = new List<double>();
            List<double> intensities2 = new List<double>();
            foreach (var aok in pepWithSetMods.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
            {
                mz2.Add(aok.ToMz(1));
                mz2.Add((aok + 1.003).ToMz(1));
                intensities2.Add(intensity);
                intensities2.Add(intensity);
            }
            var MassSpectrum2 = new MzmlMzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
            var scan2 = new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, rt + 0.01, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 1, "scan=2");
            scan2.ComputeSelectedPeakIntensity(MassSpectrum1);
            scan2.ComputeMonoisotopicPeakIntensity(MassSpectrum1);
            ScansHere.Add(scan2);
            Scans = ScansHere.ToArray();
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

        #endregion Public Properties

        #region Public Methods

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            Scans[0] = new MzmlScan(Scans[0].OneBasedScanNumber, massSpectrum, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, Scans[0].NativeId);
        }

        public override IMzmlScan GetOneBasedScan(int scanNumber)
        {
            return Scans[scanNumber - 1];
        }

        public override IEnumerable<IMzmlScan> GetMS1Scans()
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods
    }
}