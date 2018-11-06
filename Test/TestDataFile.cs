using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    internal class TestDataFile : MsDataFile
    {
        public TestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 402.18629720155.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            var mz2 = new double[] { 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=2", 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType.HCD, 1, 402.18629720155.ToMz(2)));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(double closeMassDifference) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 402.18629720155.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };
            var mz2 = new double[] { 50, 60, 70, 147.0764, 258.132 - closeMassDifference - Constants.ProtonMass, 258.132 - Constants.ProtonMass, 275.1350 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=2", 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType.HCD, 1, 402.18629720155.ToMz(2)));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(bool emptyScan) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50 };
            var intensities1 = new double[] { 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };
            var mz2 = new double[] { 1 };
            var intensities2 = new double[] { 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=2", 402.18629720155.ToMz(2), 2, 1, 402.18629720155.ToMz(2), 2, DissociationType.HCD, 1, 402.18629720155.ToMz(2)));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(string slightlyLargerDataFile) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 630.27216.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };
            var mz2 = new double[] { 50, 60, 70, 76.0393, 133.0608, 147.0764, 190.0822, 247.1037, 257.1244, 258.127, 275.1350, 385.1830, 442.2045, 630.27216 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=2", 630.27216.ToMz(2), 2, 1, 630.27216.ToMz(2), 2, DissociationType.HCD, 1, 630.27216.ToMz(2)));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(double precursor, double[] products) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { precursor.ToMz(2) };
            var intensities1 = new double[] { 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };
            var mz2 = products;
            var intensities2 = new double[products.Length];
            for (int i = 0; i < intensities2.Length; i++)
            {
                intensities2[i] = 1;
            }
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=2", precursor.ToMz(2), 2, 1, precursor.ToMz(2), 2, DissociationType.HCD, 1, precursor.ToMz(2)));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(List<PeptideWithSetModifications> pepWithSetModss, bool additionalMasses = false) 
            : base(pepWithSetModss.Count * 2, new SourceFile(@"no nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
        {
            var ScansHere = new List<MsDataScan>();
            for (int i = 0; i < pepWithSetModss.Count; i++)
            {
                var pepWithSetMods = pepWithSetModss[i];
                var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(3), (pepWithSetMods.MonoisotopicMass + Constants.C13MinusC12).ToMz(3), (pepWithSetMods.MonoisotopicMass + 2 * Constants.C13MinusC12).ToMz(3), pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + Constants.C13MinusC12).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2 * Constants.C13MinusC12).ToMz(2) };
                var intensities1 = new double[] { 1, 0.5, 0.25, 1, 0.5, 0.25 };
                var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

                ScansHere.Add(new MsDataScan(MassSpectrum1, 2 * i + 1, 1, true, Polarity.Positive, 2 * i, new MzLibUtil.MzRange(0, 10000), "gg", MZAnalyzerType.Orbitrap, 1000, 1, null, "scan=1"));

                List<double> mz2 = new List<double>();
                List<double> intensities2 = new List<double>();
                IEnumerable<double> additionalMassesArray;
                if (additionalMasses)
                    additionalMassesArray = new List<double> { 260.08307817722, 397.14199003569, 498.18966850487, 612.23259594625, 683.2697097314, 146.10552769922, 217.14264148437 };
                else
                    additionalMassesArray = new List<double>();
                foreach (var aok in pepWithSetMods.Fragment(DissociationType.HCD, FragmentationTerminus.Both))
                {
                    mz2.Add(aok.NeutralMass.ToMz(1));
                    mz2.Add((aok.NeutralMass + Constants.C13MinusC12).ToMz(1));
                    intensities2.Add(1);
                    intensities2.Add(1);
                }
                var MassSpectrum2 = new MzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
                ScansHere.Add(new MsDataScan(MassSpectrum2, 2 * i + 2, 2, true, Polarity.Positive, 2 * i + 1, new MzLibUtil.MzRange(0, 10000), "gg", MZAnalyzerType.Orbitrap, 234734, 1, null, "scan=2", pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 2 * i + 1, pepWithSetMods.MonoisotopicMass.ToMz(2)));
            }
            Scans = ScansHere.ToArray();
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods) 
            : base(2, new SourceFile(@"no nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
        {
            var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + Constants.C13MinusC12).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2 * Constants.C13MinusC12).ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            List<double> mz2 = new List<double>();
            List<double> intensities2 = new List<double>();
            foreach (var aok in pepWithSetMods.Fragment(DissociationType.HCD, FragmentationTerminus.Both))
            {
                mz2.Add(aok.NeutralMass.ToMz(1));
                mz2.Add((aok.NeutralMass + Constants.C13MinusC12).ToMz(1));
                intensities2.Add(1);
                intensities2.Add(1);
            }
            var MassSpectrum2 = new MzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
            var scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, 1, null, "scan=2", pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2));
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
                var MassSpectrum1 = new MzSpectrum(new double[] { pepWithSetMods.MonoisotopicMass.ToMz(2), (pepWithSetMods.MonoisotopicMass + Constants.C13MinusC12).ToMz(2), (pepWithSetMods.MonoisotopicMass + 2 * Constants.C13MinusC12).ToMz(2) }, new double[] { 1, 1, 1 }, false);

                List<double> mz2 = new List<double>();
                List<double> intensities2 = new List<double>();
                foreach (var aok in pepWithSetMods.Fragment(DissociationType.HCD, FragmentationTerminus.Both))
                {
                    var t1 = aok.NeutralMass.ToMz(1);
                    var c = 0.0000001;
                    mz2.Add(t1 + c * Math.Pow(t1, 2));
                    var t2 = (aok.NeutralMass + Constants.C13MinusC12).ToMz(1);
                    mz2.Add(t2 + c * Math.Pow(t2, 2));
                    intensities2.Add(1);
                    intensities2.Add(1);
                }
                var MassSpectrum2 = new MzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);

                var scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, 1, null, "scan=2", pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2));
                scan2.ComputeSelectedPeakIntensity(MassSpectrum1);
                scan2.ComputeMonoisotopicPeakIntensity(MassSpectrum1);
                var ScansHere = new List<MsDataScan>  { new MsDataScan( MassSpectrum1,1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000,1, null, "scan=1"),
               scan2 };
                Scans = ScansHere.ToArray();
            }
        }

        public TestDataFile(PeptideWithSetModifications pepWithSetMods, int charge, double intensity, double rt) : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { pepWithSetMods.MonoisotopicMass.ToMz(charge), (pepWithSetMods.MonoisotopicMass + 1.003).ToMz(charge), (pepWithSetMods.MonoisotopicMass + 2.005).ToMz(charge) };
            var intensities1 = new double[] { intensity, intensity * 10, intensity / 10 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);

            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, rt, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            List<double> mz2 = new List<double>();
            List<double> intensities2 = new List<double>();
            foreach (var aok in pepWithSetMods.Fragment(DissociationType.HCD,FragmentationTerminus.Both))
            {
                mz2.Add(aok.NeutralMass.ToMz(1));
                mz2.Add((aok.NeutralMass + 1.003).ToMz(1));
                intensities2.Add(intensity);
                intensities2.Add(intensity);
            }
            var MassSpectrum2 = new MzSpectrum(mz2.OrderBy(b => b).ToArray(), intensities2.ToArray(), false);
            var scan2 = new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, rt + 0.01, new MzLibUtil.MzRange(0, 10000), "df", MZAnalyzerType.Orbitrap, 234734, 1, null, "scan=2", pepWithSetMods.MonoisotopicMass.ToMz(2), 2, 1, pepWithSetMods.MonoisotopicMass.ToMz(2), 2, DissociationType.HCD, 1, pepWithSetMods.MonoisotopicMass.ToMz(2));
            scan2.ComputeSelectedPeakIntensity(MassSpectrum1);
            scan2.ComputeMonoisotopicPeakIntensity(MassSpectrum1);
            ScansHere.Add(scan2);
            Scans = ScansHere.ToArray();
        }

        public TestDataFile(int MS3 = 5) : base(MS3, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 764.1376.ToMz(2) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            var mz2 = new double[] { 52, 62, 72, 147.0764, 257.1244, 258.127, 275.1350, 502 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=2", 764.1376.ToMz(2), 2, 1, 764.1376.ToMz(2), 2, DissociationType.CID, 1, 764.1376.ToMz(1)));

            var mz3 = new double[] { 53, 63, 73, 148.0764, 258.1244, 259.127, 276.1350, 503 };
            var intensities3 = new double[] { 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum3 = new MzSpectrum(mz3, intensities3, false);
            ScansHere.Add(new MsDataScan(MassSpectrum3, 3, 2, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=3", 764.1376.ToMz(2), 2, 1, 764.1376.ToMz(2), 2, DissociationType.ETD, 1, 764.1376.ToMz(1)));

            var mz4 = new double[] { 54, 64, 74, 149.0764, 259.1244, 260.127, 277.1350, 504 };
            var intensities4 = new double[] { 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum4 = new MzSpectrum(mz4, intensities4, false);
            ScansHere.Add(new MsDataScan(MassSpectrum4, 4, 3, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=4", 275.1350.ToMz(1), 1, 1, 275.1350.ToMz(1), 1, DissociationType.HCD, 2, 275.1350.ToMz(1)));

            var mz5 = new double[] { 55, 65, 75, 150.0764, 260.1244, 261.127, 278.1350, 505 };
            var intensities5 = new double[] { 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum5 = new MzSpectrum(mz5, intensities5, false);
            ScansHere.Add(new MsDataScan(MassSpectrum5, 5, 3, true, Polarity.Positive, 2,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 1, null, "scan=5", 257.1244.ToMz(1), 1, 1, 257.1244.ToMz(1), 1, DissociationType.HCD, 2, 257.1244.ToMz(1)));

            Scans = ScansHere.ToArray();
        }

        public TestDataFile(double[] ms2Mz, double[] ms2Intensities, double precursorMass, int precursorZ, double rt = 1.0) : base(2, new SourceFile(null, null, null, null, null))
        {
            var ms1 = new MzSpectrum(new double[] { precursorMass.ToMz(precursorZ), (precursorMass + 1.003).ToMz(precursorZ) }, new double[] { 1, 1}, false);
            var ms2 = new MzSpectrum(ms2Mz, ms2Intensities, false);

            var ScansHere = new List<MsDataScan> {
                new MsDataScan(ms1, 1, 1, true, Polarity.Positive, rt, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1"),
                new MsDataScan(ms2, 1, 2, true, Polarity.Positive, rt + 0.01, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=2", precursorMass.ToMz(precursorZ), precursorZ, 1, precursorMass.ToMz(precursorZ), 1.0, DissociationType.HCD, 1, precursorMass.ToMz(precursorZ)) };
            
            Scans = ScansHere.ToArray();
        }

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

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, null, Scans[0].NativeId);
        }

        public override MsDataScan GetOneBasedScan(int scanNumber)
        {
            return Scans[scanNumber - 1];
        }

        public override IEnumerable<MsDataScan> GetMS1Scans()
        {
            throw new NotImplementedException();
        }
    }
}