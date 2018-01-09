using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using IO.MzML;
using MassSpectrometry;
using UsefulProteomicsDatabases;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class XLTest
    {
        #region Public Methods

        [Test]
        public static void XLTestXlPosCal()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            DigestionParams digestionParams = new DigestionParams
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                MaxMissedCleavages = 2,
                Protease = protease,
                MinPeptideLength = 1
            };
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();

            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();

            var pep = ye[0];
            Assert.AreEqual(pep.BaseSequence, "MNNNK");
            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(CrosslinkerType.DSS);
            Assert.AreEqual(crosslinker.CrosslinkerModSite, 'K');
            Assert.AreEqual(Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass, 128.09496301518999, 1e-9);
            var n = pep.CompactPeptide(TerminusType.None).NTerminalMasses;
            var c = pep.CompactPeptide(TerminusType.None).CTerminalMasses;
            Assert.AreEqual(n.Count(), 4);
            Assert.AreEqual(c.Count(), 4);
            Assert.AreEqual(c[0], 128.09496301518999, 1e-6);
            var x = PsmCross.XlPosCal(pep.CompactPeptide(TerminusType.None), crosslinker).ToArray();
            Assert.AreEqual(x[0], 4);

            var pep2 = ye[2];
            Assert.AreEqual("MNNNKQQQQ", pep2.BaseSequence);
            var n2 = pep2.CompactPeptide(TerminusType.None).NTerminalMasses;
            var c2 = pep2.CompactPeptide(TerminusType.None).CTerminalMasses;
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            Assert.AreEqual(n2[4] - n2[3], 128.09496301518999, 1e-6);
            var x2 = PsmCross.XlPosCal(pep2.CompactPeptide(TerminusType.None), crosslinker).ToArray();
            Assert.AreEqual(x2[0], 4);
        }

        [Test]
        public static void XLTestGenerateIntensityRanks()
        {
            double[] mz = new double[] { 1.0, 1.3, 1.5, 1.7, 1.9, 2.1 };
            double[] intensity = new double[] { 1.1, 1.1, 0.5, 3.2, 0.5, 6.0};
            int[] rank = PsmCross.GenerateIntensityRanks(mz, intensity);
            int[] Rank = new int[] { 4, 3, 6, 2, 5, 1 };
            Assert.AreEqual(rank, Rank);
        }

        [Test]
        public static void XLTestCalculateTotalProductMassesMightHave()
        {
            var CommonParameters = new CommonParameters();
            var proteinList = new List<Protein> { new Protein("CASIQKFGERLCVLHEKTPVSEK", null) };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);

            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Carbamidomethyl of C", "Common Fixed", motif2, TerminusLocalization.Any, 57.02146372068994);
            var variableModifications = new List<ModificationWithMass>() { mod1 };
            var fixedModifications = new List<ModificationWithMass>() { mod2 };
            var localizeableModifications = new List<ModificationWithMass>();

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y };
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

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, lp, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());

            var results = (IndexingResults)engine.Run();

            var digestedList = proteinList[0].Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            foreach (var fdfd in digestedList)
            {
                fdfd.CompactPeptide(TerminusType.None);
                //Assert.Contains(fdfd.CompactPeptide(TerminusType.None), results.PeptideIndex);
            }

            var productMasses = digestedList[3].CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });


            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(CrosslinkerType.DSS);
            var x = PsmCross.XlPosCal(digestedList[3].CompactPeptide(TerminusType.None), crosslinker).ToArray();
            Assert.AreEqual(x[0], 5);

            var myMsDataFile = new XLTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();


            var psmCrossAlpha = new PsmCross(digestedList[3].CompactPeptide(TerminusType.None), 0, 0, i, listOfSortedms2Scans[0]);
            var psmCrossBeta = new PsmCross(digestedList[5].CompactPeptide(TerminusType.None), 0, 0, i, listOfSortedms2Scans[0]);

            var modMassAlpha1 = psmCrossBeta.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass;

            //Another method to calculate modification mass of cross-linked peptides
            //var modMassAlpha2 = listOfSortedms2Scans[0].PrecursorMass - psmCrossAlpha.compactPeptide.MonoisotopicMassIncludingFixedMods;

            var productMassesAlphaList = PsmCross.XLCalculateTotalProductMasses(psmCrossAlpha, modMassAlpha1, crosslinker, lp, true, false);

            Assert.AreEqual(productMassesAlphaList[0].ProductMz.Length, 35);
            Assert.AreEqual(productMassesAlphaList[0].ProductMz[26], 2312.21985342336);
        }


        #endregion Public Methods
    }

    internal class XLTestDataFile : MsDataFile<IMzmlScan>
    {
        #region Public Constructors

        public XLTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 50, 60, 70, 80, 90, 2871.45827313843.ToMz(3) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };

            var mz2 = new double[] { 50, 60, 70, 76.0393, 133.0608, 147.0764, 190.0822, 247.1037, 257.1244, 258.127, 275.1350, 385.1830, 442.2045, 630.27216, 900 };
            var intensities2 = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 5.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 100000, 2871.45827313843.ToMz(3), 3, 1, 2871.45827313843.ToMz(3), 2, DissociationType.HCD, 1, 2871.45827313843.ToMz(3), 1, "scan=2"));

            Scans = ScansHere.ToArray();
        }

        #endregion Public Constructors

        #region Public Properties

        public string FilePath
        {
            get
            {
                return "XLTestDataFile";
            }
        }

        public string Name
        {
            get
            {
                return "XLTestDataFile";
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

        #endregion Public Methods
    }
}