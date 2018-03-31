using Chemistry;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class FdrTest
    {
        #region Public Methods

        [Test]
        public static void FdrTestMethod()
        {
            MassDiffAcceptor searchModes = new DotMassDiffAcceptor(null, new List<double> { 0, 1.0029 }, new PpmTolerance(5));
            List<string> nestedIds = new List<string>();

            Protein p = new Protein("MNKNNKNNNKNNNNK", null);
            DigestionParams digestionParams = new DigestionParams();
            var digested = p.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();

            PeptideWithSetModifications pep1 = digested[0];
            PeptideWithSetModifications pep2 = digested[1];
            PeptideWithSetModifications pep3 = digested[2];
            PeptideWithSetModifications pep4 = digested[3];

            TestDataFile t = new TestDataFile(new List<PeptideWithSetModifications> { pep1, pep2, pep3 });

            CompactPeptide peptide1 = new CompactPeptide(pep1, TerminusType.None);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan1 = t.GetOneBasedScan(2) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, peptide1.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(peptide1, 0, 3, 0, scan1);

            CompactPeptide peptide2 = new CompactPeptide(pep2, TerminusType.None);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan2 = t.GetOneBasedScan(4) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, peptide2.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(peptide2, 1, 2, 1, scan2);

            CompactPeptide peptide3 = new CompactPeptide(pep3, TerminusType.None);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan3 = t.GetOneBasedScan(6) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, peptide3.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm3 = new PeptideSpectralMatch(peptide3, 0, 1, 2, scan3);

            CompactPeptide peptide4 = new CompactPeptide(pep4, TerminusType.None);
            psm3.AddOrReplace(peptide4, 1, 1, true);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {
                    peptide1, new HashSet<PeptideWithSetModifications>{ pep1 }
                },
                {
                    peptide2, new HashSet<PeptideWithSetModifications>{ pep2 }
                },
                {
                    peptide3, new HashSet<PeptideWithSetModifications>{ pep3 }
                },
                {
                    peptide4, new HashSet<PeptideWithSetModifications>{ pep4 }
                },
            };
            psm1.MatchToProteinLinkedPeptides(matching);
            psm2.MatchToProteinLinkedPeptides(matching);
            psm3.MatchToProteinLinkedPeptides(matching);

            var newPsms = new List<PeptideSpectralMatch> { psm1, psm2, psm3 };
            CommonParameters cp = new CommonParameters
            {
                CalculateDeltaScore = false,
                CalculateEValue = true
            };
            FdrAnalysisEngine fdr = new FdrAnalysisEngine(newPsms, searchModes.NumNotches, cp, nestedIds);

            fdr.Run();

            Assert.AreEqual(2, searchModes.NumNotches);
            Assert.AreEqual(0, newPsms[0].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[0].FdrInfo.CumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[1].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[1].FdrInfo.CumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[2].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[2].FdrInfo.CumulativeTargetNotch);

            Assert.AreEqual(0, newPsms[0].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(1, newPsms[0].FdrInfo.CumulativeTarget);
            Assert.AreEqual(0, newPsms[1].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(2, newPsms[1].FdrInfo.CumulativeTarget);
            Assert.AreEqual(0, newPsms[2].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(3, newPsms[2].FdrInfo.CumulativeTarget);
        }

        #endregion Public Methods
    }
}