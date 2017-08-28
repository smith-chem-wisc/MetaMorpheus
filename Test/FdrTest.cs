using Chemistry;
using EngineLayer;
using EngineLayer.Analysis;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class FdrTest
    {
        #region Public Methods

        [Test]
        public static void FdrTestMethod()
        {
            List<MassDiffAcceptor> searchModes = new List<MassDiffAcceptor> { new DotMassDiffAcceptor(null, new List<double> { 0, 1.0029 }, new PpmTolerance(5)) };
            List<string> nestedIds = new List<string>();
            List<Psm>[] newPsms = new List<Psm>[1];

            Protein p = new Protein("MNKNNKNNNKNNNNK", null);
            DigestionParams digestionParams = new DigestionParams();
            var digested = p.Digest(digestionParams, new List<ModificationWithMass>()).ToList();

            PeptideWithSetModifications pep1 = digested[0].GetPeptidesWithSetModifications(digestionParams, new List<ModificationWithMass>()).First();
            PeptideWithSetModifications pep2 = digested[1].GetPeptidesWithSetModifications(digestionParams, new List<ModificationWithMass>()).First();
            PeptideWithSetModifications pep3 = digested[2].GetPeptidesWithSetModifications(digestionParams, new List<ModificationWithMass>()).First();
            PeptideWithSetModifications pep4 = digested[3].GetPeptidesWithSetModifications(digestionParams, new List<ModificationWithMass>()).First();

            TestDataFile t = new TestDataFile(new List<PeptideWithSetModifications> { pep1, pep2, pep3 });

            CompactPeptide peptide1 = new CompactPeptide(pep1, TerminusType.None);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan1 = t.GetOneBasedScan(2) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, new MzPeak(peptide1.MonoisotopicMassIncludingFixedMods.ToMz(1), 1), 1, null);
            Psm psm1 = new Psm(peptide1, 0, 3, 0, scan1);

            CompactPeptide peptide2 = new CompactPeptide(pep2, TerminusType.None);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan2 = t.GetOneBasedScan(4) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, new MzPeak(peptide2.MonoisotopicMassIncludingFixedMods.ToMz(1), 1), 1, null);
            Psm psm2 = new Psm(peptide2, 1, 2, 1, scan2);

            CompactPeptide peptide3 = new CompactPeptide(pep3, TerminusType.None);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan3 = t.GetOneBasedScan(6) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, new MzPeak(peptide3.MonoisotopicMassIncludingFixedMods.ToMz(1), 1), 1, null);
            Psm psm3 = new Psm(peptide3, 0, 1, 2, scan3);

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

            newPsms[0] = new List<Psm> { psm1, psm2, psm3 };
            FdrAnalysisEngine fdr = new FdrAnalysisEngine(newPsms, searchModes, nestedIds);

            fdr.Run();

            Assert.AreEqual(2, searchModes[0].NumNotches);
            Assert.AreEqual(0, newPsms[0][0].FdrInfo.cumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[0][0].FdrInfo.cumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[0][1].FdrInfo.cumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[0][1].FdrInfo.cumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[0][2].FdrInfo.cumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[0][2].FdrInfo.cumulativeTargetNotch);

            Assert.AreEqual(0, newPsms[0][0].FdrInfo.cumulativeDecoy);
            Assert.AreEqual(1, newPsms[0][0].FdrInfo.cumulativeTarget);
            Assert.AreEqual(0, newPsms[0][1].FdrInfo.cumulativeDecoy);
            Assert.AreEqual(2, newPsms[0][1].FdrInfo.cumulativeTarget);
            Assert.AreEqual(0, newPsms[0][2].FdrInfo.cumulativeDecoy);
            Assert.AreEqual(3, newPsms[0][2].FdrInfo.cumulativeTarget);
        }

        #endregion Public Methods
    }
}