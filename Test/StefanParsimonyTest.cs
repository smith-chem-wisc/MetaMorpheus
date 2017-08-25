using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class StefanParsimonyTest
    {
        #region Public Methods

        [Test]
        public static void ParsimonyVariableTreatAsUnique()
        {
            bool localizeable = false;
            var hah = GetInfo(localizeable);

            var newPsms = hah.Item1;
            var compactPeptideToProteinPeptideMatching = hah.Item2;
            var massDiffAcceptors = hah.Item3;
            var noOneHitWonders = hah.Item4;
            var compactPeptide1 = hah.Item5;
            var compactPeptide2 = hah.Item6;

            Assert.AreEqual(3, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            bool modPeptidesAreUnique = true;
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique, new List<string>());
            pae.Run();

            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            Assert.AreEqual(2, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide1].Select(b => b.Protein)).Count);
            Assert.AreEqual(2, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide2].Select(b => b.Protein)).Count);
        }

        [Test]
        public static void ParsimonyVariableDontTreatAsUnique()
        {
            bool localizeable = false;
            var hah = GetInfo(localizeable);

            var newPsms = hah.Item1;
            var compactPeptideToProteinPeptideMatching = hah.Item2;
            var massDiffAcceptors = hah.Item3;
            var noOneHitWonders = hah.Item4;
            var compactPeptide1 = hah.Item5;
            var compactPeptide2 = hah.Item6;

            Assert.AreEqual(3, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            bool modPeptidesAreUnique = false;
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique, new List<string>());
            pae.Run();

            Assert.AreEqual(4, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(4, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            Assert.AreEqual(2, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide1].Select(b => b.Protein)).Count);
            Assert.AreEqual(2, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide2].Select(b => b.Protein)).Count);
        }

        [Test]
        public static void ParsimonyLocalizeableTreatAsUnique()
        {
            bool localizeable = true;
            var hah = GetInfo(localizeable);

            var newPsms = hah.Item1;
            var compactPeptideToProteinPeptideMatching = hah.Item2;
            var massDiffAcceptors = hah.Item3;
            var noOneHitWonders = hah.Item4;
            var compactPeptide1 = hah.Item5;
            var compactPeptide2 = hah.Item6;

            Assert.AreEqual(3, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            bool modPeptidesAreUnique = true;
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique, new List<string>());
            pae.Run();

            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            Assert.AreEqual(1, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide1].Select(b => b.Protein)).Count);
            Assert.AreEqual(1, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide2].Select(b => b.Protein)).Count);
        }

        [Test]
        public static void ParsimonyLocalizeableDontTreatAsUnique()
        {
            bool localizeable = true;
            var hah = GetInfo(localizeable);

            var newPsms = hah.Item1;
            var compactPeptideToProteinPeptideMatching = hah.Item2;
            var massDiffAcceptors = hah.Item3;
            var noOneHitWonders = hah.Item4;
            var compactPeptide1 = hah.Item5;
            var compactPeptide2 = hah.Item6;

            Assert.AreEqual(3, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            bool modPeptidesAreUnique = false;
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique, new List<string>());
            pae.Run();

            Assert.AreEqual(4, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(4, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            Assert.AreEqual(2, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide1].Select(b => b.Protein)).Count);
            Assert.AreEqual(2, new HashSet<Protein>(compactPeptideToProteinPeptideMatching[compactPeptide2].Select(b => b.Protein)).Count);
        }

        [Test]
        public static void ParsimonyWeirdCatch()
        {
            Protein protein1 = new Protein("MATSIK", "protein1", isDecoy: true);
            Protein protein2 = new Protein("MATSLK", "protein2");
            Protein protein3 = new Protein("MTASIK", "protein3");

            IEnumerable<ModificationWithMass> allKnownFixedModifications = new List<ModificationWithMass>();
            DigestionParams digestionParams = new DigestionParams();
            var from1 = protein1.Digest(digestionParams, allKnownFixedModifications).First();
            var from2 = protein2.Digest(digestionParams, allKnownFixedModifications).First();
            var from3 = protein3.Digest(digestionParams, allKnownFixedModifications).First();

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            PeptideWithSetModifications pep1 = from1.GetPeptidesWithSetModifications(digestionParams,variableModifications).First();
            PeptideWithSetModifications pep2 = from2.GetPeptidesWithSetModifications(digestionParams,variableModifications).First();
            PeptideWithSetModifications pep3 = from3.GetPeptidesWithSetModifications(digestionParams,variableModifications).First();

            CompactPeptide compactPeptide1 = pep1.CompactPeptide(TerminusType.None);
            CompactPeptide compactPeptide2 = pep2.CompactPeptide(TerminusType.None);
            CompactPeptide compactPeptide3 = pep3.CompactPeptide(TerminusType.None);

            Assert.AreEqual(compactPeptide1, compactPeptide2);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {compactPeptide1, new HashSet<PeptideWithSetModifications>{pep1, pep2} },
                {compactPeptide3, new HashSet<PeptideWithSetModifications>{pep3} }
            };

            var cool = (ProteinParsimonyResults)new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, false, new List<string>()).Run();

            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching.Count);

            // Only 1 because the target is removed!!!
            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide3].Count);
        }

        #endregion Public Methods

        #region Private Methods

        private static Tuple<List<Psm>[], Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>, List<MassDiffAcceptor>, bool, CompactPeptideBase, CompactPeptideBase> GetInfo(bool localizeable)
        {
            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams
                {
                    MinPeptideLength = null,
                    MaxMissedCleavages = 0,
                    InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                    MaxModsForPeptide = 1,
                    MaxModificationIsoforms = 2,
                },
                ScoreCutoff = 1
            };

            // Alanine = Glycine + CH2
            Protein protein1 = new Protein("MA", "protein1");
            Protein protein2 = new Protein("MG", "protein2");
            Protein protein3;
            double monoisotopicMass = Chemistry.ChemicalFormula.ParseFormula("CH2").MonoisotopicMass;
            ModificationMotif.TryGetMotif("G", out ModificationMotif motif1);
            ModificationMotif.TryGetMotif("A", out ModificationMotif motif2);
            TerminusLocalization modificationSites = TerminusLocalization.Any;
            List<ModificationWithMass> allKnownFixedModifications = new List<ModificationWithMass>
            {
                new ModificationWithMass("CH2 on Glycine", null, motif1, modificationSites, monoisotopicMass)
            };
            List<ModificationWithMass> variableModifications;

            ModificationWithMass alanineMod = new ModificationWithMass("CH2 on Alanine", null, motif2, modificationSites, monoisotopicMass);

            if (localizeable)
            {
                variableModifications = new List<ModificationWithMass>();
                IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{alanineMod} }
                };
                protein3 = new Protein("MA", "protein3", oneBasedModifications: oneBasedModifications);
            }
            else
            {
                variableModifications = new List<ModificationWithMass>();
                variableModifications = new List<ModificationWithMass> { alanineMod };
                protein3 = new Protein("MA", "protein3");
            }

            var prot1List = protein1.Digest(CommonParameters.DigestionParams, allKnownFixedModifications);
            PeptideWithPossibleModifications pepWithPossibleModifications1 = prot1List.First();
            var pep1list = pepWithPossibleModifications1.GetPeptidesWithSetModifications(CommonParameters.DigestionParams, variableModifications);
            PeptideWithSetModifications pepWithSetModifications1 = pep1list.First();

            var prot2List = protein2.Digest(CommonParameters.DigestionParams, allKnownFixedModifications);
            PeptideWithPossibleModifications pepWithPossibleModifications2 = prot2List.First();
            var pep2list = pepWithPossibleModifications2.GetPeptidesWithSetModifications(CommonParameters.DigestionParams, variableModifications);
            PeptideWithSetModifications pepWithSetModifications2 = pep2list.First();

            var prot3List = protein3.Digest(CommonParameters.DigestionParams, allKnownFixedModifications);
            PeptideWithPossibleModifications pepWithPossibleModifications3 = prot3List.First();
            var pep3list = pepWithPossibleModifications3.GetPeptidesWithSetModifications(CommonParameters.DigestionParams, variableModifications).ToList();
            PeptideWithSetModifications pepWithSetModifications3 = pep3list.Last();

            CompactPeptide compactPeptide1 = new CompactPeptide(pepWithSetModifications1, TerminusType.None);
            CompactPeptide compactPeptideDuplicate = new CompactPeptide(pepWithSetModifications2, TerminusType.None);
            Assert.AreEqual(compactPeptide1, compactPeptideDuplicate);
            CompactPeptide compactPeptide2 = new CompactPeptide(pepWithSetModifications3, TerminusType.None);

            List<Psm>[] newPsms = new List<Psm>[1];
            string fullFilePath = null;
            double intensity = 0;
            double mz = 0;
            IMzPeak precursorMonoisotopicPeak = new MzPeak(mz, intensity);
            int precursorCharge = 0;
            TestDataFile testDataFile = new TestDataFile();
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan = testDataFile.GetOneBasedScan(2) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(mzLibScan, precursorMonoisotopicPeak, precursorCharge, fullFilePath);
            int scanIndex = 0;
            double score = 0;
            int notch = 0;
            Psm psm1 = new Psm(compactPeptide1, notch, score, scanIndex, scan);
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0);
            Psm psm2 = new Psm(compactPeptide1, notch, score, scanIndex, scan);
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0);
            Psm psm3 = new Psm(compactPeptide2, notch, score, scanIndex, scan);
            psm3.SetFdrValues(0, 0, 0, 0, 0, 0);
            newPsms[0] = new List<Psm>
            {
                psm1,
                psm2,
                psm3
            };

            List<MassDiffAcceptor> massDiffAcceptors = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            SequencesToActualProteinPeptidesEngine stappe = new SequencesToActualProteinPeptidesEngine(newPsms, new List<Protein> { protein1, protein2, protein3 }, allKnownFixedModifications, variableModifications, TerminusType.None, new List<DigestionParams> { CommonParameters.DigestionParams }, new List<string>());

            var haha = (SequencesToActualProteinPeptidesEngineResults)stappe.Run();
            var compactPeptideToProteinPeptideMatching = haha.CompactPeptideToProteinPeptideMatching;

            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching.Count);

            psm1.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            bool noOneHitWonders = false;

            return new Tuple<List<Psm>[], Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>, List<MassDiffAcceptor>, bool, CompactPeptideBase, CompactPeptideBase>
            (
                newPsms, compactPeptideToProteinPeptideMatching, massDiffAcceptors, noOneHitWonders, compactPeptide1, compactPeptide2
            );
        }

        #endregion Private Methods
    }
}