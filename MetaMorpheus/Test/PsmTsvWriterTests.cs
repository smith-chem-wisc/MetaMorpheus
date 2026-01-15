using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using Omics.Digestion;
using Omics.Modifications;
using Easy.Common.Extensions;
using EngineLayer.SpectrumMatch;
using Readers;
using EngineLayer.FdrAnalysis;
using System.Linq;
using System.Reflection;

namespace Test
{
    internal static class PsmTsvWriterTests
    {
        [Test]
        public static void ResolveModificationsTest()
        {
            double mass = 12.0 + new PeptideWithSetModifications(new Protein("LNLDLDND", "prot1"),new DigestionParams(),1,8,CleavageSpecificity.Full,"",0, new Dictionary<int, Modification>(),0,null).MonoisotopicMass.ToMz(1);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), mass, 1, "", new CommonParameters());

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);

            Dictionary<DissociationType, List<double>> NeutralLosses = new Dictionary<DissociationType, List<double>>();
            NeutralLosses.Add(DissociationType.HCD, new List<double> { 0 });

            Modification modFormula_C1 = new Modification(_originalId: "modC", _accession: "", _modificationType: "mt", _featureType: "", _target: motif1, _locationRestriction: "Anywhere.", _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("C1")), null, null, null, null, _neutralLosses: NeutralLosses, null, null);
            Modification modFormula_H1 = new Modification(_originalId: "modH", _accession: "", _modificationType: "mt", _featureType: "", _target: motif1, _locationRestriction: "Anywhere.", _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("H1")), null, null, null, null, _neutralLosses: NeutralLosses, null, null);

            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification>{ modFormula_C1, modFormula_H1 }},
            };
            Protein protein1 = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: oneBasedModifications);

            Dictionary<int, Modification> allModsOneIsNterminus1 = new Dictionary<int, Modification>
            {
                {2, modFormula_C1},
            };

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus1, 0);

            Dictionary<int, Modification> allModsOneIsNterminus2 = new Dictionary<int, Modification>
            {
                {2,modFormula_H1 },
            };

            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus2, 0);

            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1), scoreCutoff: 1);
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", CommonParameters));

            List<MatchedFragmentIon> mfi = new List<MatchedFragmentIon>();

            //we're adding a neutral loss of 5 to the product to make sure we hit the right spot in the unit test to add that loss to the product ion string
            Product p = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 5);
            mfi.Add(new MatchedFragmentIon(p, 1, 1, 1));
            SpectralMatch myPsm = new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, new CommonParameters(), mfi);

            myPsm.AddOrReplace(pwsm2, 10, 0, true, mfi);

            myPsm.ResolveAllAmbiguities();

            //Here we have a situation where there are two mods at the same position with different chemical formuala. They cannot be resolved and so the return value is null.
            Assert.That(myPsm.ModsChemicalFormula, Is.Null);
            var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');

            string myPsmString = myPsm.ToString();
            string[] myPsmStringSplit = myPsmString.Split('\t');
            var ppmErrorIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.MassDiffPpm);
            string ppmErrorString = myPsmStringSplit[ppmErrorIndex];

            //The two different mods produce two separate mass errors, which are both then reported
            Assert.That(ppmErrorString, Is.EqualTo("0.00000|11801.30000"));

            //Make sure we see produt ion neutral losses in the output.
            var matchedIonSeriesIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.MatchedIonSeries);
            string matchedIonSeries = myPsmStringSplit[matchedIonSeriesIndex];
            Assert.That(matchedIonSeries, Is.EqualTo("[(b1-5.00)+1]"));

            //removing one of the peptides to reset for the next test
            var tentativeSpectralMatch = new SpectralMatchHypothesis(0, pwsm2, mfi, myPsm.Score);
            myPsm.RemoveThisAmbiguousPeptide(tentativeSpectralMatch);

            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus1, 0);
            myPsm.AddOrReplace(pwsm3, 10, 0, true, mfi);

            myPsm.ResolveAllAmbiguities();

            //Now we have removed one of the peptides with a different chemical formual and replaced it with a mod that has the same chemical formula as the remaining original best peptide
            //Here we have a situation where there are two mods at the same position have the same chemical formuala and they can be resolved and so the return value the chemical formual of the mod.
            Assert.That(myPsm.ModsChemicalFormula.Formula.ToString(), Is.EqualTo("C"));

            myPsmString = myPsm.ToString();
            myPsmStringSplit = myPsmString.Split('\t');
            ppmErrorString = myPsmStringSplit[ppmErrorIndex];

            Assert.That(ppmErrorString, Is.EqualTo("0"));
        }

        /// <summary>
        /// Test Case 1: Verifies that when peptide is null, all output fields are empty strings
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_PeptideNull()
        {
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, null, writePeptideLevelFdr: false);
            
            Assert.That(dict[SpectrumMatchFromTsvHeader.SpectralAngle], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.LocalizedScores], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.ImprovementPossible], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTarget], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoy], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValue], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTargetNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValueNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.PEP], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.PEP_QValue], Is.EqualTo(" "));
        }

        /// <summary>
        /// Test Case 2: Tests localized scores handling when null
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_LocalizedScoresNull()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.LocalizedScores = null;
            
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: false);
            
            Assert.That(dict[SpectrumMatchFromTsvHeader.LocalizedScores], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.ImprovementPossible], Is.EqualTo(" "));
        }

        /// <summary>
        /// Test Case 3: Tests localized scores calculation when present
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_LocalizedScoresPresent()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.LocalizedScores = new List<double> { 8.5, 9.0, 10.5 };
            
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: false);
            
            Assert.That(dict[SpectrumMatchFromTsvHeader.LocalizedScores], Is.Not.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.ImprovementPossible], Is.EqualTo("0.500"));
        }

        /// <summary>
        /// Test Case 4: Verifies FDR fields are empty when FdrInfo is null
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_FdrInfoNull()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.PsmFdrInfo = null;
            psm.PeptideFdrInfo = null;
            
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: false);
            
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTarget], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoy], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValue], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTargetNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValueNotch], Is.EqualTo(" "));
        }

        /// <summary>
        /// Test Case 5: Notch ambiguous with QValueNotch > 1, PSM level FDR, min == null
        /// Tests the case where BestMatchingBioPolymersWithSetMods returns null from MinBy
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_NotchAmbiguous_PsmLevel_MinNull()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.PsmFdrInfo = new FdrInfo
            {
                CumulativeTarget = 100,
                CumulativeDecoy = 5,
                QValue = 0.05,
                QValueNotch = 2.0,
                CumulativeTargetNotch = 0,
                CumulativeDecoyNotch = 0,
                PEP = 0.01,
                PEP_QValue = 0.02
            };
            psm.ResolveAllAmbiguities();

            // Use reflection to set _BestMatchingBioPolymersWithSetMods to an empty collection
            var fieldInfo = typeof(SpectralMatch).GetField("_BestMatchingBioPolymersWithSetMods",
                BindingFlags.NonPublic | BindingFlags.Instance);
            fieldInfo.SetValue(psm, new List<SpectralMatchHypothesis>());
            var notchFieldInfo = typeof(SpectralMatch).GetProperty("Notch",
                BindingFlags.Public | BindingFlags.Instance);
            notchFieldInfo.SetValue(psm, null);

            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: false);

            // When min == null (no BestMatchingBioPolymersWithSetMods), values stay as " "
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTargetNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValueNotch], Is.EqualTo(" "));
        }

        /// <summary>
        /// Test Case 6: Notch ambiguous with QValueNotch > 1, PSM level FDR, min != null but QValueNotch.HasValue == false
        /// Tests when hypothesis exists but QValueNotch is null
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_NotchAmbiguous_PsmLevel_QValueNotchNull()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            var peptide2 = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.AddOrReplace(peptide2, 10, 1, true, mfi);
            
            psm.PsmFdrInfo = new FdrInfo
            {
                CumulativeTarget = 100,
                CumulativeDecoy = 5,
                QValue = 0.05,
                QValueNotch = 2.0,
                PEP = 0.01,
                PEP_QValue = 0.02
            };
            
            // Don't set QValueNotch on hypotheses - they should be null
            psm.ResolveAllAmbiguities();
            
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: false);
            
            // When min.QValueNotch.HasValue == false, values stay as " "
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTargetNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch], Is.EqualTo(" "));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValueNotch], Is.EqualTo(" "));
        }

        /// <summary>
        /// Test Case 7: Notch ambiguous with QValueNotch > 1, PSM level FDR, successful resolution
        /// Tests when min != null and QValueNotch.HasValue == true
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_NotchAmbiguous_PsmLevel_SuccessfulResolution()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            var peptide2 = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.AddOrReplace(peptide2, 10, 1, true, mfi);
            
            psm.PsmFdrInfo = new FdrInfo
            {
                CumulativeTarget = 100,
                CumulativeDecoy = 5,
                QValue = 0.05,
                QValueNotch = 2.0,
                PEP = 0.01,
                PEP_QValue = 0.02
            };
            
            // Set QValueNotch on hypotheses
            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                hypothesis.QValueNotch = 0.03;
                hypothesis.CumulativeTargetNotch = 95;
                hypothesis.CumulativeDecoyNotch = 3;
            }
            
            psm.ResolveAllAmbiguities();
            
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: false);

            // Values should be populated from the hypothesis with minimum QValueNotch
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTargetNotch], Is.EqualTo("95.000000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch], Is.EqualTo("3.000000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValueNotch], Is.EqualTo("0.030000"));
        }

        /// <summary>
        /// Test Case 8: Notch ambiguous with QValueNotch > 1, PEPTIDE level FDR, successful resolution
        /// Tests when min != null and PeptideQValueNotch.HasValue == true
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_NotchAmbiguous_PeptideLevel_SuccessfulResolution()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            var peptide2 = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.AddOrReplace(peptide2, 10, 1, true, mfi);
            
            psm.PeptideFdrInfo = new FdrInfo
            {
                CumulativeTarget = 80,
                CumulativeDecoy = 4,
                QValue = 0.04,
                QValueNotch = 2.5,
                PEP = 0.015,
                PEP_QValue = 0.025
            };
            
            // Set PeptideQValueNotch on hypotheses
            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                hypothesis.PeptideQValueNotch = 0.025;
                hypothesis.PeptideCumulativeTargetNotch = 75;
                hypothesis.PeptideCumulativeDecoyNotch = 2;
            }
            
            psm.ResolveAllAmbiguities();
            
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: true);
            
            // Values should be populated from the hypothesis with minimum PeptideQValueNotch
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTargetNotch], Is.EqualTo("75.000000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch], Is.EqualTo("2.000000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValueNotch], Is.EqualTo("0.025000"));
        }

        /// <summary>
        /// Test Case 9: Tests the else branch when Notch is NOT ambiguous
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_NotchNotAmbiguous()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.PsmFdrInfo = new FdrInfo
            {
                CumulativeTarget = 100,
                CumulativeDecoy = 5,
                QValue = 0.05,
                QValueNotch = 0.06,
                CumulativeTargetNotch = 98,
                CumulativeDecoyNotch = 4,
                PEP = 0.01,
                PEP_QValue = 0.02
            };
            psm.ResolveAllAmbiguities();
            
            var dict = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchScoreData(dict, psm, writePeptideLevelFdr: false);
            
            // Values should come directly from FdrInfo
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTarget], Is.EqualTo("100"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoy], Is.EqualTo("5"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValue], Is.EqualTo("0.050000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeTargetNotch], Is.EqualTo("98.000000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch], Is.EqualTo("4.000000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.QValueNotch], Is.EqualTo("0.060000"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.PEP], Is.EqualTo("0.01"));
            Assert.That(dict[SpectrumMatchFromTsvHeader.PEP_QValue], Is.EqualTo("0.02"));
        }

        /// <summary>
        /// Test Case 10: Complete end-to-end integration test with ToString method
        /// Tests notch ambiguity resolution through the full ToString pipeline
        /// </summary>
        [Test]
        public static void TestAddMatchScoreData_ToStringIntegration_NotchAmbiguity()
        {
            var protein = new Protein("PEPTIDESEQUENCE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            var peptide2 = new PeptideWithSetModifications(protein, digestionParams, 1, 7, 
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            
            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);
            psm.AddOrReplace(peptide2, 10, 1, true, mfi);
            psm.LocalizedScores = new List<double> { 8.0, 9.5 };
            
            psm.PsmFdrInfo = new FdrInfo
            {
                CumulativeTarget = 100,
                CumulativeDecoy = 5,
                QValue = 0.05,
                QValueNotch = 2.0,
                PEP = 0.01,
                PEP_QValue = 0.02
            };
            
            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                hypothesis.QValueNotch = 0.04;
                hypothesis.CumulativeTargetNotch = 90;
                hypothesis.CumulativeDecoyNotch = 3;
            }
            
            psm.ResolveAllAmbiguities();
            
            var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
            string psmString = psm.ToString(new Dictionary<string, int>(), writePeptideLevelFdr: false);
            string[] psmStringSplit = psmString.Split('\t');
            
            var cumulativeTargetNotchIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.CumulativeTargetNotch);
            var qValueNotchIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.QValueNotch);
            var localizedScoresIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.LocalizedScores);
            
            Assert.That(psmStringSplit[cumulativeTargetNotchIndex], Is.EqualTo("90.000000"));
            Assert.That(psmStringSplit[qValueNotchIndex], Is.EqualTo("0.040000"));
            Assert.That(psmStringSplit[localizedScoresIndex], Contains.Substring("8.000"));
            Assert.That(psmStringSplit[localizedScoresIndex], Contains.Substring("9.500"));
        }
    }
}