using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using Omics.Digestion;
using Omics.Modifications;
using Easy.Common.Extensions;
using EngineLayer.SpectrumMatch;
using Readers;
using Readers.ProForma;
using EngineLayer.FdrAnalysis;
using GuiFunctions;
using System.Linq;
using System.Reflection;
using Transcriptomics;
using Transcriptomics.Digestion;

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

        #region Resolve Null Handling

        /// <summary>
        /// Resolve(IEnumerable&lt;string&gt;) with no usable value at all (all null) resolves to a null
        /// ResolvedString - the caller then emits an empty cell rather than crashing.
        /// </summary>
        [Test]
        public static void Resolve_String_AllNullValues_ReturnsNullResolvedString()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string> { null, null, null });
            Assert.That(resolvedString, Is.Null);
            Assert.That(resolvedValue, Is.Null);
        }

        [Test]
        public static void Resolve_String_EmptyCollection_ReturnsNullResolvedString()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string>());
            Assert.That(resolvedString, Is.Null);
            Assert.That(resolvedValue, Is.Null);
        }

        [Test]
        public static void Resolve_String_SingleNullElement_ReturnsNullResolvedString()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string> { null });
            Assert.That(resolvedString, Is.Null);
            Assert.That(resolvedValue, Is.Null);
        }

        /// <summary>
        /// A null element among equal non-null values cannot be collapsed (first.Equals(null) is
        /// false), so Resolve joins the raw list, rendering the null element as an empty cell.
        /// </summary>
        [Test]
        public static void Resolve_String_NullAmongSameValues_JoinsWithEmptyCell()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string> { "a", null, "a" });
            Assert.That(resolvedString, Is.EqualTo("a||a"));
            Assert.That(resolvedValue, Is.Null);
        }

        [Test]
        public static void Resolve_String_NullAmongDistinctValues_JoinsWithEmptyCell()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string> { "a", null, "b" });
            Assert.That(resolvedString, Is.EqualTo("a||b"));
            Assert.That(resolvedValue, Is.Null);
        }

        /// <summary>
        /// The ambiguousIfNull overload: an all-null list still resolves to a null ResolvedString even
        /// when a fallback is supplied - the fallback only selects the joining strategy, it does not
        /// rescue an all-null input.
        /// </summary>
        [Test]
        public static void Resolve_StringWithFallback_AllNullValues_ReturnsNullResolvedString()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string> { null, null }, "fallback");
            Assert.That(resolvedString, Is.Null);
            Assert.That(resolvedValue, Is.Null);
        }

        /// <summary>
        /// The ambiguousIfNull overload joins only distinct values when a fallback is present, with the
        /// null element rendered as an empty cell.
        /// </summary>
        [Test]
        public static void Resolve_StringWithFallback_NullAmongValues_JoinsDistinctWithEmptyCell()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string> { "b", null, "a", "b" }, "fallback");
            Assert.That(resolvedString, Is.EqualTo("b||a"));
            Assert.That(resolvedValue, Is.Null);
        }

        /// <summary>
        /// The ambiguousIfNull overload keeps duplicated values (and the empty cell from the null
        /// element) when no fallback is supplied.
        /// </summary>
        [Test]
        public static void Resolve_StringWithoutFallback_NullAmongSameValues_JoinsDuplicatesWithEmptyCell()
        {
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(new List<string> { "a", null, "a" }, ambiguousIfNull: null);
            Assert.That(resolvedString, Is.EqualTo("a||a"));
            Assert.That(resolvedValue, Is.Null);
        }

        /// <summary>
        /// ModsChemicalFormulas / ModsCombinedChemicalFormula null-path: any hypothesis whose
        /// modification list contains a null Modification resolves to "unknown" rather than throwing.
        /// </summary>
        [Test]
        public static void Resolve_Modifications_NullModificationElement_ReturnsUnknown()
        {
            var enumerable = new List<List<Modification>>
            {
                new List<Modification> { null }
            };
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(enumerable);
            Assert.That(resolvedString, Is.EqualTo("unknown"));
            Assert.That(resolvedValue, Is.Null);
        }

        /// <summary>
        /// A modification whose ChemicalFormula is null hits the same guard and also resolves to
        /// "unknown" rather than throwing.
        /// </summary>
        [Test]
        public static void Resolve_Modifications_NullChemicalFormula_ReturnsUnknown()
        {
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif);
            Modification modWithoutFormula = new Modification(_originalId: "noFormula", _modificationType: "mt",
                _target: motif, _locationRestriction: "Anywhere.");
            Assert.That(modWithoutFormula.ChemicalFormula, Is.Null);

            var enumerable = new List<List<Modification>>
            {
                new List<Modification> { modWithoutFormula }
            };
            var (resolvedString, resolvedValue) = PsmTsvWriter.Resolve(enumerable);
            Assert.That(resolvedString, Is.EqualTo("unknown"));
            Assert.That(resolvedValue, Is.Null);
        }

        #endregion

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

        /// <summary>
        /// Top-down (AnalyteType.Proteoform): the ProForma column is emitted immediately after Full
        /// Sequence in both header and data row, and carries the proteoform's ProForma 2.0 string.
        /// Both AllPSMs and AllProteoforms go through this same GetTabSeparatedHeader/ToString path,
        /// so this covers the column in both output files.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void ProForma_TopDownSearch_ColumnFollowsFullSequenceAndCarriesProFormaString()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Proteoform;
            try
            {
                ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
                Modification oxidation = new Modification(_originalId: "Oxidation", _modificationType: "Common Variable",
                    _target: motif, _locationRestriction: "Anywhere.",
                    _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("O1")));

                var allModsOneIsNterminus = new Dictionary<int, Modification> { { 2, oxidation } }; // residue 1 (M)
                Protein protein = new Protein("MPEPTIDEK", "prot_td");
                PeptideWithSetModifications proteoform = new PeptideWithSetModifications(
                    protein, new DigestionParams(), 1, 9, CleavageSpecificity.Full, "", 0, allModsOneIsNterminus, 0);

                double mass = 12.0 + proteoform.MonoisotopicMass.ToMz(1);
                var scan = new Ms2ScanWithSpecificMass(
                    new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                        0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                    mass, 1, "", new CommonParameters());

                var psm = new PeptideSpectralMatch(proteoform, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();

                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                int fullSeqIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.FullSequence);
                int proFormaIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ProForma);

                // Column is present and sits immediately after Full Sequence.
                Assert.That(proFormaIndex, Is.EqualTo(fullSeqIndex + 1));

                string[] rowSplits = psm.ToString(new Dictionary<string, int>()).Split('\t');

                // Header and data row stay in sync - the invariant the single analyte-type gate exists to protect.
                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));

                // Compared against a literal rather than ToProFormaString(), so the assertion cannot move
                // with the production code it is checking.
                Assert.That(rowSplits[proFormaIndex], Is.EqualTo("M[Oxidation]PEPTIDEK"));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// Non-top-down runs: the ProForma column must appear in neither the header nor the data row.
        /// Asserting the row's field count matches the header's is what actually defends the invariant -
        /// a gate that diverged between the two would shift every downstream column by one.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        [TestCase(AnalyteType.Peptide)]
        [TestCase(AnalyteType.Oligo)]
        public static void ProForma_NonTopDownSearch_ColumnAbsentFromHeader(AnalyteType analyteType)
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = analyteType;
            try
            {
                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                Assert.That(headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ProForma), Is.EqualTo(-1));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// Bottom-up: the column is absent from the header AND the data row, and the two stay aligned.
        /// Only AnalyteType.Peptide is exercised with a real row - a PeptideSpectralMatch written under
        /// AnalyteType.Oligo is not a state any run produces (oligo runs match OSMs), and the pre-existing
        /// sequence-variation gate writes columns for it that the oligo header omits.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void ProForma_BottomUpSearch_ColumnAbsentFromRowAndRowStaysAligned()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            try
            {
                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                var psm = BuildProFormaTestPsm(out PeptideWithSetModifications peptide);
                string[] rowSplits = psm.ToString(new Dictionary<string, int>()).Split('\t');

                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));
                Assert.That(rowSplits, Does.Not.Contain(peptide.ToProFormaString()));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// The PR claims the ProForma value is resolved across ambiguous matches with the same Resolve(...)
        /// idiom as the neighbouring sequence columns. With two proteoforms whose ProForma strings differ,
        /// that idiom joins them with '|' - the branch a single-proteoform test never enters.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void ProForma_AmbiguousProteoforms_ValuesJoinedByPipe()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Proteoform;
            try
            {
                ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
                Modification oxidation = new Modification(_originalId: "Oxidation", _modificationType: "Common Variable",
                    _target: motif, _locationRestriction: "Anywhere.",
                    _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("O1")));

                Protein protein = new Protein("MPEPTIDEKM", "prot_td_ambig");
                // Same base sequence, oxidation on a different M - so the two ProForma strings differ.
                var modsOnFirstM = new Dictionary<int, Modification> { { 2, oxidation } };
                var modsOnLastM = new Dictionary<int, Modification> { { 11, oxidation } };

                PeptideWithSetModifications proteoformA = new PeptideWithSetModifications(
                    protein, new DigestionParams(), 1, 10, CleavageSpecificity.Full, "", 0, modsOnFirstM, 0);
                PeptideWithSetModifications proteoformB = new PeptideWithSetModifications(
                    protein, new DigestionParams(), 1, 10, CleavageSpecificity.Full, "", 0, modsOnLastM, 0);

                double mass = 12.0 + proteoformA.MonoisotopicMass.ToMz(1);
                var scan = new Ms2ScanWithSpecificMass(
                    new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                        0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                    mass, 1, "", new CommonParameters());

                var psm = new PeptideSpectralMatch(proteoformA, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
                psm.AddOrReplace(proteoformB, 10, 0, true, new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();

                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                int proFormaIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ProForma);
                string[] rowSplits = psm.ToString(new Dictionary<string, int>()).Split('\t');

                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));
                Assert.That(rowSplits[proFormaIndex], Does.Contain("|"));
                Assert.That(rowSplits[proFormaIndex], Does.Contain(proteoformA.ToProFormaString()));
                Assert.That(rowSplits[proFormaIndex], Does.Contain(proteoformB.ToProFormaString()));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// Builds the single-peptide PSM shared by the ProForma column tests.
        /// </summary>
        private static PeptideSpectralMatch BuildProFormaTestPsm(out PeptideWithSetModifications peptide)
        {
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            Modification oxidation = new Modification(_originalId: "Oxidation", _modificationType: "Common Variable",
                _target: motif, _locationRestriction: "Anywhere.",
                _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("O1")));

            var allModsOneIsNterminus = new Dictionary<int, Modification> { { 2, oxidation } }; // residue 1 (M)
            Protein protein = new Protein("MPEPTIDEK", "prot_td");
            peptide = new PeptideWithSetModifications(
                protein, new DigestionParams(), 1, 9, CleavageSpecificity.Full, "", 0, allModsOneIsNterminus, 0);

            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            return psm;
        }

        /// <summary>
        /// The MetaDraw grid keys its ProForma column off the loaded results rather than
        /// GlobalVariables.AnalyteType. These cover the predicate that makes that decision.
        /// </summary>
        [Test]
        public static void ShouldShowProFormaColumn_NullOrEmptyCollection_ReturnsFalse()
        {
            Assert.That(MetaDrawLogic.ShouldShowProFormaColumn(null), Is.False);
            Assert.That(MetaDrawLogic.ShouldShowProFormaColumn(new List<SpectrumMatchFromTsv>()), Is.False);
        }

        #region Oligo Terminus Columns

        /// <summary>
        /// Builds a single-hypothesis OSM for the oligo terminus column tests.
        /// The oligo carries explicit, non-default 5'- and 3'-terminus formulas so the test
        /// asserts the exact strings written, not the reader's fallbacks.
        /// </summary>
        private static OligoSpectralMatch BuildTerminiTestOsm(out OligoWithSetMods oligo)
        {
            ChemicalFormula fivePrime = ChemicalFormula.ParseFormula("O-3P-1");
            ChemicalFormula threePrime = ChemicalFormula.ParseFormula("H2O4P");
            oligo = BuildTestOligo("AAA", fivePrime, threePrime);
            return BuildOsmFromOligo(oligo);
        }

        private static OligoSpectralMatch BuildOsmFromOligo(OligoWithSetMods oligo)
        {
            double mass = 12.0 + oligo.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                mass, 1, "", new CommonParameters());

            var osm = new OligoSpectralMatch(oligo, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            osm.ResolveAllAmbiguities();
            return osm;
        }

        /// <summary>
        /// Builds an oligo over a parent with a non-null GeneNames list, so the full
        /// AddPeptideSequenceData path (the geneString column) does not throw during serialization.
        /// </summary>
        private static OligoWithSetMods BuildTestOligo(string baseSequence, IHasChemicalFormula fivePrime, IHasChemicalFormula threePrime)
        {
            var rna = new RNA(baseSequence, "test_rna", geneNames: new List<System.Tuple<string, string>>());
            return new OligoWithSetMods(baseSequence, n: rna, fivePrimeTerminus: fivePrime, threePrimeTerminus: threePrime);
        }

        /// <summary>
        /// Oligo runs: the 5'-Terminus column is emitted immediately after Essential Sequence,
        /// followed by 3'-Terminus, in both header and data row, and both carry the terminus
        /// chemical-formula strings. Both AllOSMs and AllOligos go through this same
        /// GetTabSeparatedHeader/ToString path, so this covers the columns in both output files.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void Termini_OligoSearch_ColumnsFollowEssentialSequenceAndCarryTerminusFormulas()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            try
            {
                var osm = BuildTerminiTestOsm(out OligoWithSetMods oligo);

                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                int essentialSequenceIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.EssentialSequence);
                int fivePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.FivePrimeTerminus);
                int threePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ThreePrimeTerminus);

                // 5' first, then 3', immediately after Essential Sequence.
                Assert.That(fivePrimeIndex, Is.EqualTo(essentialSequenceIndex + 1));
                Assert.That(threePrimeIndex, Is.EqualTo(fivePrimeIndex + 1));

                string[] rowSplits = osm.ToString(new Dictionary<string, int>()).Split('\t');

                // Header and data row stay in sync - the invariant the single analyte-type gate exists to protect.
                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));

                // Compared against the formula strings directly, so the assertions cannot move with
                // the production code they are checking.
                Assert.That(rowSplits[fivePrimeIndex], Is.EqualTo(oligo.FivePrimeTerminus.ThisChemicalFormula.Formula));
                Assert.That(rowSplits[threePrimeIndex], Is.EqualTo(oligo.ThreePrimeTerminus.ThisChemicalFormula.Formula));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// Non-oligo runs: the terminus columns must appear in neither the header nor the data row.
        /// Asserting the row's field count matches the header's is what actually defends the invariant -
        /// a gate that diverged between the two would shift every downstream column by one.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        [TestCase(AnalyteType.Peptide)]
        [TestCase(AnalyteType.Proteoform)]
        public static void Termini_NonOligoSearch_ColumnsAbsent(AnalyteType analyteType)
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = analyteType;
            try
            {
                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');

                Assert.That(headerSplits.IndexOf(SpectrumMatchFromTsvHeader.FivePrimeTerminus), Is.EqualTo(-1));
                Assert.That(headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ThreePrimeTerminus), Is.EqualTo(-1));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// Bottom-up (AnalyteType.Peptide): the terminus columns are absent from the header AND the
        /// data row, and the two stay aligned - so a PSM file never carries a trailing empty pair of
        /// cells and never shifts any existing column.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void Termini_BottomUpSearch_ColumnsAbsentFromRowAndRowStaysAligned()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            try
            {
                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                var psm = BuildProFormaTestPsm(out PeptideWithSetModifications peptide);
                string[] rowSplits = psm.ToString(new Dictionary<string, int>()).Split('\t');

                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));
                Assert.That(rowSplits, Does.Not.Contain(SpectrumMatchFromTsvHeader.FivePrimeTerminus));
                Assert.That(rowSplits, Does.Not.Contain(SpectrumMatchFromTsvHeader.ThreePrimeTerminus));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// The terminus columns resolve across ambiguous matches with the same Resolve(...) idiom as
        /// the neighbouring sequence columns. With two oligos whose terminus formulas differ, that
        /// idiom joins them with '|' - the branch a single-oligo test never enters.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void Termini_AmbiguousOligos_ValuesJoinedByPipe()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            try
            {
                // Same base sequence, different 5'-terminus - so the two terminus strings differ.
                ChemicalFormula minusPhosphateFivePrime = ChemicalFormula.ParseFormula("O-3P-1");
                ChemicalFormula otherFivePrime = ChemicalFormula.ParseFormula("PO4");
                ChemicalFormula otherThreePrime = ChemicalFormula.ParseFormula("HO");
                OligoWithSetMods oligoA = BuildTestOligo("AAA", minusPhosphateFivePrime, otherThreePrime);
                OligoWithSetMods oligoB = BuildTestOligo("AAA", otherFivePrime, otherThreePrime);

                double mass = 12.0 + oligoA.MonoisotopicMass.ToMz(1);
                var scan = new Ms2ScanWithSpecificMass(
                    new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                        0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                    mass, 1, "", new CommonParameters());

                var ambiguousOsm = new OligoSpectralMatch(oligoA, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
                ambiguousOsm.AddOrReplace(oligoB, 10, 0, true, new List<MatchedFragmentIon>());
                ambiguousOsm.ResolveAllAmbiguities();

                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                int fivePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.FivePrimeTerminus);
                int threePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ThreePrimeTerminus);
                string[] rowSplits = ambiguousOsm.ToString(new Dictionary<string, int>()).Split('\t');

                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));
                Assert.That(rowSplits[fivePrimeIndex], Does.Contain("|"));
                Assert.That(rowSplits[fivePrimeIndex], Does.Contain(minusPhosphateFivePrime.Formula));
                Assert.That(rowSplits[fivePrimeIndex], Does.Contain(otherFivePrime.Formula));
                // The 3'-terminus is identical across both hypotheses, so it resolves to a single value.
                Assert.That(rowSplits[threePrimeIndex], Is.EqualTo(otherThreePrime.Formula));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void Termini_Oligo_RoundTripWritesAndReadsTerminusFormulas()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "TerminiRoundTrip_" + System.Guid.NewGuid().ToString("N") + ".osmtsv");
            try
            {
                // The same base sequence, with explicitly different terminus formulas per row. The
                // formulas are deliberately not the NucleicAcid/Rnase defaults so that, if the reader
                // silently used its PreviousResidue/NextResidue fallback instead of the written
                // columns, the assertions below would fail.
                ChemicalFormula rowAFive = ChemicalFormula.ParseFormula("HO");
                ChemicalFormula rowAThree = ChemicalFormula.ParseFormula("PO4");
                ChemicalFormula rowBFive = ChemicalFormula.ParseFormula("PO3");
                ChemicalFormula rowBThree = ChemicalFormula.ParseFormula("H2O");

                var osmA = BuildOsmFromOligo(BuildTestOligo("AAA", rowAFive, rowAThree));
                var osmB = BuildOsmFromOligo(BuildTestOligo("AAA", rowBFive, rowBThree));

                // Real OSM output always carries FDR values; the mzLib reader requires the QValue cell
                // (written from PsmFdrInfo, null for a fresh synthetic match) to be a parseable number.
                foreach (var osm in new[] { osmA, osmB })
                {
                    osm.PsmFdrInfo = new FdrInfo();
                }

                using (StreamWriter writer = new StreamWriter(filePath))
                {
                    writer.WriteLine(SpectralMatch.GetTabSeparatedHeader());
                    writer.WriteLine(osmA.ToString(new Dictionary<string, int>()));
                    writer.WriteLine(osmB.ToString(new Dictionary<string, int>()));
                }

                List<OsmFromTsv> readBack = SpectrumMatchTsvReader.ReadOsmTsv(filePath, out List<string> warnings);

                Assert.That(warnings, Is.Empty, "reading the round-tripped file produced warnings");
                Assert.That(readBack, Has.Count.EqualTo(2));

                // Each row must come back with exactly the terminus formulas its own writer row carried.
                Assert.That(readBack[0].FivePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(rowAFive.Formula));
                Assert.That(readBack[0].ThreePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(rowAThree.Formula));
                Assert.That(readBack[1].FivePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(rowBFive.Formula));
                Assert.That(readBack[1].ThreePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(rowBThree.Formula));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
                if (File.Exists(filePath))
                {
                    File.Delete(filePath);
                }
            }
        }

        /// <summary>
        /// The terminus columns are gated on AnalyteType.Oligo but the values are pulled from the
        /// matched hypothesis via a defensive cast to OligoWithSetMods. An OSM whose hypothesis is not
        /// an oligo exercises that null path end-to-end: the terminus cells resolve to null and are
        /// written as empty cells, with the rest of the row untouched.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void Termini_Oligo_NonOligoHypothesis_WritesEmptyCells()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            try
            {
                Protein protein = new Protein("MPEPTIDEK", "prot_td");
                PeptideWithSetModifications peptide = new PeptideWithSetModifications(
                    protein, new DigestionParams(), 1, 9, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

                double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
                var scan = new Ms2ScanWithSpecificMass(
                    new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                        0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""),
                    mass, 1, "", new CommonParameters());

                var osm = new OligoSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
                osm.ResolveAllAmbiguities();

                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                int fivePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.FivePrimeTerminus);
                int threePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ThreePrimeTerminus);
                string[] rowSplits = osm.ToString(new Dictionary<string, int>()).Split('\t');

                Assert.That(fivePrimeIndex, Is.GreaterThanOrEqualTo(0));
                Assert.That(threePrimeIndex, Is.GreaterThanOrEqualTo(0));
                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));
                // Defensive cast hits null -> Resolve(all-null) -> empty cells, not a crash.
                Assert.That(rowSplits[fivePrimeIndex], Is.EqualTo(string.Empty));
                Assert.That(rowSplits[threePrimeIndex], Is.EqualTo(string.Empty));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// The terminal-cap columns participate in the same null contract as every other sequence
        /// column: when there is no spectral match at all (pepWithModsIsNull, the header row),
        /// AddPeptideSequenceData emits a single blank cell for each cap column instead of invoking
        /// Resolve - so the header and data rows stay column-aligned.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void Termini_Header_SmNull_WritesBlankCells()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            try
            {
                var dict = new Dictionary<string, string>();
                PsmTsvWriter.AddPeptideSequenceData(dict, null, null);

                Assert.That(dict[SpectrumMatchFromTsvHeader.FivePrimeTerminus], Is.EqualTo(" "));
                Assert.That(dict[SpectrumMatchFromTsvHeader.ThreePrimeTerminus], Is.EqualTo(" "));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// The terminal-cap writer guards each null hop on the way to the formula string. A terminus
        /// whose ThisChemicalFormula is null (never produced by a real digestion, but defensively
        /// supported) must yield an empty cell for that block and leave the other cap untouched.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public static void Termini_Oligo_NullTerminusChemicalFormula_WritesEmptyCell()
        {
            var previousAnalyteType = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            try
            {
                ChemicalFormula realThreePrime = ChemicalFormula.ParseFormula("H2O4P");
                OligoWithSetMods oligo = BuildTestOligo("AAA", new NullFormulaTerminus(), realThreePrime);

                var osm = BuildOsmFromOligo(oligo);

                var headerSplits = SpectralMatch.GetTabSeparatedHeader().Split('\t');
                int fivePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.FivePrimeTerminus);
                int threePrimeIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.ThreePrimeTerminus);
                string[] rowSplits = osm.ToString(new Dictionary<string, int>()).Split('\t');

                Assert.That(fivePrimeIndex, Is.GreaterThanOrEqualTo(0));
                Assert.That(threePrimeIndex, Is.GreaterThanOrEqualTo(0));
                Assert.That(rowSplits.Length, Is.EqualTo(headerSplits.Length));
                // Null ThisChemicalFormula -> null formula string -> Resolve(all null) -> empty cell.
                Assert.That(rowSplits[fivePrimeIndex], Is.EqualTo(string.Empty));
                // The other cap still writes normally.
                Assert.That(rowSplits[threePrimeIndex], Is.EqualTo(realThreePrime.Formula));
            }
            finally
            {
                GlobalVariables.AnalyteType = previousAnalyteType;
            }
        }

        /// <summary>
        /// A terminus object with a monoisotopic mass but no chemical formula, exercising the
        /// writer's ThisChemicalFormula null guard. Never produced by real digestion.
        /// </summary>
        private class NullFormulaTerminus : IHasChemicalFormula
        {
            public double MonoisotopicMass => 0.0;

            public ChemicalFormula ThisChemicalFormula => null;
        }

        #endregion

        #region Collisional Energy Tests

        [Test]
        public static void TestCollisionalEnergy_ParsedFromScan()
        {
            var protein = new Protein("PEPTIDE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, "",
                    hcdEnergy: "42.00"),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);

            Assert.That(psm.CollisionalEnergy, Is.EqualTo(42.0));
        }

        [Test]
        public static void TestCollisionalEnergy_NullWhenNoHcdEnergy()
        {
            var protein = new Protein("PEPTIDE", "TestProtein");
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

            Assert.That(psm.CollisionalEnergy, Is.Null);
        }

        [Test]
        public static void TestCollisionalEnergy_NullWhenInvalidHcdEnergy()
        {
            var protein = new Protein("PEPTIDE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, "",
                    hcdEnergy: "N/A"),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);

            Assert.That(psm.CollisionalEnergy, Is.Null);
        }

        [Test]
        public static void TestCollisionalEnergy_ColumnInOutputWhenFlagSet()
        {
            var protein = new Protein("PEPTIDE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, "",
                    hcdEnergy: "28"),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);

            var headerSplits = SpectralMatch.GetTabSeparatedHeader(includeCollisionalEnergyColumn: true).Split('\t');
            var psmString = psm.ToString(new Dictionary<string, int>(), includeCollisionalEnergyColumn: true);
            var psmSplits = psmString.Split('\t');

            var collisionEnergyIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.CollisionEnergy);
            Assert.That(collisionEnergyIndex, Is.GreaterThanOrEqualTo(0));
            Assert.That(psmSplits[collisionEnergyIndex], Is.EqualTo("28.00"));
        }

        [Test]
        public static void TestCollisionalEnergy_ColumnNotIncludedByDefault()
        {
            var protein = new Protein("PEPTIDE", "TestProtein");
            var digestionParams = new DigestionParams();
            var peptide = new PeptideWithSetModifications(protein, digestionParams, 1, 7,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            double mass = 12.0 + peptide.MonoisotopicMass.ToMz(1);
            var scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                    0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, "",
                    hcdEnergy: "28"),
                mass, 1, "", new CommonParameters());

            var mfi = new List<MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), mfi);

            var headerDefault = SpectralMatch.GetTabSeparatedHeader().Split('\t');
            var headerWithCol = SpectralMatch.GetTabSeparatedHeader(includeCollisionalEnergyColumn: true).Split('\t');

            Assert.That(headerDefault, Does.Not.Contain(SpectrumMatchFromTsvHeader.CollisionEnergy));
            Assert.That(headerWithCol, Does.Contain(SpectrumMatchFromTsvHeader.CollisionEnergy));
        }

        [Test]
        public static void TestCollisionalEnergy_WritesNotAvailable_WhenMissing()
        {
            var protein = new Protein("PEPTIDE", "TestProtein");
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

            var headerSplits = SpectralMatch.GetTabSeparatedHeader(includeCollisionalEnergyColumn: true).Split('\t');
            var psmString = psm.ToString(new Dictionary<string, int>(), includeCollisionalEnergyColumn: true);
            var psmSplits = psmString.Split('\t');

            var collisionEnergyIndex = headerSplits.IndexOf(SpectrumMatchFromTsvHeader.CollisionEnergy);
            Assert.That(collisionEnergyIndex, Is.GreaterThanOrEqualTo(0));
            Assert.That(psmSplits[collisionEnergyIndex], Is.EqualTo("N/A"));
        }

        #endregion
    }
}