using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using Omics.Digestion;
using Omics.Modifications;
using TaskLayer;
using Omics;

namespace Test
{
    [TestFixture]
    public static class MultiProteaseParsimonyTest
    {
        /// <summary>
        /// three proteins proteases will be trypsina and argC two proteins will produce the peptide by both arge C and trypsin (cleave at R)
        /// the  last one can only make it by trypsin. There is unique peptide also for the trypsin only protein.
        /// </summary>
        [Test]
        public static void MultiProteaseTest()
        {
            string[] sequences = {
                "ABCKPEPR",
                "BRPEPR",
                "ARPEPR"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            CommonParameters commonParameters_ArgC = new CommonParameters(digestionParams: new DigestionParams(protease: "Arg-C", minPeptideLength: 1));

            PeptideWithSetModifications pepA_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 5, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_3T = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2A = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_ArgC.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_3A = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_ArgC.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmPEPR_T = new PeptideSpectralMatch(pepA_1T, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmPEPR_T.AddOrReplace(pepA_2T, 10, 0, true, new List<MatchedFragmentIon>());
            psmPEPR_T.AddOrReplace(pepA_3T, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmPEPR_A = new PeptideSpectralMatch(pepA_2A, 0, 10, 0, scan, commonParameters_ArgC, new List<MatchedFragmentIon>());
            psmPEPR_A.AddOrReplace(pepA_3A, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmABCK_T = new PeptideSpectralMatch(pepB_1T, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());

            List<SpectralMatch> psms = new List<SpectralMatch> { psmPEPR_T, psmPEPR_A, psmABCK_T };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));
            FilteredPsms filteredPsms = FilteredPsms.Filter(psms, new CommonParameters(), includeDecoys: true, includeContaminants: true, includeAmbiguous: true, includeAmbiguousMods: true, includeHighQValuePsms: true);

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters_Tryp.DigestionParams);
            digestionParamsList.Add(commonParameters_ArgC.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(filteredPsms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, filteredPsms.FilteredPsmsList, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.That(proteinGroups.Count, Is.EqualTo(2));
            var proteinGroup1 = proteinGroups.Where(h => h.ProteinGroupName == "1").First();
            Assert.That(proteinGroup1.UniquePeptides.Count, Is.EqualTo(1));
            Assert.That(proteinGroup1.AllPeptides.Count, Is.EqualTo(2));
            var proteinGroup2 = proteinGroups.Where(h => h.ProteinGroupName == "2|3").First();
            Assert.That(proteinGroup2.UniquePeptides.Count, Is.EqualTo(0));
            Assert.That(proteinGroup2.AllPeptides.Count, Is.EqualTo(4));
        }

        /// <summary>
        /// In this test, we are simulating a single peptide (pepA) with the same base sequence coming from two different protease digestions.
        /// The proteases cleave at the same spots (mimicking the overlap of peptides between Trypsin and either ArgC or LysC).
        /// So pepA is a shared peptide between protein 1 and protein 2 for both proteases.
        /// With only pepA being observed a protein group containign protien1|protein2 would exist.
        /// If for only one protease we observe a unique peptide (pepB) for protein 2 we would then know for certain that protein2 exists in our sample.
        /// </summary>
        [Test]
        public static void MultiProteaseSamePeptideSameProteinsDifferentProteases()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFG-ABC",
            };
            //both proteases are cleaving at the same spots to simulate trypsin and argC producing the same peptides

            List<DigestionMotif> motifs = new List<DigestionMotif>
            {
                new DigestionMotif("-", null, 1, null),
                new DigestionMotif("Z", null, 1, null),
            };
            var protease1 = new Protease("proteaseA1", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("proteaseB1", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            //a hashset of peptideWithSetModifications from all proteases
            var pwsmList = new HashSet<PeptideWithSetModifications>();

            var proteinList = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                proteinList.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters1 = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: commonParameters1.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters1.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters1.DigestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_1Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, commonParameters1, new List<MatchedFragmentIon>());
            psmABC_Dp1.AddOrReplace(pepA_2Dp1, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepA_1Dp2, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            psmABC_Dp2.AddOrReplace(pepA_2Dp2, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmEFG_Dp1 = new PeptideSpectralMatch(pepB_2Dp1, 0, 10, 0, scan, commonParameters1, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>() { psmABC_Dp1, psmABC_Dp2, psmEFG_Dp1 };

            psms.ForEach(p => p.ResolveAllAmbiguities());
            psms.ForEach(p => p.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));
            FilteredPsms filteredPsms = FilteredPsms.Filter(psms, new CommonParameters(), includeDecoys: true, includeContaminants: true, includeAmbiguous: true, includeAmbiguousMods: true, includeHighQValuePsms: true);

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters1.DigestionParams);
            digestionParamsList.Add(commonParameters2.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(filteredPsms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, filteredPsms.FilteredPsmsList, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            // should result in 1 protein group (protein2)
            Assert.That(proteinGroups.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("2"));
        }

        /// <summary>
        /// In this test, we are ensuring that although two peptides may have the same base sequence (ABC) if they result from only a single protein in the "proteome"
        /// when digested with a protease they should be considered unique.
        /// Therefore, ABC should be a unique peptide for protein 1 with protease A and for protein 2 with protease B.
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_SharedSequenceCanBeUniquePeptide()
        {
            string[] sequences = {
                "-XYZ--ABC",      // Protein 1: 9 residues
                "-XYZ-EFGABC",    // Protein 2: 11 residues
            };

            // Protease A cleaves after - and G
            // Protease B cleaves after G only
            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("-", null, 1, null), new DigestionMotif("G", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease1 = new Protease("proteaseA", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("proteaseB", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            var proteinList = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                proteinList.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            
            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            // Protein 1 with Protease A: -XYZ--ABC -> cleaves after - and G -> XYZ (2-4), ABC (7-9)
            // Protein 2 with Protease A: -XYZ-EFGABC -> cleaves after - and G -> XYZ (2-4), EFG (6-8), ABC (9-11)
            // Protein 1 with Protease B: -XYZ--ABC -> cleaves after G only -> no G, so no internal cleavage
            // Protein 2 with Protease B: -XYZ-EFGABC -> cleaves after G -> ABC (9-11)

            // XYZ peptides - shared between protein 1 and 2 with protease A
            PeptideWithSetModifications pepXYZ_P1_proteaseA = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZ", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZ_P2_proteaseA = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZ", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            
            // ABC from Protein 1 with Protease A - unique to protein 1 for this protease
            PeptideWithSetModifications pepABC_P1_proteaseA = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            
            // EFGABC from Protein 2 with Protease A - unique to protein 2
            PeptideWithSetModifications pepEFGABC_P2_proteaseA = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFGABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            
            // ABC from Protein 2 with Protease B - unique to protein 2 for this protease
            PeptideWithSetModifications pepABC_P2_proteaseB = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            // PSM for XYZ - shared between protein 1 and 2 with protease A
            SpectralMatch psmXYZ_proteaseA = new PeptideSpectralMatch(pepXYZ_P1_proteaseA, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            psmXYZ_proteaseA.AddOrReplace(pepXYZ_P2_proteaseA, 10, 0, true, new List<MatchedFragmentIon>());
            
            // PSM for ABC - unique to protein 1 with protease A
            SpectralMatch psmABC_proteaseA = new PeptideSpectralMatch(pepABC_P1_proteaseA, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            
            // PSM for EFGABC - unique to protein 2 with protease A
            SpectralMatch psmEFGABC_proteaseA = new PeptideSpectralMatch(pepEFGABC_P2_proteaseA, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            
            // PSM for ABC - unique to protein 2 with protease B
            SpectralMatch psmABC_proteaseB = new PeptideSpectralMatch(pepABC_P2_proteaseB, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>() { psmXYZ_proteaseA, psmABC_proteaseA, psmEFGABC_proteaseA, psmABC_proteaseB };

            psms.ForEach(p => p.ResolveAllAmbiguities());
            psms.ForEach(p => p.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0));
            FilteredPsms filteredPsms = FilteredPsms.Filter(psms, new CommonParameters(), includeDecoys: true, includeContaminants: true, includeAmbiguous: true, includeAmbiguousMods: true, includeHighQValuePsms: true);

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(filteredPsms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, filteredPsms.FilteredPsmsList, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            
            // Should result in 2 protein groups
            Assert.That(proteinGroups.Count, Is.EqualTo(2));

            var proteinGroup1 = proteinGroups.Where(p => p.ProteinGroupName == "1").First();
            Assert.That(proteinGroup1.AllPeptides.Count, Is.EqualTo(2)); // XYZ and ABC from proteaseA
            Assert.That(proteinGroup1.UniquePeptides.Count, Is.EqualTo(1)); // ABC from proteaseA is unique
            Assert.That(proteinGroup1.UniquePeptides.First().BaseSequence, Is.EqualTo("ABC"));
            Assert.That(proteinGroup1.UniquePeptides.First().DigestionParams.DigestionAgent.Name, Is.EqualTo("proteaseA"));

            var proteinGroup2 = proteinGroups.Where(p => p.ProteinGroupName == "2").First();
            Assert.That(proteinGroup2.AllPeptides.Count, Is.EqualTo(3)); // XYZ from proteaseA, EFGABC from proteaseA, ABC from proteaseB
            Assert.That(proteinGroup2.UniquePeptides.Count, Is.EqualTo(2)); // EFGABC from proteaseA and ABC from proteaseB are unique
            var uniquePeptideSequences = proteinGroup2.UniquePeptides.Select(p => p.BaseSequence).ToList();
            Assert.That(uniquePeptideSequences.Contains("ABC"));
            Assert.That(uniquePeptideSequences.Contains("EFGABC"));
        }

        /// <summary>
        /// In this test, we want to ensure that protein groups that are actually distinguishable because of multiprotease data are not being merged.
        /// Without taking into account the protease peptides would result from, these two proteins (1 and 2) would have the same peptide base sequences supporting them.
        /// If that was all that was considered for merging indistinguishable protein groups then we would end up with "1|2", but this would be incorrect.
        /// Because of multiple proteases, these two proteins are actually distinguishable: ABC can only come from protein 1 with protease A and only from protein 2 with protease B.
        /// The same can be said for EFG. Therefore, we should end up with a protein list containing both protein 1 and protein 2 supported by 2 unique peptides for each!
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_IndistringuishableProteinsNowDistinguishable()
        {
            string[] sequences = {
                "ABCEFG",
                "EFGABC",
            };

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("C", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease1 = new Protease("testA", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("testB", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            // Protease A (testA) cleaves after C:
            // Protein 1 (ABCEFG): ABC (1-3), EFG (4-6) - both unique to protein 1 for this protease
            // Protein 2 (EFGABC): EFGABC (no internal C until position 6) - no ABC or EFG peptides produced separately
            
            // Protease B (testB) cleaves after G:
            // Protein 1 (ABCEFG): ABCEFG (cleaves after G at position 6) - no separate peptides
            // Protein 2 (EFGABC): EFG (1-3), ABC (4-6) - both unique to protein 2 for this protease

            // Peptides from Protein 1 with testA (cleaves after C)
            PeptideWithSetModifications pep_ABC_P1_testA = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pep_EFG_P1_testA = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            // Peptides from Protein 2 with testB (cleaves after G)
            PeptideWithSetModifications pep_EFG_P2_testB = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pep_ABC_P2_testB = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            // PSM for ABC - only from Protein 1 with testA (unique peptide for protein 1)
            SpectralMatch psmABC_testA = new PeptideSpectralMatch(pep_ABC_P1_testA, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            
            // PSM for EFG - only from Protein 1 with testA (unique peptide for protein 1)
            SpectralMatch psmEFG_testA = new PeptideSpectralMatch(pep_EFG_P1_testA, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            
            // PSM for ABC - only from Protein 2 with testB (unique peptide for protein 2)
            SpectralMatch psmABC_testB = new PeptideSpectralMatch(pep_ABC_P2_testB, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            
            // PSM for EFG - only from Protein 2 with testB (unique peptide for protein 2)
            SpectralMatch psmEFG_testB = new PeptideSpectralMatch(pep_EFG_P2_testB, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>() { psmABC_testA, psmEFG_testA, psmABC_testB, psmEFG_testB };

            psms.ForEach(p => p.ResolveAllAmbiguities());
            psms.ForEach(p => p.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0));

            FilteredPsms filteredPsms = FilteredPsms.Filter(psms, commonParameters, includeDecoys: true, includeContaminants: true, includeAmbiguous: true, includeAmbiguousMods: true, includeHighQValuePsms: true);
            
            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(filteredPsms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, filteredPsms.FilteredPsmsList, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            
            // Should result in 2 protein groups - protein 1 and protein 2, each with 2 unique peptides
            Assert.That(proteinGroups.Count, Is.EqualTo(2));
            
            var proteinGroup1 = proteinGroups.Where(pg => pg.ProteinGroupName == "1").First();
            Assert.That(proteinGroup1.AllPeptides.Count, Is.EqualTo(2));
            Assert.That(proteinGroup1.UniquePeptides.Count, Is.EqualTo(2));
            Assert.That(proteinGroup1.AllPeptides.Select(pep => pep.BaseSequence).Contains("ABC"));
            Assert.That(proteinGroup1.AllPeptides.Select(pep => pep.BaseSequence).Contains("EFG"));
            Assert.That(proteinGroup1.AllPeptides.All(pep => pep.DigestionParams.DigestionAgent.Name == "testA"));

            var proteinGroup2 = proteinGroups.Where(pg => pg.ProteinGroupName == "2").First();
            Assert.That(proteinGroup2.AllPeptides.Count, Is.EqualTo(2));
            Assert.That(proteinGroup2.UniquePeptides.Count, Is.EqualTo(2));
            Assert.That(proteinGroup2.AllPeptides.Select(pep => pep.BaseSequence).Contains("ABC"));
            Assert.That(proteinGroup2.AllPeptides.Select(pep => pep.BaseSequence).Contains("EFG"));
            Assert.That(proteinGroup2.AllPeptides.All(pep => pep.DigestionParams.DigestionAgent.Name == "testB"));
        }

        /// <summary>
        /// In this test, we are showing that peptide ABC although coming from the same protein, same location can be 2 separate unique peptides
        /// because of the digestion of the two proteases. The resultant Protein 1 protein group should contain 2 unique peptides both of which have sequence ABC.
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_SameAminoAcidsResultInTwoUniquePeptidesForOneProtein()
        {
            string[] sequences = {
                "-XYZ-EFGABC"  // 11 residues
            };

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("A", null, 0, null), new DigestionMotif("Z", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease1 = new Protease("test1", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("test2", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            // ABC peptide from protein 1 with test1 (positions 9-11)
            PeptideWithSetModifications pepA = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            
            // ABC peptide from protein 1 with test2 (positions 9-11)
            PeptideWithSetModifications pepB = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            var peptideList = new HashSet<PeptideWithSetModifications>
            {
                pepA,
                pepB
            };

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>();
            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            // creates psms for specific PWSM
            foreach (var pwsm in peptideList)
            {
                if (pwsm.DigestionParams == commonParameters.DigestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(pwsm, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>()));
                }
                if (pwsm.DigestionParams == commonParameters2.DigestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(pwsm, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>()));
                }
                psms.Last().SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                psms.Last().ResolveAllAmbiguities();
            }
            
            FilteredPsms filteredPsms = FilteredPsms.Filter(psms, new CommonParameters(), includeDecoys: true, includeContaminants: true, includeAmbiguous: true, includeAmbiguousMods: true, includeHighQValuePsms: true);
            
            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(filteredPsms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, filteredPsms.FilteredPsmsList, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            // Should result in 1 protein group with 2 unique peptides (both ABC, but from different proteases)
            Assert.That(proteinGroups.Count, Is.EqualTo(1));

            var proteinGroup1 = proteinGroups.Where(pg => pg.ProteinGroupName == "1").First();
            Assert.That(proteinGroup1.AllPeptides.Count, Is.EqualTo(2));
            Assert.That(proteinGroup1.UniquePeptides.Count, Is.EqualTo(2));
            Assert.That(proteinGroup1.AllPeptides.All(pep => pep.BaseSequence == "ABC"));
            
            // Verify that the two unique peptides come from different proteases
            var proteaseNames = proteinGroup1.UniquePeptides.Select(pep => pep.DigestionParams.DigestionAgent.Name).ToList();
            Assert.That(proteaseNames.Contains("test1"));
            Assert.That(proteaseNames.Contains("test2"));
        }

        /// <summary>
        /// In this test, we are ensuring that FDR calculation at the PSM level is working when the same spearman rank peptide level FDR is assigned for all peptides
        /// Takes advantage of the fact that in the second example, the only difference is that the second PSM has a lower score (awarded the same FDR rank)
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_TestingSpearmanRankFDR()
        {
            CommonParameters commonParameters_tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin"));
            CommonParameters commonParameters_gluC = new CommonParameters(digestionParams: new DigestionParams(protease: "Glu-C"));

            // target or decoy protein
            Protein t = new Protein("P", "1");
            Protein d = new Protein("P", "2", isDecoy: true);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[1], new double[1], false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 0, null);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());
            List<MatchedFragmentIon> f = new List<MatchedFragmentIon>();

            List<SpectralMatch> psms = new List<SpectralMatch>
            {
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_tryp.DigestionParams, p: t), 0, 20, 1, scan, commonParameters_tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_gluC.DigestionParams, p: t), 0, 19, 1, scan, commonParameters_gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_tryp.DigestionParams, p: t), 0, 18, 1, scan, commonParameters_tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_gluC.DigestionParams, p: t), 0, 17, 1, scan, commonParameters_gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_gluC.DigestionParams, p: d), 0, 16, 1, scan, commonParameters_gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_gluC.DigestionParams, p: t), 0, 15, 1, scan, commonParameters_gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_tryp.DigestionParams, p: t), 0, 14, 1, scan, commonParameters_tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_tryp.DigestionParams, p: d), 0, 13, 1, scan, commonParameters_tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_tryp.DigestionParams, p: d), 0, 12, 1, scan, commonParameters_tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: commonParameters_tryp.DigestionParams, p: t), 0, 11, 1, scan, commonParameters_tryp, f),
            };

            psms.ForEach(p => p.ResolveAllAmbiguities());

            List<(string fileName, CommonParameters fileSpecificParameters)> fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("filename", new CommonParameters()) };

            new FdrAnalysisEngine(psms, 0, new CommonParameters(), fsp, new List<string>()).Run();
            psms = psms.OrderByDescending(p => p).ToList();

            //q-value is computed   as targetCount / (decoyCount + targetCount) for each protease separately
            //once a higher q-value is found, it is used for all subsequent PSMs with the same protease even if increasing number of targets would lower the q-value

            //	Row	t/d	score	protease	targetCount	decoyCount	q-value
            //	0	t	20	tryp	1	0	0
            //	1	t	19	gluC	1	0	0
            //	2	t	18	tryp	2	0	0
            //	3	t	17	gluC	2	0	0
            //	4	d	16	gluC	2	1	0.5
            //	5	t	15	gluC	3	1	0.5
            //	6	t	14	tryp	3	0	0
            //	7	d	13	tryp	3	1	0.333333333
            //	8	d	12	tryp	3	2	0.666666667
            //	9	t	11	tryp	4	2	0.666666667

            Assert.That(Math.Round(psms[0].FdrInfo.QValue, 2), Is.EqualTo(0.00));
            Assert.That(Math.Round(psms[1].FdrInfo.QValue, 2), Is.EqualTo(0.00));
            Assert.That(Math.Round(psms[2].FdrInfo.QValue, 2), Is.EqualTo(0.00));
            Assert.That(Math.Round(psms[3].FdrInfo.QValue, 2), Is.EqualTo(0.00));
            Assert.That(Math.Round(psms[4].FdrInfo.QValue, 2), Is.EqualTo(0.50));
            Assert.That(Math.Round(psms[5].FdrInfo.QValue, 2), Is.EqualTo(0.50));
            Assert.That(Math.Round(psms[6].FdrInfo.QValue, 2), Is.EqualTo(0.00));
            Assert.That(Math.Round(psms[7].FdrInfo.QValue, 2), Is.EqualTo(0.33));
            Assert.That(Math.Round(psms[8].FdrInfo.QValue, 2), Is.EqualTo(0.67));
            Assert.That(Math.Round(psms[9].FdrInfo.QValue, 2), Is.EqualTo(0.67));

        }
    }
}