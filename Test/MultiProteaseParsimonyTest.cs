using EngineLayer;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

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

            DigestionParams digestionParams_Tryp = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            DigestionParams digestionParams_ArgC = new DigestionParams(protease: "Arg-C", minPeptideLength: 1);

            PeptideWithSetModifications pepA_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 5, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_3T = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams_Tryp, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2A = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams_ArgC, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_3A = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams_ArgC, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmPEPR_T = new PeptideSpectralMatch(pepA_1T, 0, 10, 0, scan, digestionParams_Tryp, new List<MatchedFragmentIon>());
            psmPEPR_T.AddOrReplace(pepA_2T, 10, 0, true, new List<MatchedFragmentIon>(),0);
            psmPEPR_T.AddOrReplace(pepA_3T, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmPEPR_A = new PeptideSpectralMatch(pepA_2A, 0, 10, 0, scan, digestionParams_ArgC, new List<MatchedFragmentIon>());
            psmPEPR_A.AddOrReplace(pepA_3A, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmABCK_T = new PeptideSpectralMatch(pepB_1T, 0, 10, 0, scan, digestionParams_Tryp, new List<MatchedFragmentIon>());

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmPEPR_T, psmPEPR_A, psmABCK_T };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams_Tryp);
            digestionParamsList.Add(digestionParams_ArgC);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.AreEqual(2, proteinGroups.Count);
            var proteinGroup1 = proteinGroups.Where(h => h.ProteinGroupName == "1").First();
            Assert.AreEqual(1, proteinGroup1.UniquePeptides.Count);
            Assert.AreEqual(2, proteinGroup1.AllPeptides.Count);
            var proteinGroup2 = proteinGroups.Where(h => h.ProteinGroupName == "2|3").First();
            Assert.AreEqual(0, proteinGroup2.UniquePeptides.Count);
            Assert.AreEqual(4, proteinGroup2.AllPeptides.Count);
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

            DigestionParams digestionParams = new DigestionParams(protease: protease1.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_1Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            psmABC_Dp1.AddOrReplace(pepA_2Dp1, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepA_1Dp2, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmABC_Dp2.AddOrReplace(pepA_2Dp2, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmEFG_Dp1 = new PeptideSpectralMatch(pepB_2Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>() { psmABC_Dp1, psmABC_Dp2, psmEFG_Dp1 };

            psms.ForEach(p => p.ResolveAllAmbiguities());
            psms.ForEach(p => p.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            // should result in 1 protein group (protein2)
            Assert.AreEqual(1, proteinGroups.Count);
            Assert.AreEqual("2", proteinGroups.ElementAt(0).ProteinGroupName);
        }

        /// <summary>
        /// In this test, we are ensuring that although two peptides may have the same base sequence (ABC) if they result from only a single protein in the "proteome"
        /// when digested with a  protease they should be considered unique.
        /// Therefore, ABC should be a unique peptide for protein 1 with protease A and for protein 2 with protease B.
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_SharedSequenceCanBeUniquePeptide()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
            };

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("-", null, 1, null), new DigestionMotif("G", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease1 = new Protease("proteaseA", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("proteaseB", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            //a hashset of peptideWithSetModifications from all proteases
            var pwsmList = new HashSet<PeptideWithSetModifications>();

            var proteinList = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                proteinList.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease1.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZ", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZ", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepC_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFGABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_2Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepB_1Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepB_2Dp2, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmEFGABC_Dp1 = new PeptideSpectralMatch(pepC_2Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmXYZ_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            psmXYZ_Dp1.AddOrReplace(pepA_2Dp1, 10, 0, true, new List<MatchedFragmentIon>(),0);

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>() { psmABC_Dp1, psmABC_Dp2, psmEFGABC_Dp1, psmXYZ_Dp1 };

            psms.ForEach(p => p.ResolveAllAmbiguities());
            psms.ForEach(p => p.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.AreEqual(2, proteinGroups.Count);

            var proteinGroup1 = proteinGroups.Where(p => p.ProteinGroupName == "1").First();
            Assert.AreEqual(2, proteinGroup1.AllPeptides.Count);
            Assert.AreEqual(1, proteinGroup1.UniquePeptides.Count);
            var pg1pep1 = proteinGroup1.AllPeptides.Where(p => p.BaseSequence == "XYZ").First();
            Assert.That(pg1pep1.DigestionParams.Protease.Name == "proteaseA");
            var pg1pep2 = proteinGroup1.AllPeptides.Where(p => p.BaseSequence == "ABC").First();
            Assert.That(pg1pep2.DigestionParams.Protease.Name == "proteaseA");
            Assert.That(proteinGroup1.UniquePeptides.First().BaseSequence.Equals("ABC"));

            var proteinGroup2 = proteinGroups.Where(p => p.ProteinGroupName == "2").First();
            Assert.AreEqual(3, proteinGroup2.AllPeptides.Count);
            Assert.AreEqual(2, proteinGroup2.UniquePeptides.Count);
            var pg2pep1 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "XYZ").First();
            Assert.That(pg2pep1.DigestionParams.Protease.Name == "proteaseA");
            var pg2pep2 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "ABC").First();
            Assert.That(pg2pep2.DigestionParams.Protease.Name == "proteaseB");
            var pg2pep3 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "EFGABC").First();
            Assert.That(pg2pep3.DigestionParams.Protease.Name == "proteaseA");
            var uniquePeptideSequences = proteinGroup2.UniquePeptides.Select(p => p.BaseSequence).ToList();
            Assert.That(uniquePeptideSequences.Contains("ABC"));
            Assert.That(uniquePeptideSequences.Contains("EFGABC"));
        }

        /// <summary>
        /// In this test, we want to ensure that protein groups that are actually distinguishable becasue of multiprotease data are not being merged.
        /// Without taking into account the protease peptides would result from, these two proteins (1 and2) would have the same peptide base sequences supporting them.
        /// If that was all that was considered for merging indistinguishable protein groups then we would end up with "1|2", but this would be incorrect.
        /// Becasue of multiple proteases, these two proteins are actually distinguishable ABC can only come from protein 1 with protease A and only from protein2 wiht proteaseB.
        /// The same can be said for EFG. Therefore, we should end up with a protein list contianing both protein 1 and protein 2 supported by 2 unique peptides for each!
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_IndistringuishableProteinsNowDistinguishable()
        {
            string[] sequences = {
                "ABCEFG",
                "EFGABC",
            };

            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage = new List<Tuple<string, FragmentationTerminus>> {
                new Tuple<string, FragmentationTerminus>("C", FragmentationTerminus.C) };
            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage2 = new List<Tuple<string, FragmentationTerminus>> {
                new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C) };

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("C", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease = new Protease("testA", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("testB", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepA_2Dp2, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmEFG_Dp1 = new PeptideSpectralMatch(pepB_1Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmEFG_Dp2 = new PeptideSpectralMatch(pepB_2Dp2, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>() { psmABC_Dp1, psmABC_Dp2, psmEFG_Dp1, psmEFG_Dp2 };

            psms.ForEach(h => h.ResolveAllAmbiguities());
            psms.ForEach(h => h.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.AreEqual(2, proteinGroups.Count);

            // check first protein group
            ProteinGroup pg1 = proteinGroups.Where(v => v.ProteinGroupName == "1").First();
            PeptideWithSetModifications pg1pep1 = pg1.AllPeptides.Where(v => v.BaseSequence == "ABC").First();
            PeptideWithSetModifications pg1pep2 = pg1.AllPeptides.Where(v => v.BaseSequence == "EFG").First();
            Assert.That(pg1.UniquePeptides.Contains(pg1pep1));
            Assert.That(pg1pep1.DigestionParams.Protease.Name == "testA");
            Assert.That(pg1.UniquePeptides.Contains(pg1pep2));
            Assert.That(pg1pep2.DigestionParams.Protease.Name == "testA");
            Assert.That(pg1.AllPeptides.Count == 2);
            Assert.That(pg1.UniquePeptides.Count == 2);

            // check second protein group
            ProteinGroup pg2 = proteinGroups.Where(v => v.ProteinGroupName == "2").First();
            PeptideWithSetModifications pg2pep1 = pg2.AllPeptides.Where(v => v.BaseSequence == "ABC").First();
            PeptideWithSetModifications pg2pep2 = pg2.AllPeptides.Where(v => v.BaseSequence == "EFG").First();
            Assert.That(pg2.UniquePeptides.Contains(pg2pep1));
            Assert.That(pg2pep1.DigestionParams.Protease.Name == "testB");
            Assert.That(pg2.UniquePeptides.Contains(pg2pep2));
            Assert.That(pg2pep2.DigestionParams.Protease.Name == "testB");
            Assert.That(pg2.AllPeptides.Count == 2);
            Assert.That(pg2.UniquePeptides.Count == 2);
        }

        /// <summary>
        /// In this test, we are showing that peptide ABC although coming from the same protein, same location can be 2 separate unique peptides
        /// because of the digestion of the two proteases. The resultant Protein 1 protien group should contian 2 unique peptides both of which have sequence ABC.
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_SameAminoAcidsResultInTwoUniquePeptidesForOneProtein()
        {
            string[] sequences = {
                "-XYZ-EFGABC"
            };

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("A", null, 0, null), new DigestionMotif("Z", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease = new Protease("test1", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test2", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            PeptideWithSetModifications pepA = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            var peptideList = new HashSet<PeptideWithSetModifications>
            {
                pepA,
                pepB
            };

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();
            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            // creates psms for specific PWSM
            foreach (var pwsm in peptideList)
            {
                if (pwsm.DigestionParams == digestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(pwsm, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                }
                if (pwsm.DigestionParams == digestionParams2)
                {
                    psms.Add(new PeptideSpectralMatch(pwsm, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>()));
                }
                psms.Last().SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
                psms.Last().ResolveAllAmbiguities();
            }

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.AreEqual(1, proteinGroups.Count);
            Assert.AreEqual("1", proteinGroups.ElementAt(0).ProteinGroupName);
            Assert.AreEqual(2, proteinGroups.ElementAt(0).UniquePeptides.Count);
        }

        /// <summary>
        /// In this test, the peptide sequence ABC  results in a unique peptide for protein 1 when the sample is digested with protease test5.
        /// But when the sample is digested with protease test6 the base sequece ABC is a shared peptide between protein 2 and 3.
        /// The protein list should contain protein 1 and the protein group protein2|protein3
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_TestingPeptideBaseSequenceCanBeBothSharedAndUnique()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
                "-XYZ-GABC"
            };

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("-", null, 1, null), new DigestionMotif("Z", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease = new Protease("test5", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test6", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_3Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepA_2Dp2, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmABC_Dp2.AddOrReplace(pepA_3Dp2, 10, 0, true, new List<MatchedFragmentIon>(),0);

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>() { psmABC_Dp1, psmABC_Dp2 };

            psms.ForEach(h => h.ResolveAllAmbiguities());
            psms.ForEach(h => h.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.AreEqual(2, proteinGroups.Count);
            if (proteinGroups.ElementAt(0).ProteinGroupName == "1")
            {
                Assert.AreEqual("1", proteinGroups.ElementAt(0).ProteinGroupName);
                Assert.AreEqual(1, proteinGroups.ElementAt(0).UniquePeptides.Count);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).FullSequence);
                Assert.AreEqual("2|3", proteinGroups.ElementAt(1).ProteinGroupName);
                Assert.AreEqual(0, proteinGroups.ElementAt(1).UniquePeptides.Count);
                Assert.AreEqual(2, proteinGroups.ElementAt(1).AllPeptides.Count);
            }
            else
            {
                Assert.AreEqual("1", proteinGroups.ElementAt(1).ProteinGroupName);
                Assert.AreEqual(1, proteinGroups.ElementAt(1).UniquePeptides.Count);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(0).FullSequence);
                Assert.AreEqual("2|3", proteinGroups.ElementAt(0).ProteinGroupName);
                Assert.AreEqual(0, proteinGroups.ElementAt(0).UniquePeptides.Count);
                Assert.AreEqual(2, proteinGroups.ElementAt(0).AllPeptides.Count);
            }
        }

        /// <summary>
        /// In this test, like the previous test, ABC base sewuence can be either unique or shared depending on the protease.
        /// Unlike the previous test, only a psm corresponding to the test3 protease digesting, producing the unique peptide occurs.
        /// Therefore, only protein 1 is listed in the protein list. This test ensures in another manner that unique peptide assignment
        /// is operating correctly.
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_BaseSequenceCanBeSharedOrUniqueButOnlyUnqiuePSMSeen()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
                "-XYZ-GABC"
            };

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("-", null, 1, null), new DigestionMotif("Z", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };


            var protease = new Protease("test3", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test4", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            var peptideList = new List<PeptideWithSetModifications>();
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                    }
                }
            }

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();
            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            // creates psms for specific PWSM
            foreach (var pwsm in peptideList)
            {
                if (pwsm.DigestionParams == digestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(pwsm, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                }
                psms.Last().SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
                psms.Last().ResolveAllAmbiguities();
            }

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.AreEqual(1, proteinGroups.Count);
            Assert.AreEqual("1", proteinGroups.ElementAt(0).ProteinGroupName);
            Assert.AreEqual(1, proteinGroups.ElementAt(0).UniquePeptides.Count);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).FullSequence);
        }

        /// <summary>
        /// In this test, generated PSMs are set to various good or bad FDR levels (above or below 1% FDR) and the filterng
        /// of PSMs is tested prior to parsimony. Only PSMs with FDR at or below 1% should make it through to parsimony.
        /// </summary>
        [Test]
        public static void TestPSMFdrFiltering_Simulated()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
            };
            
            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("-", null, 1, null), new DigestionMotif("Z", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease = new Protease("testC", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("testD", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new List<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        case "XYZ": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "-XYZ-EFG": peptideList.Add(peptide); break;
                    }
                }
            }

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();
            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            double goodFdr = 0.00100;
            double badFdr = 0.0200;

            // creates psms for specific PWSM
            for (int i = 0; i < peptideList.Count(); i++)
            {
                if (peptideList[i].DigestionParams == digestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(peptideList[i], 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                }
                if (peptideList[i].DigestionParams == digestionParams2)
                {
                    psms.Add(new PeptideSpectralMatch(peptideList[i], 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>()));
                }

                switch (i)
                {
                    case 0:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, badFdr, 0, 0, 0, false);
                        break;

                    case 1:
                        psms.Last().SetFdrValues(0, 0, badFdr, 0, 0, goodFdr, 0, 0, 0, false);
                        break;

                    case 2:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
                        break;

                    case 3:
                        psms.Last().SetFdrValues(0, 0, badFdr, 0, 0, badFdr, 0, 0, 0, false);
                        break;

                    case 4:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
                        break;

                    case 5:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
                        break;
                }
                psms.Last().ResolveAllAmbiguities();
            }

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            psms.ElementAt(5).SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);

            //this iscopy of code that filteres psms in PostSearch Analysis Task
            var fdrFilteredPsms = new List<PeptideSpectralMatch>();
            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm != null && psm.FdrInfo.QValue <= 0.0100 && psm.FdrInfo.QValueNotch <= 0.0100)
                {
                    fdrFilteredPsms.Add(psm);
                }
            }

            Assert.AreEqual(3, fdrFilteredPsms.Count);

            var test1 = fdrFilteredPsms.Contains(psms.ElementAt(2));
            var test2 = fdrFilteredPsms.Contains(psms.ElementAt(4));
            var test3 = fdrFilteredPsms.Contains(psms.ElementAt(5));
            var test4 = fdrFilteredPsms.Contains(psms.ElementAt(0));
            var test5 = fdrFilteredPsms.Contains(psms.ElementAt(1));
            var test6 = fdrFilteredPsms.Contains(psms.ElementAt(3));
            Assert.AreEqual(true, test1);
            Assert.AreEqual(true, test2);
            Assert.AreEqual(true, test3);
            Assert.AreEqual(false, test4);
            Assert.AreEqual(false, test5);
            Assert.AreEqual(false, test6);
        }

        /// <summary>
        /// This test again check the PSM FDR filtering capabilities prior to parsimony, unlike the previous test,
        /// this test utilizes a sample file with real speectra generating FDR from target decoy analysis instead of simulated values.
        /// This test actually checks the code inside the post search analysis test.
        /// </summary>

        [Test]
        public static void TestPSMFdrFiltering_RealFile()
        {
            SearchTask Task1 = new SearchTask
            {
                CommonParameters = new CommonParameters
                (
                    qValueOutputFilter: 1
                ),

                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchTarget = true,
                    WritePrunedDatabase = true,
                    SearchType = SearchType.Classic
                }
            };
            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPSMFdrFiltering_RealFileTest");

            var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("TestPSMFdrFiltering_RealFile", Task1) }, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            var thisTaskOutputFolder = Path.Combine(MySetUpClass.outputFolder, @"TestPSMFdrFiltering_RealFile");

            var psms = Path.Combine(thisTaskOutputFolder, "AllPSMs.psmtsv");

            Assert.AreEqual(12, File.ReadLines(psms).Count());
            var protGroups = Path.Combine(thisTaskOutputFolder, "AllProteinGroups.tsv");

            Assert.AreEqual(7, File.ReadLines(protGroups).Count());
            Directory.Delete(outputFolder, true);
        }

        /// <summary>
        /// In this test, the peptide sequence ABC  results in a unique peptide for protein 1 when the sample is digested with protease alpha.
        /// But when the sample is digested with protease beta the base sequece ABC is a shared peptide between protein 2 and 4.
        /// Peptide EFG is shared between protein 3 and 4. This is a more complex testing set to ensure that Parsing of shared peptides when unique proteins
        /// are present is being handled correctly.
        /// The protein list should contain protein 1 and the protein 4.
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_TestingSameBaseSequenceSharedandUniqueMoreComplexSample()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-XYGABC",
                "-ABGEFG-XYZ",
                "-XYZ-GEFGABC",
            };
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("-", null, 1, null), new DigestionMotif("Z", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease = new Protease("proteaseAlpha", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("proteaseBeta", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            PeptideWithSetModifications pepABC_1Alpha = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABC_2Beta = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABC_4Beta = new PeptideWithSetModifications(protein: p.ElementAt(3), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepEFG_4Beta = new PeptideWithSetModifications(protein: p.ElementAt(3), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepEFG_3Beta = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 5, oneBasedEndResidueInProtein: 7, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmABC_Alpha = new PeptideSpectralMatch(pepABC_1Alpha, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmABC_Beta = new PeptideSpectralMatch(pepABC_2Beta, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmABC_Beta.AddOrReplace(pepABC_4Beta, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmEFG_Beta = new PeptideSpectralMatch(pepEFG_3Beta, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmEFG_Beta.AddOrReplace(pepEFG_4Beta, 10, 0, true, new List<MatchedFragmentIon>(),0);

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmABC_Alpha, psmABC_Beta, psmEFG_Beta };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.AreEqual(2, proteinGroups.Count);
            var proteinGroup1 = proteinGroups.Where(h => h.ProteinGroupName == "1").First();
            Assert.AreEqual(1, proteinGroup1.UniquePeptides.Count);
            Assert.AreEqual(1, proteinGroup1.AllPeptides.Count);
            var proteinGroup2 = proteinGroups.Where(h => h.ProteinGroupName == "4").First();
            Assert.AreEqual(0, proteinGroup2.UniquePeptides.Count);
            Assert.AreEqual(2, proteinGroup2.AllPeptides.Count);
        }

        /// <summary>
        /// This test is to ensure that proteins, even with multiprotease are truly indistinguishable.
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_TestingActuallyIndistinguisableProteins()
        {
            string[] sequences = {
                "KABCKXYZK",
                "KXYZKABCK",
            };
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: "Lys-C (don't cleave before proline)", minPeptideLength: 1);

            PeptideWithSetModifications pepABCK_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABCK_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABCK_1L = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABCK_2L = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_1L = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_2L = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmABCK_T = new PeptideSpectralMatch(pepABCK_1T, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            psmABCK_T.AddOrReplace(pepABCK_2T, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmABCK_L = new PeptideSpectralMatch(pepABCK_1L, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmABCK_L.AddOrReplace(pepABCK_2L, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmXYZK_T = new PeptideSpectralMatch(pepXYZK_1T, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmXYZK_T.AddOrReplace(pepXYZK_2T, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmXYZK_L = new PeptideSpectralMatch(pepXYZK_1L, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmXYZK_L.AddOrReplace(pepXYZK_2L, 10, 0, true, new List<MatchedFragmentIon>(),0);

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmABCK_T, psmABCK_L, psmXYZK_T, psmXYZK_L };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.AreEqual(1, proteinGroups.Count);
            Assert.AreEqual("1|2", proteinGroups.ElementAt(0).ProteinGroupName);
            Assert.AreEqual(8, proteinGroups.ElementAt(0).AllPeptides.Count);
        }

        /// <summary>
        /// This test ensures that having the main portion of the greedy algorithm be agnostic of protease is valid, and ensures that when needed the algorithm uses protease to break ties
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_TestingGreedyAlgorithm()
        {
            string[] sequences = {
                "-ABC-XYZ",
                "-ABC-GABC",
                "-XYZ-GABC"
            };
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            
            List<DigestionMotif> motifs1 = new List<DigestionMotif> { new DigestionMotif("-", null, 0, null), new DigestionMotif("-", null, 1, null) };
            List<DigestionMotif> motifs2 = new List<DigestionMotif> { new DigestionMotif("G", null, 1, null) };

            var protease = new Protease("proteaseDash", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("proteaseG", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            PeptideWithSetModifications pepABC_1Dash = new PeptideWithSetModifications(p.ElementAt(0), digestionParams, 2, 4, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pepABC_2Dash = new PeptideWithSetModifications(p.ElementAt(1), digestionParams, 7, 9, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pepABC_2G = new PeptideWithSetModifications(p.ElementAt(1), digestionParams2, 2, 4, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pepABC_3G = new PeptideWithSetModifications(p.ElementAt(2), digestionParams2, 7, 9, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmABC_Dash = new PeptideSpectralMatch(pepABC_1Dash, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            psmABC_Dash.AddOrReplace(pepABC_2Dash, 10, 0, true, new List<MatchedFragmentIon>(),0);
            PeptideSpectralMatch psmABC_G = new PeptideSpectralMatch(pepABC_2G, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>());
            psmABC_G.AddOrReplace(pepABC_3G, 10, 0, true, new List<MatchedFragmentIon>(),0);

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmABC_Dash, psmABC_G };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, double.NaN, double.NaN, double.NaN, false));

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.AreEqual(1, proteinGroups.Count);
            Assert.AreEqual("2", proteinGroups.ElementAt(0).ProteinGroupName);
            Assert.AreEqual(2, proteinGroups.ElementAt(0).AllPeptides.Count);
        }

        /// <summary>
        /// This test ensures that FDR for each psm is calculated accoriding to its protease
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_TestingProteaseSpecificFDRCalculations()
        {
            // two protease options
            DigestionParams tryp = new DigestionParams(protease: "trypsin");
            DigestionParams gluC = new DigestionParams(protease: "Glu-C");

            // target or decoy protein
            Protein t = new Protein("P", "1");
            Protein d = new Protein("P", "2", isDecoy: true);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[1], new double[1], false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 0, null);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());
            List<MatchedFragmentIon> f = new List<MatchedFragmentIon>();

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>
            {
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: tryp, p: t), 0, 20, 1, scan, tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: gluC, p: t), 0, 19, 1, scan, gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: tryp, p: t), 0, 18, 1, scan, tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: gluC, p: t), 0, 17, 1, scan, gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: gluC, p: d), 0, 16, 1, scan, gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: gluC, p: t), 0, 15, 1, scan, gluC, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: tryp, p: t), 0, 14, 1, scan, tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: tryp, p: d), 0, 13, 1, scan, tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: tryp, p: d), 0, 12, 1, scan, tryp, f),
                new PeptideSpectralMatch(new PeptideWithSetModifications("P", null, digestionParams: tryp, p: t), 0, 11, 1, scan, tryp, f),
            };

            psms.ForEach(p => p.ResolveAllAmbiguities());

            new FdrAnalysisEngine(psms, 0, new CommonParameters(), new List<string>()).Run();
            psms = psms.OrderByDescending(p => p.Score).ToList();

            Assert.AreEqual(0.00, Math.Round(psms[0].FdrInfo.QValue, 2));
            Assert.AreEqual(0.00, Math.Round(psms[1].FdrInfo.QValue, 2));
            Assert.AreEqual(0.00, Math.Round(psms[2].FdrInfo.QValue, 2));
            Assert.AreEqual(0.00, Math.Round(psms[3].FdrInfo.QValue, 2));
            Assert.AreEqual(0.33, Math.Round(psms[4].FdrInfo.QValue, 2));
            Assert.AreEqual(0.33, Math.Round(psms[5].FdrInfo.QValue, 2));
            Assert.AreEqual(0.00, Math.Round(psms[6].FdrInfo.QValue, 2));
            Assert.AreEqual(0.33, Math.Round(psms[7].FdrInfo.QValue, 2));
            Assert.AreEqual(0.50, Math.Round(psms[8].FdrInfo.QValue, 2));
            Assert.AreEqual(0.50, Math.Round(psms[9].FdrInfo.QValue, 2));
        }
    }
}