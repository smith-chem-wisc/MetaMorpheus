﻿using EngineLayer;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters_Tryp.DigestionParams);
            digestionParamsList.Add(commonParameters_ArgC.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
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

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters1.DigestionParams);
            digestionParamsList.Add(commonParameters2.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            // should result in 1 protein group (protein2)
            Assert.That(proteinGroups.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("2"));
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
            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZ", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZ", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_1Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepC_2Dp1 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFGABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_2Dp2 = new PeptideWithSetModifications(protein: proteinList.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepB_1Dp1, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            SpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepB_2Dp2, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            SpectralMatch psmEFGABC_Dp1 = new PeptideSpectralMatch(pepC_2Dp1, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            SpectralMatch psmXYZ_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            psmXYZ_Dp1.AddOrReplace(pepA_2Dp1, 10, 0, true, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>() { psmABC_Dp1, psmABC_Dp2, psmEFGABC_Dp1, psmXYZ_Dp1 };

            psms.ForEach(p => p.ResolveAllAmbiguities());
            psms.ForEach(p => p.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters.DigestionParams);
            digestionParamsList.Add(commonParameters2.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.That(proteinGroups.Count, Is.EqualTo(2));

            var proteinGroup1 = proteinGroups.Where(p => p.ProteinGroupName == "1").First();
            Assert.That(proteinGroup1.AllPeptides.Count, Is.EqualTo(2));
            Assert.That(proteinGroup1.UniquePeptides.Count, Is.EqualTo(1));
            var pg1pep1 = proteinGroup1.AllPeptides.Where(p => p.BaseSequence == "XYZ").First();
            Assert.That(pg1pep1.DigestionParams.DigestionAgent.Name == "proteaseA");
            var pg1pep2 = proteinGroup1.AllPeptides.Where(p => p.BaseSequence == "ABC").First();
            Assert.That(pg1pep2.DigestionParams.DigestionAgent.Name == "proteaseA");
            Assert.That(proteinGroup1.UniquePeptides.First().BaseSequence.Equals("ABC"));

            var proteinGroup2 = proteinGroups.Where(p => p.ProteinGroupName == "2").First();
            Assert.That(proteinGroup2.AllPeptides.Count, Is.EqualTo(3));
            Assert.That(proteinGroup2.UniquePeptides.Count, Is.EqualTo(2));
            var pg2pep1 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "XYZ").First();
            Assert.That(pg2pep1.DigestionParams.DigestionAgent.Name == "proteaseA");
            var pg2pep2 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "ABC").First();
            Assert.That(pg2pep2.DigestionParams.DigestionAgent.Name == "proteaseB");
            var pg2pep3 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "EFGABC").First();
            Assert.That(pg2pep3.DigestionParams.DigestionAgent.Name == "proteaseA");
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

            var protease1 = new Protease("testA", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("testB", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 4, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepB_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            SpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepA_2Dp2, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            SpectralMatch psmEFG_Dp1 = new PeptideSpectralMatch(pepB_1Dp1, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            SpectralMatch psmEFG_Dp2 = new PeptideSpectralMatch(pepB_2Dp2, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>() { psmABC_Dp1, psmABC_Dp2, psmEFG_Dp1, psmEFG_Dp2 };

            psms.ForEach(h => h.ResolveAllAmbiguities());
            psms.ForEach(h => h.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.That(proteinGroups.Count, Is.EqualTo(2));

            // check first protein group
            ProteinGroup pg1 = proteinGroups.Where(v => v.ProteinGroupName == "1").First();
            IBioPolymerWithSetMods pg1pep1 = pg1.AllPeptides.Where(v => v.BaseSequence == "ABC").First();
            IBioPolymerWithSetMods pg1pep2 = pg1.AllPeptides.Where(v => v.BaseSequence == "EFG").First();
            Assert.That(pg1.UniquePeptides.Contains(pg1pep1));
            Assert.That(pg1pep1.DigestionParams.DigestionAgent.Name == "testA");
            Assert.That(pg1.UniquePeptides.Contains(pg1pep2));
            Assert.That(pg1pep2.DigestionParams.DigestionAgent.Name == "testA");
            Assert.That(pg1.AllPeptides.Count == 2);
            Assert.That(pg1.UniquePeptides.Count == 2);

            // check second protein group
            ProteinGroup pg2 = proteinGroups.Where(v => v.ProteinGroupName == "2").First();
            IBioPolymerWithSetMods pg2pep1 = pg2.AllPeptides.Where(v => v.BaseSequence == "ABC").First();
            IBioPolymerWithSetMods pg2pep2 = pg2.AllPeptides.Where(v => v.BaseSequence == "EFG").First();
            Assert.That(pg2.UniquePeptides.Contains(pg2pep1));
            Assert.That(pg2pep1.DigestionParams.DigestionAgent.Name == "testB");
            Assert.That(pg2.UniquePeptides.Contains(pg2pep2));
            Assert.That(pg2pep2.DigestionParams.DigestionAgent.Name == "testB");
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

            PeptideWithSetModifications pepA = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
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

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.That(proteinGroups.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("1"));
            Assert.That(proteinGroups.ElementAt(0).UniquePeptides.Count, Is.EqualTo(2));
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

            var protease1 = new Protease("test5", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("test6", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            PeptideWithSetModifications pepA_1Dp1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_2Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepA_3Dp2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmABC_Dp1 = new PeptideSpectralMatch(pepA_1Dp1, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            SpectralMatch psmABC_Dp2 = new PeptideSpectralMatch(pepA_2Dp2, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            psmABC_Dp2.AddOrReplace(pepA_3Dp2, 10, 0, true, new List<MatchedFragmentIon>());

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>() { psmABC_Dp1, psmABC_Dp2 };

            psms.ForEach(h => h.ResolveAllAmbiguities());
            psms.ForEach(h => h.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.That(proteinGroups.Count, Is.EqualTo(2));
            if (proteinGroups.ElementAt(0).ProteinGroupName == "1")
            {
                Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("1"));
                Assert.That(proteinGroups.ElementAt(0).UniquePeptides.Count, Is.EqualTo(1));
                Assert.That(proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).FullSequence, Is.EqualTo("ABC"));
                Assert.That(proteinGroups.ElementAt(1).ProteinGroupName, Is.EqualTo("2|3"));
                Assert.That(proteinGroups.ElementAt(1).UniquePeptides.Count, Is.EqualTo(0));
                Assert.That(proteinGroups.ElementAt(1).AllPeptides.Count, Is.EqualTo(2));
            }
            else
            {
                Assert.That(proteinGroups.ElementAt(1).ProteinGroupName, Is.EqualTo("1"));
                Assert.That(proteinGroups.ElementAt(1).UniquePeptides.Count, Is.EqualTo(1));
                Assert.That(proteinGroups.ElementAt(1).UniquePeptides.ElementAt(0).FullSequence, Is.EqualTo("ABC"));
                Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("2|3"));
                Assert.That(proteinGroups.ElementAt(0).UniquePeptides.Count, Is.EqualTo(0));
                Assert.That(proteinGroups.ElementAt(0).AllPeptides.Count, Is.EqualTo(2));
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


            var protease1 = new Protease("test3", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("test4", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            var peptideList = new List<IBioPolymerWithSetMods>();
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                    }
                }
            }

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
                psms.Last().SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                psms.Last().ResolveAllAmbiguities();
            }

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            Assert.That(proteinGroups.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("1"));
            Assert.That(proteinGroups.ElementAt(0).UniquePeptides.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).FullSequence, Is.EqualTo("ABC"));
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

            var protease1 = new Protease("testC", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("testD", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new List<IBioPolymerWithSetMods>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        case "XYZ": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(commonParameters2.DigestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "-XYZ-EFG": peptideList.Add(peptide); break;
                    }
                }
            }

            // builds psm list to match to peptides
            List<SpectralMatch> psms = new List<SpectralMatch>();
            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            double goodFdr = 0.00100;
            double badFdr = 0.0200;

            // creates psms for specific PWSM
            for (int i = 0; i < peptideList.Count(); i++)
            {
                if (peptideList[i].DigestionParams == commonParameters.DigestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(peptideList[i], 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>()));
                }
                if (peptideList[i].DigestionParams == commonParameters2.DigestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(peptideList[i], 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>()));
                }

                switch (i)
                {
                    case 0:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, badFdr, 0, 0);
                        break;

                    case 1:
                        psms.Last().SetFdrValues(0, 0, badFdr, 0, 0, goodFdr, 0, 0);
                        break;

                    case 2:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0);
                        break;

                    case 3:
                        psms.Last().SetFdrValues(0, 0, badFdr, 0, 0, badFdr, 0, 0);
                        break;

                    case 4:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0);
                        break;

                    case 5:
                        psms.Last().SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0);
                        break;
                }
                psms.Last().ResolveAllAmbiguities();
            }

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;

            psms.ElementAt(5).SetFdrValues(0, 0, goodFdr, 0, goodFdr, 0, 0, 0);

            //this iscopy of code that filteres psms in PostSearch Analysis Task
            var fdrFilteredPsms = new List<SpectralMatch>();
            foreach (SpectralMatch psm in psms)
            {
                if (psm != null && psm.FdrInfo.QValue <= 0.0100 && psm.FdrInfo.QValueNotch <= 0.0100)
                {
                    fdrFilteredPsms.Add(psm);
                }
            }

            Assert.That(fdrFilteredPsms.Count, Is.EqualTo(3));

            var test1 = fdrFilteredPsms.Contains(psms.ElementAt(2));
            var test2 = fdrFilteredPsms.Contains(psms.ElementAt(4));
            var test3 = fdrFilteredPsms.Contains(psms.ElementAt(5));
            var test4 = fdrFilteredPsms.Contains(psms.ElementAt(0));
            var test5 = fdrFilteredPsms.Contains(psms.ElementAt(1));
            var test6 = fdrFilteredPsms.Contains(psms.ElementAt(3));
            Assert.That(test1, Is.EqualTo(true));
            Assert.That(test2, Is.EqualTo(true));
            Assert.That(test3, Is.EqualTo(true));
            Assert.That(test4, Is.EqualTo(false));
            Assert.That(test5, Is.EqualTo(false));
            Assert.That(test6, Is.EqualTo(false));
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

            Assert.That(File.ReadLines(psms).Count(), Is.EqualTo(11));
            var protGroups = Path.Combine(thisTaskOutputFolder, "AllQuantifiedProteinGroups.tsv");

            Assert.That(File.ReadLines(protGroups).Count(), Is.EqualTo(7));
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

            var protease1 = new Protease("proteaseAlpha", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("proteaseBeta", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            PeptideWithSetModifications pepABC_1Alpha = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABC_2Beta = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABC_4Beta = new PeptideWithSetModifications(protein: p.ElementAt(3), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 12, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepEFG_4Beta = new PeptideWithSetModifications(protein: p.ElementAt(3), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 7, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepEFG_3Beta = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters2.DigestionParams, oneBasedStartResidueInProtein: 5, oneBasedEndResidueInProtein: 7, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "EFG", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmABC_Alpha = new PeptideSpectralMatch(pepABC_1Alpha, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            SpectralMatch psmABC_Beta = new PeptideSpectralMatch(pepABC_2Beta, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            psmABC_Beta.AddOrReplace(pepABC_4Beta, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmEFG_Beta = new PeptideSpectralMatch(pepEFG_3Beta, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            psmEFG_Beta.AddOrReplace(pepEFG_4Beta, 10, 0, true, new List<MatchedFragmentIon>());

            List<SpectralMatch> psms = new List<SpectralMatch> { psmABC_Alpha, psmABC_Beta, psmEFG_Beta };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters.DigestionParams);
            digestionParamsList.Add(commonParameters2.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.That(proteinGroups.Count, Is.EqualTo(2));
            var proteinGroup1 = proteinGroups.Where(h => h.ProteinGroupName == "1").First();
            Assert.That(proteinGroup1.UniquePeptides.Count, Is.EqualTo(1));
            Assert.That(proteinGroup1.AllPeptides.Count, Is.EqualTo(1));
            var proteinGroup2 = proteinGroups.Where(h => h.ProteinGroupName == "4").First();
            Assert.That(proteinGroup2.UniquePeptides.Count, Is.EqualTo(0));
            Assert.That(proteinGroup2.AllPeptides.Count, Is.EqualTo(2));
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

            CommonParameters commonParameters_tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));
            CommonParameters commonParameters_LysC = new CommonParameters(digestionParams: new DigestionParams(protease: "Lys-C (don't cleave before proline)", minPeptideLength: 1));

            PeptideWithSetModifications pepABCK_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_tryp.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABCK_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_tryp.DigestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABCK_1L = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_LysC.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepABCK_2L = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_LysC.DigestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_1T = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_tryp.DigestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_2T = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_tryp.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_1L = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_LysC.DigestionParams, oneBasedStartResidueInProtein: 6, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepXYZK_2L = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_LysC.DigestionParams, oneBasedStartResidueInProtein: 2, oneBasedEndResidueInProtein: 4, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "XYZK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmABCK_T = new PeptideSpectralMatch(pepABCK_1T, 0, 10, 0, scan, commonParameters_tryp, new List<MatchedFragmentIon>());
            psmABCK_T.AddOrReplace(pepABCK_2T, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmABCK_L = new PeptideSpectralMatch(pepABCK_1L, 0, 10, 0, scan, commonParameters_LysC, new List<MatchedFragmentIon>());
            psmABCK_L.AddOrReplace(pepABCK_2L, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmXYZK_T = new PeptideSpectralMatch(pepXYZK_1T, 0, 10, 0, scan, commonParameters_LysC, new List<MatchedFragmentIon>());
            psmXYZK_T.AddOrReplace(pepXYZK_2T, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmXYZK_L = new PeptideSpectralMatch(pepXYZK_1L, 0, 10, 0, scan, commonParameters_LysC, new List<MatchedFragmentIon>());
            psmXYZK_L.AddOrReplace(pepXYZK_2L, 10, 0, true, new List<MatchedFragmentIon>());

            List<SpectralMatch> psms = new List<SpectralMatch> { psmABCK_T, psmABCK_L, psmXYZK_T, psmXYZK_L };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters_tryp.DigestionParams);
            digestionParamsList.Add(commonParameters_LysC.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.That(proteinGroups.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("1|2"));
            Assert.That(proteinGroups.ElementAt(0).AllPeptides.Count, Is.EqualTo(8));
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

            var protease1 = new Protease("proteaseDash", CleavageSpecificity.Full, null, null, motifs1);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("proteaseG", CleavageSpecificity.Full, null, null, motifs2);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease1.Name, minPeptideLength: 1));
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(protease: protease2.Name, minPeptideLength: 1));

            PeptideWithSetModifications pepABC_1Dash = new PeptideWithSetModifications(p.ElementAt(0), commonParameters.DigestionParams, 2, 4, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pepABC_2Dash = new PeptideWithSetModifications(p.ElementAt(1), commonParameters.DigestionParams, 7, 9, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pepABC_2G = new PeptideWithSetModifications(p.ElementAt(1), commonParameters2.DigestionParams, 2, 4, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pepABC_3G = new PeptideWithSetModifications(p.ElementAt(2), commonParameters2.DigestionParams, 7, 9, CleavageSpecificity.Unknown, "ABCK", 0, new Dictionary<int, Modification>(), 0);

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            SpectralMatch psmABC_Dash = new PeptideSpectralMatch(pepABC_1Dash, 0, 10, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            psmABC_Dash.AddOrReplace(pepABC_2Dash, 10, 0, true, new List<MatchedFragmentIon>());
            SpectralMatch psmABC_G = new PeptideSpectralMatch(pepABC_2G, 0, 10, 0, scan, commonParameters2, new List<MatchedFragmentIon>());
            psmABC_G.AddOrReplace(pepABC_3G, 10, 0, true, new List<MatchedFragmentIon>());

            List<SpectralMatch> psms = new List<SpectralMatch> { psmABC_Dash, psmABC_G };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            var digestionParamsList = new HashSet<IDigestionParams>();
            digestionParamsList.Add(commonParameters.DigestionParams);
            digestionParamsList.Add(commonParameters2.DigestionParams);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            List<ProteinGroup> proteinGroups = results.SortedAndScoredProteinGroups;
            Assert.That(proteinGroups.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.ElementAt(0).ProteinGroupName, Is.EqualTo("2"));
            Assert.That(proteinGroups.ElementAt(0).AllPeptides.Count, Is.EqualTo(2));
        }

        /// <summary>
        /// This test ensures that FDR for each psm is calculated according to its protease
        /// </summary>
        [Test]
        public static void MultiProteaseParsimony_TestingProteaseSpecificFDRCalculations()
        {
            // two protease options
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