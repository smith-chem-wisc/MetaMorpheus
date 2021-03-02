using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Test
{
    [TestFixture]
    public static class RescueAndResolveTest
    {
        [Test]
        public static void rescueCase1_true()
        {

            string[] sequences = {
                "PEPTIDEKPEPTIDESRPEPTIDESK",
                "PEPTIDESRPEPTIDESK"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 30.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 30.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_unique = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 18, oneBasedEndResidueInProtein: 26, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot2_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmShared2 = new PeptideSpectralMatch(prot1_shared2, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared2.AddOrReplace(prot2_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmUnique1 = new PeptideSpectralMatch(prot1_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmShared2, psmUnique1 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(2, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(3, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
            Assert.AreEqual("2", results.SortedAndScoredProteinGroups[1].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[1].AllPeptides.Count);
        }

        [Test]
        public static void rescueCase1_false()
        {

            string[] sequences = {
                "PEPTIDEKPEPTIDESRPEPTIDESK",
                "PEPTIDESRPEPTIDESK"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 5.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 5.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_unique = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 18, oneBasedEndResidueInProtein: 26, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot2_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmShared2 = new PeptideSpectralMatch(prot1_shared2, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared2.AddOrReplace(prot2_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmUnique1 = new PeptideSpectralMatch(prot1_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmShared2, psmUnique1 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(1, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(3, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
        }

        [Test]
        public static void rescueCase1_2_true()
        {
            string[] sequences = {
                "PEPTIDEKPEPTIDESR",
                "PEPTIDESR",
                "PEPTIDESRPEPTIDESK"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 25.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 25.0);
            LongReadInfo prot3 = new LongReadInfo(p.ElementAt(2).Accession, 25.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);
            proteinToProteogenomicInfo.Add(p.ElementAt(2).Accession, prot3);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_unique = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_unique = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot2_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            psmShared1.AddOrReplace(prot3_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmUnique3 = new PeptideSpectralMatch(prot3_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmUnique1 = new PeptideSpectralMatch(prot1_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmUnique3, psmUnique1 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(3, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
            Assert.AreEqual("2", results.SortedAndScoredProteinGroups[1].ProteinGroupName);
            Assert.AreEqual(1, results.SortedAndScoredProteinGroups[1].AllPeptides.Count);
            Assert.AreEqual("3", results.SortedAndScoredProteinGroups[2].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[2].AllPeptides.Count);
        }

        [Test]
        public static void rescueCase1_2_false()
        {
            string[] sequences = {
                "PEPTIDEKPEPTIDESR",
                "PEPTIDESR",
                "PEPTIDESRPEPTIDESK"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 25.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 5.0);
            LongReadInfo prot3 = new LongReadInfo(p.ElementAt(2).Accession, 25.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);
            proteinToProteogenomicInfo.Add(p.ElementAt(2).Accession, prot3);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_unique = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_unique = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot2_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            psmShared1.AddOrReplace(prot3_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmUnique3 = new PeptideSpectralMatch(prot3_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmUnique1 = new PeptideSpectralMatch(prot1_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmUnique3, psmUnique1 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(2, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
            Assert.AreEqual("3", results.SortedAndScoredProteinGroups[1].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[1].AllPeptides.Count);

        }

        [Test]
        public static void rescueCase2_true()
        {
            string[] sequences = {
                "PEPTIDEKPEPTIDESR",
                "PEPTIDESRPAPTIDESR",
                "PAPTIDESRPEPTADESK"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 25.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 25.0);
            LongReadInfo prot3 = new LongReadInfo(p.ElementAt(2).Accession, 25.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);
            proteinToProteogenomicInfo.Add(p.ElementAt(2).Accession, prot3);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_unique = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_unique = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot2_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmShared2 = new PeptideSpectralMatch(prot2_shared2, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared2.AddOrReplace(prot3_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmUnique3 = new PeptideSpectralMatch(prot3_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmUnique1 = new PeptideSpectralMatch(prot1_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmUnique3, psmUnique1, psmShared2 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(3, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
            Assert.AreEqual("2", results.SortedAndScoredProteinGroups[1].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[1].AllPeptides.Count);
            Assert.AreEqual("3", results.SortedAndScoredProteinGroups[2].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[2].AllPeptides.Count);
        }

        [Test]
        public static void rescueCase2_false()
        {
            string[] sequences = {
                "PEPTIDEKPEPTIDESR",
                "PEPTIDESRPAPTIDESR",
                "PAPTIDESRPEPTADESK"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 25.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 10.0);
            LongReadInfo prot3 = new LongReadInfo(p.ElementAt(2).Accession, 25.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);
            proteinToProteogenomicInfo.Add(p.ElementAt(2).Accession, prot3);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_unique = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_unique = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot2_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmShared2 = new PeptideSpectralMatch(prot2_shared2, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared2.AddOrReplace(prot3_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            PeptideSpectralMatch psmUnique3 = new PeptideSpectralMatch(prot3_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            PeptideSpectralMatch psmUnique1 = new PeptideSpectralMatch(prot1_unique, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmUnique3, psmUnique1, psmShared2 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(2, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
            Assert.AreEqual("3", results.SortedAndScoredProteinGroups[1].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[1].AllPeptides.Count);
        }

        [Test]
        public static void rescueCase3_true()
        {
            string[] sequences = {
                "PEPTIDEKPEPTIDESRPEPTIDESK",
                "PEPTIDESRPEPTIDESK",
                "PEPTIDEKPEPTIDESR"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 25.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 25.0);
            LongReadInfo prot3 = new LongReadInfo(p.ElementAt(2).Accession, 25.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);
            proteinToProteogenomicInfo.Add(p.ElementAt(2).Accession, prot3);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared3 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 18, oneBasedEndResidueInProtein: 26, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared3 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared2 = new PeptideSpectralMatch(prot1_shared2, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared2.AddOrReplace(prot2_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            psmShared2.AddOrReplace(prot3_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot3_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);

            PeptideSpectralMatch psmShared3 = new PeptideSpectralMatch(prot1_shared3, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared3.AddOrReplace(prot2_shared3, 10, 0, true, new List<MatchedFragmentIon>(), 0);


            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmShared3, psmShared2 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(3, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(3, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
            Assert.AreEqual("3", results.SortedAndScoredProteinGroups[1].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[1].AllPeptides.Count);
            Assert.AreEqual("2", results.SortedAndScoredProteinGroups[2].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[2].AllPeptides.Count);
        }

        [Test]
        public static void rescueCase3_true_remove()
        {
            string[] sequences = {
                "PEPTIDEKPEPTIDESRPEPTIDESK",
                "PEPTIDESRPEPTIDESK",
                "PEPTIDEKPEPTIDESR"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 5.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 25.0);
            LongReadInfo prot3 = new LongReadInfo(p.ElementAt(2).Accession, 25.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);
            proteinToProteogenomicInfo.Add(p.ElementAt(2).Accession, prot3);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared3 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 18, oneBasedEndResidueInProtein: 26, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared3 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared2 = new PeptideSpectralMatch(prot1_shared2, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared2.AddOrReplace(prot2_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            psmShared2.AddOrReplace(prot3_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot3_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);

            PeptideSpectralMatch psmShared3 = new PeptideSpectralMatch(prot1_shared3, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared3.AddOrReplace(prot2_shared3, 10, 0, true, new List<MatchedFragmentIon>(), 0);


            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmShared3, psmShared2 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(2, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("3", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
            Assert.AreEqual("2", results.SortedAndScoredProteinGroups[1].ProteinGroupName);
            Assert.AreEqual(2, results.SortedAndScoredProteinGroups[1].AllPeptides.Count);
        }

        [Test]
        public static void rescueCase3_false()
        {
            string[] sequences = {
                "PEPTIDEKPEPTIDESRPEPTIDESK",
                "PEPTIDESRPEPTIDESK",
                "PEPTIDEKPEPTIDESR"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }
            LongReadInfo prot1 = new LongReadInfo(p.ElementAt(0).Accession, 25.0);
            LongReadInfo prot2 = new LongReadInfo(p.ElementAt(1).Accession, 5.0);
            LongReadInfo prot3 = new LongReadInfo(p.ElementAt(2).Accession, 5.0);

            Dictionary<string, LongReadInfo> proteinToProteogenomicInfo = new Dictionary<string, LongReadInfo>();
            proteinToProteogenomicInfo.Add(p.ElementAt(0).Accession, prot1);
            proteinToProteogenomicInfo.Add(p.ElementAt(1).Accession, prot2);
            proteinToProteogenomicInfo.Add(p.ElementAt(2).Accession, prot3);

            GlobalVariables.ProteinToProteogenomicInfo = proteinToProteogenomicInfo;

            CommonParameters commonParameters_Tryp = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications prot1_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot1_shared3 = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 18, oneBasedEndResidueInProtein: 26, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 9, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot2_shared3 = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 10, oneBasedEndResidueInProtein: 18, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared2 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 17, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDESR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications prot3_shared1 = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParameters_Tryp.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDEK", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);


            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            PeptideSpectralMatch psmShared2 = new PeptideSpectralMatch(prot1_shared2, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared2.AddOrReplace(prot2_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);
            psmShared2.AddOrReplace(prot3_shared2, 10, 0, true, new List<MatchedFragmentIon>(), 0);

            PeptideSpectralMatch psmShared1 = new PeptideSpectralMatch(prot1_shared1, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared1.AddOrReplace(prot3_shared1, 10, 0, true, new List<MatchedFragmentIon>(), 0);

            PeptideSpectralMatch psmShared3 = new PeptideSpectralMatch(prot1_shared3, 0, 10, 0, scan, commonParameters_Tryp, new List<MatchedFragmentIon>());
            psmShared3.AddOrReplace(prot2_shared3, 10, 0, true, new List<MatchedFragmentIon>(), 0);


            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmShared1, psmShared3, psmShared2 };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, true, 25.0, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            // score protein groups and merge indistinguishable ones
            ProteinScoringAndFdrEngine proteinScoringEngine = new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms, false, true, true, new CommonParameters(), null, new List<string>());
            var results = (ProteinScoringAndFdrResults)proteinScoringEngine.Run();

            Assert.AreEqual(1, results.SortedAndScoredProteinGroups.Count);
            Assert.AreEqual("1", results.SortedAndScoredProteinGroups[0].ProteinGroupName);
            Assert.AreEqual(3, results.SortedAndScoredProteinGroups[0].AllPeptides.Count);
        }

    }
}
