﻿using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Text;
using EngineLayer;
using MassSpectrometry;
using Proteomics;
using System.Linq;
using IO.MzML;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using System.IO;

namespace Test
{
    [TestFixture]
    public static class MultiProteaseParsimonyTest
    {

        [Test]
        public static void MultiProteaseOneProteinProduceSameUniquePeptideFor2Proteases()
        {
            string[] sequences = {
                "-XYZ-EFGABC"                
            };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("A", TerminusType.N), new Tuple<string, TerminusType>("Z", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("test1", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test2", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
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

            PeptideWithSetModifications pepA = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams:digestionParams, oneBasedStartResidueInProtein:9, oneBasedEndResidueInProtein:11, peptideDescription:"ABC", missedCleavages:0,allModsOneIsNterminus:new Dictionary<int, ModificationWithMass>(),numFixedMods:0);
            PeptideWithSetModifications pepB = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: digestionParams2, oneBasedStartResidueInProtein: 9, oneBasedEndResidueInProtein: 11, peptideDescription: "ABC", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, ModificationWithMass>(), numFixedMods: 0);
            peptideList.Add(pepA);

            peptideList.Add(pepB);
            
            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, List<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new List<PeptideWithSetModifications> { pepA, pepB });
           
            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            
                            else { break; }

                       
                    }
                }
            }
            foreach (var psm in psms)
            {
                psm.SetFdrValues(0, 0, 0.00, 0, 0, 0.00, 0, 0, 0, false);
            }
            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Oxidation of M", "Common Variable", motif2, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, true, new CommonParameters(), null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;
            
            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;

            foreach (var PSM in psms)
            {
                PSM.MatchToProteinLinkedPeptides(CompactPeptidesToProteinPeptidesMatching);
            }

            var PSMS = psms;
            var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms,
                    false, false, true, new CommonParameters(), new List<string> { "Search" }).Run();
            var tests = psms;
            var ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            Assert.AreEqual(1, ProteinGroups.Count);
            Assert.AreEqual("1", ProteinGroups.ElementAt(0).ProteinGroupName);
            Assert.AreEqual(1, ProteinGroups.ElementAt(0).UniquePeptides.Count);

        }

        [Test]
        public static void MultiProteaseUniqueTest()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
                "-XYZ-GABC"
            };

            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("-", TerminusType.C), new Tuple<string, TerminusType>("Z", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("test7", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test8", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
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

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        case "XYZ": peptideList.Add(peptide); break;
                        case "GABC": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "-XYZ-EFG": peptideList.Add(peptide); break;
                        case "-XYZ-G": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[2], PWSM[6] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1], PWSM[5], PWSM[9] });
            dictionary.Add(peptides[3], new HashSet<PeptideWithSetModifications> { PWSM[3] });
            dictionary.Add(peptides[4], new HashSet<PeptideWithSetModifications> { PWSM[4] });
            dictionary.Add(peptides[7], new HashSet<PeptideWithSetModifications> { PWSM[7] });
            dictionary.Add(peptides[8], new HashSet<PeptideWithSetModifications> { PWSM[8] });

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }

                        case "EFGABC": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;

                        case "GABC": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;

                        case "XYZ": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;

                        case "-XYZ-EFG": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2)); break;

                        case "-XYZ-G": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2)); break;
                    }
                }
            }
            foreach (var psm in psms)
            {
                psm.SetFdrValues(0, 0, 0.00, 0, 0, 0.00, 0, 0, 0, false);
            }
            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Oxidation of M", "Common Variable", motif2, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, true, new CommonParameters(), null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;
            Assert.AreEqual(6, CompactPeptidesToProteinPeptidesMatching.Count);

            Assert.AreEqual(3, CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.Count);
            Assert.AreEqual("XYZ", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("XYZ", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(1).BaseSequence);
            Assert.AreEqual("XYZ", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(2).BaseSequence);

            Assert.AreEqual(3, CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.Count);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(1).BaseSequence);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(2).BaseSequence);

            Assert.AreEqual(1, CompactPeptidesToProteinPeptidesMatching.ElementAt(2).Value.Count);
            Assert.AreEqual("EFGABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(2).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("2", CompactPeptidesToProteinPeptidesMatching.ElementAt(2).Value.ElementAt(0).Protein.Accession);

            Assert.AreEqual(1, CompactPeptidesToProteinPeptidesMatching.ElementAt(3).Value.Count);
            Assert.AreEqual("-XYZ-EFG", CompactPeptidesToProteinPeptidesMatching.ElementAt(3).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("2", CompactPeptidesToProteinPeptidesMatching.ElementAt(3).Value.ElementAt(0).Protein.Accession);

            Assert.AreEqual(1, CompactPeptidesToProteinPeptidesMatching.ElementAt(4).Value.Count);
            Assert.AreEqual("GABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(4).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("3", CompactPeptidesToProteinPeptidesMatching.ElementAt(4).Value.ElementAt(0).Protein.Accession);

            Assert.AreEqual(1, CompactPeptidesToProteinPeptidesMatching.ElementAt(5).Value.Count);
            Assert.AreEqual("-XYZ-G", CompactPeptidesToProteinPeptidesMatching.ElementAt(5).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("3", CompactPeptidesToProteinPeptidesMatching.ElementAt(5).Value.ElementAt(0).Protein.Accession);

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;

            foreach (var PSM in psms)
            {
                PSM.MatchToProteinLinkedPeptides(CompactPeptidesToProteinPeptidesMatching);
            }

            var PSMS = psms;
            var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms,
                    false, false, true, new CommonParameters(), new List<string> { "Search" }).Run();
            var tests = psms;
            var ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;
            Assert.AreEqual(3, ProteinGroups.Count);
            foreach (var pg in ProteinGroups)
            {
                if (pg.ProteinGroupName == "1")
                {
                    Assert.AreEqual(2, pg.AllPeptides.Count);
                    Assert.AreEqual(1, pg.UniquePeptides.Count);
                }
                if (pg.ProteinGroupName == "2")
                {
                    Assert.AreEqual(4, pg.AllPeptides.Count);
                    Assert.AreEqual(2, pg.UniquePeptides.Count);
                }
                if (pg.ProteinGroupName == "3")
                {
                    Assert.AreEqual(4, pg.AllPeptides.Count);
                    Assert.AreEqual(2, pg.UniquePeptides.Count);
                }
            }
        }
        [Test]
        public static void MultiProteaseUniqueAndSharedBothPresent()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
                "-XYZ-GABC"
            };

            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("-", TerminusType.C), new Tuple<string, TerminusType>("Z", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("test5", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test6", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
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

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        case "XYZ": peptideList.Add(peptide); break;
                        case "GABC": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "-XYZ-EFG": peptideList.Add(peptide); break;
                        case "-XYZ-G": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[2], PWSM[6] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1], PWSM[5], PWSM[9] });
            dictionary.Add(peptides[3], new HashSet<PeptideWithSetModifications> { PWSM[3] });
            dictionary.Add(peptides[4], new HashSet<PeptideWithSetModifications> { PWSM[4] });
            dictionary.Add(peptides[7], new HashSet<PeptideWithSetModifications> { PWSM[7] });
            dictionary.Add(peptides[8], new HashSet<PeptideWithSetModifications> { PWSM[8] });

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }

                    }
                }
            }
            foreach (var psm in psms)
            {
                psm.SetFdrValues(0, 0, 0.00, 0, 0, 0.00, 0, 0, 0, false);
            }
            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Oxidation of M", "Common Variable", motif2, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, true, new CommonParameters(), null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;
            

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;

            foreach (var PSM in psms)
            {
                PSM.MatchToProteinLinkedPeptides(CompactPeptidesToProteinPeptidesMatching);
            }

            var PSMS = psms;
            var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms,
                    false, false, true, new CommonParameters(), new List<string> { "Search" }).Run();
            var tests = psms;
            var ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;
            Assert.AreEqual(2, ProteinGroups.Count);
            if (ProteinGroups.ElementAt(0).ProteinGroupName == "1")
            {
                Assert.AreEqual("1", ProteinGroups.ElementAt(0).ProteinGroupName);
                Assert.AreEqual(1, ProteinGroups.ElementAt(0).UniquePeptides.Count);
                Assert.AreEqual("ABC", ProteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).Sequence);
                Assert.AreEqual("2|3", ProteinGroups.ElementAt(1).ProteinGroupName);
                Assert.AreEqual(0, ProteinGroups.ElementAt(1).UniquePeptides.Count);
            }
            else
            {
                Assert.AreEqual("1", ProteinGroups.ElementAt(1).ProteinGroupName);
                Assert.AreEqual(1, ProteinGroups.ElementAt(1).UniquePeptides.Count);
                Assert.AreEqual("ABC", ProteinGroups.ElementAt(1).UniquePeptides.ElementAt(0).Sequence);
                Assert.AreEqual("2|3", ProteinGroups.ElementAt(0).ProteinGroupName);
                Assert.AreEqual(0, ProteinGroups.ElementAt(0).UniquePeptides.Count);
            }
           
        }

        [Test]
        public static void MultiProteaseUniqueAndSharedUniquePresent()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
                "-XYZ-GABC"
            };

            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("-", TerminusType.C), new Tuple<string, TerminusType>("Z", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("test3", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test4", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
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

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        case "XYZ": peptideList.Add(peptide); break;
                        case "GABC": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "-XYZ-EFG": peptideList.Add(peptide); break;
                        case "-XYZ-G": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[2], PWSM[6] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1], PWSM[5], PWSM[9] });
            dictionary.Add(peptides[3], new HashSet<PeptideWithSetModifications> { PWSM[3] });
            dictionary.Add(peptides[4], new HashSet<PeptideWithSetModifications> { PWSM[4] });
            dictionary.Add(peptides[7], new HashSet<PeptideWithSetModifications> { PWSM[7] });
            dictionary.Add(peptides[8], new HashSet<PeptideWithSetModifications> { PWSM[8] });

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            
                            else { break; }

                    }
                }
            }
            foreach (var psm in psms)
            {
                psm.SetFdrValues(0, 0, 0.00, 0, 0, 0.00, 0, 0, 0, false);
            }
            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Oxidation of M", "Common Variable", motif2, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, true, new CommonParameters(), null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;
            
            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;

            foreach (var PSM in psms)
            {
                PSM.MatchToProteinLinkedPeptides(CompactPeptidesToProteinPeptidesMatching);
            }
                        
            var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psms,
                    false, false, true, new CommonParameters(), new List<string> { "Search" }).Run();            
            var ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;
            Assert.AreEqual(1, ProteinGroups.Count);
            Assert.AreEqual("1", ProteinGroups.ElementAt(0).ProteinGroupName);
            Assert.AreEqual(1, ProteinGroups.ElementAt(0).UniquePeptides.Count);
            Assert.AreEqual("ABC", ProteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).Sequence);
        }

        [Test]
        public static void MultiProteaseIndistiguishableTest()
        {
            string[] sequences = {
                "ABCEFG",
                "EFGABC",
            };

            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("C", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("testA", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("testB", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
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

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {

                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[3] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1], PWSM[2] });

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }

                        case "EFG":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }
                    }
                }
            }

            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, true, new CommonParameters(), null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;
            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.Count);

            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.Count);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(1).BaseSequence);

            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.Count);
            Assert.AreEqual("EFG", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("EFG", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(1).BaseSequence);

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;
            Assert.AreEqual(2, proteinGroups.Count);

            Assert.AreEqual(2, proteinGroups.ElementAt(0).AllPeptides.Count);
            Assert.AreEqual(2, proteinGroups.ElementAt(0).UniquePeptides.Count);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("testA", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).BaseSequence);
            Assert.AreEqual("testA", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(1).BaseSequence);

            Assert.AreEqual(2, proteinGroups.ElementAt(1).AllPeptides.Count);
            Assert.AreEqual(2, proteinGroups.ElementAt(1).UniquePeptides.Count);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("testB", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).BaseSequence);
            Assert.AreEqual("testB", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(1).BaseSequence);
        }
    
        [Test]
        public static void FdrFilteringPsmsTest()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
            };

            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("-", TerminusType.C), new Tuple<string, TerminusType>("Z", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("testC", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("testD", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
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

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        case "XYZ": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "-XYZ-EFG": peptideList.Add(peptide); break;
                    }
                }
            }
            
            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[2] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1], PWSM[5] });
            dictionary.Add(peptides[3], new HashSet<PeptideWithSetModifications> { PWSM[3] });
            dictionary.Add(peptides[4], new HashSet<PeptideWithSetModifications> { PWSM[4] });

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }

                        case "EFGABC": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;

                        case "XYZ": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;

                        case "-XYZ-EFG": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2)); break;
                    }
                }
            }
            double goodFdr = 0.00100;
            double badFdr = 0.0200;
                        
            psms.ElementAt(0).SetFdrValues(0, 0, goodFdr, 0, 0, badFdr, 0, 0, 0, false);
            psms.ElementAt(1).SetFdrValues(0, 0, badFdr, 0, 0, goodFdr, 0, 0, 0, false);
            psms.ElementAt(2).SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
            psms.ElementAt(3).SetFdrValues(0, 0, badFdr, 0, 0, badFdr, 0, 0, 0, false);
            psms.ElementAt(4).SetFdrValues(0, 0, goodFdr, 0, 0, goodFdr, 0, 0, 0, false);
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

        [Test]
        public static void FdrFilteredParsimonyTest()
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
            var results = Task1.RunTask(Environment.CurrentDirectory, new List<DbForTask>() { new DbForTask(fastaName, false) }, new List<string>() { mzmlName }, "test");

            var thisTaskOutputFolder = MySetUpClass.outputFolder;
                        
            var psms = Path.Combine(thisTaskOutputFolder,  "AllPSMs.psmtsv");
            
            Assert.AreEqual(12, File.ReadLines(psms).Count());
            var protGroups = Path.Combine(thisTaskOutputFolder, "AllProteinGroups.tsv");

            Assert.AreEqual(7, File.ReadLines(protGroups).Count());
        }
    }
}