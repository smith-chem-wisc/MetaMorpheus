using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

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

            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("-", FragmentationTerminus.C), new Tuple<string, FragmentationTerminus>("Z", FragmentationTerminus.C) };
            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage2 = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C) };

            var protease1 = new Protease("test1", sequencesInducingCleavage, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("test2", sequencesInducingCleavage2, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
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

            foreach (var protein in proteinList)
            {
                //using protease #1
                foreach (var pwsm in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (pwsm.BaseSequence)
                    {
                        case "ABC": pwsmList.Add(pwsm); break;
                        case "EFGABC": pwsmList.Add(pwsm); break;
                        case "XYZ": pwsmList.Add(pwsm); break;
                    }
                }
                //using protease #2
                foreach (var pwsm in protein.Digest(digestionParams2, new List<Modification>(), new List<Modification>()))
                {
                    switch (pwsm.BaseSequence)
                    {
                        case "ABC": pwsmList.Add(pwsm); break;
                        case "-XYZ-EFG": pwsmList.Add(pwsm); break;
                    }
                }
            }

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            // dummy scan
            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var peptide in pwsmList)
            {
                switch (peptide.BaseSequence)
                {
                    case "ABC":
                        if (peptide.DigestionParams == digestionParams)
                        {
                            psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                            break;
                        }
                        if (peptide.DigestionParams == digestionParams2)
                        {
                            psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>()));
                            break;
                        }
                        else { break; }

                    case "EFGABC": psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;

                    case "XYZ": psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;

                    case "-XYZ-EFG": psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>())); break;
                }
            }

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
            Assert.That(pg1pep1.DigestionParams.Protease.Name == "test1");
            var pg1pep2 = proteinGroup1.AllPeptides.Where(p => p.BaseSequence == "ABC").First();
            Assert.That(pg1pep2.DigestionParams.Protease.Name == "test1");
            Assert.That(proteinGroup1.UniquePeptides.First().BaseSequence.Equals("ABC"));

            var proteinGroup2 = proteinGroups.Where(p => p.ProteinGroupName == "2").First();
            Assert.AreEqual(4, proteinGroup2.AllPeptides.Count);
            Assert.AreEqual(3, proteinGroup2.UniquePeptides.Count);
            var pg2pep1 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "XYZ").First();
            Assert.That(pg2pep1.DigestionParams.Protease.Name == "test1");
            var pg2pep2 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "ABC").First();
            Assert.That(pg2pep2.DigestionParams.Protease.Name == "test2");
            var pg2pep3 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "EFGABC").First();
            Assert.That(pg2pep3.DigestionParams.Protease.Name == "test1");
            var pg2pep4 = proteinGroup2.AllPeptides.Where(p => p.BaseSequence == "-XYZ-EFG").First();
            Assert.That(pg2pep4.DigestionParams.Protease.Name == "test2");
            var uniquePeptideSequences = proteinGroup2.UniquePeptides.Select(p => p.BaseSequence).ToList();
            Assert.That(uniquePeptideSequences.Contains("ABC"));
            Assert.That(uniquePeptideSequences.Contains("EFGABC"));
            Assert.That(uniquePeptideSequences.Contains("-XYZ-EFG"));
        }

        /// <summary>
        /// These protein groups would normally be indistinguishable but not with multiprotease
        /// We expect 2 protein groups out at the end!
        /// </summary>
        [Test]
        public static void MultiProteaseIndistiguishableTest()
        {
            string[] sequences = {
                "ABCEFG",
                "EFGABC",
            };

            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("C", FragmentationTerminus.C) };
            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage2 = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C) };

            var protease = new Protease("testA", sequencesInducingCleavage, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("testB", sequencesInducingCleavage2, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
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

            //generates HashSet of PeptidesWithSetMods
            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
            }

            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            // creates psms for specific PWSM
            foreach (var PWSM in peptideList)
            {
                if (PWSM.DigestionParams == digestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(PWSM, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                }
                if (PWSM.DigestionParams == digestionParams2)
                {
                    psms.Add(new PeptideSpectralMatch(PWSM, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>()));
                }
                psms.Last().SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
                psms.Last().ResolveAllAmbiguities();
            }

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);

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