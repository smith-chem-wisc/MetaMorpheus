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
        public static void MultiProteaseUniqueTest()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
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
    }
}