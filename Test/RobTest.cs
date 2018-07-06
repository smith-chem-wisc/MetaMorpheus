using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class RobTest
    {
        [Test]
        public static void TestParsimony()
        {
            // creates some test proteins and digests them (simulating a protein database)
            string[] sequences = { "AB--------",   // 1: contains unique
                                   "--C-------",   // 2: one hit wonder
                                   "---D---HHH--", // 3: subset
                                   "-B-D---HHH--", // 4: D should go to 4, not 3 (3 is subset)
                                   "-B--E-----",   // 5: subsumable
                                   "----EFG---",   // 6: indistinguishable from 8 (J will not be a "detected" PSM)
                                   "-----F----",   // 7: lone pep shared w/ decoy
                                   "--------I-",   // 8: I should go to 9, not 8
                                   "-B------I-",   // 9: I should go to 9, not 8
                                   "----EFG--J"    // 10: indistinguishable from 6 (J will not be a "detected" PSM)
                                   };

            IEnumerable<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("A", TerminusType.C), new Tuple<string, TerminusType>("B", TerminusType.C), new Tuple<string, TerminusType>("C", TerminusType.C), new Tuple<string, TerminusType>("D", TerminusType.C), new Tuple<string, TerminusType>("E", TerminusType.C), new Tuple<string, TerminusType>("F", TerminusType.C), new Tuple<string, TerminusType>("G", TerminusType.C), new Tuple<string, TerminusType>("H", TerminusType.C), new Tuple<string, TerminusType>("I", TerminusType.C), new Tuple<string, TerminusType>("J", TerminusType.C), new Tuple<string, TerminusType>("-", TerminusType.C) };
            var protease = new Protease("test", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            var peptideList = new HashSet<PeptideWithSetModifications>();
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            p.Add(new Protein("-----F----*", "D1", null, gn, new Dictionary<int, List<Modification>>(), isDecoy: true));
            p.Add(new Protein("-----F----**", "C1", null, gn, new Dictionary<int, List<Modification>>(), isContaminant: true));
            p.Add(new Protein("----E----**", "C2", null, gn, new Dictionary<int, List<Modification>>(), isContaminant: true));

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "A": peptideList.Add(peptide); break;
                        case "B": peptideList.Add(peptide); break;
                        case "C": peptideList.Add(peptide); break;
                        case "D": peptideList.Add(peptide); break;
                        case "E": peptideList.Add(peptide); break;
                        case "F": peptideList.Add(peptide); break;
                        case "G": peptideList.Add(peptide); break;
                        case "H": peptideList.Add(peptide); break;
                        case "I": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];
            HashSet<PeptideWithSetModifications>[] virtualPeptideSets = new HashSet<PeptideWithSetModifications>[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
            }

            // creates protein list
            for (int i = 0; i < virtualPeptideSets.Length; i++)
            {
                virtualPeptideSets[i] = new HashSet<PeptideWithSetModifications>();

                foreach (var virtualPeptide in peptideList)
                {
                    string peptideBaseSequence = string.Join("", peptideList.ElementAt(i).BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (virtualPeptide.BaseSequence.Contains(peptideBaseSequence))
                    {
                        virtualPeptideSets[i].Add(virtualPeptide);
                    }
                }
            }

            // populates initial peptide-virtualpeptide dictionary
            for (int i = 0; i < peptides.Length; i++)
            {
                if (!dictionary.ContainsKey(peptides[i]))
                {
                    dictionary.Add(peptides[i], virtualPeptideSets[i]);
                }
            }

            // copy for comparison later
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> initialDictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            foreach (var kvp in dictionary)
            {
                CompactPeptideBase cp = kvp.Key;
                HashSet<PeptideWithSetModifications> peps = new HashSet<PeptideWithSetModifications>();
                foreach (var pep in kvp.Value)
                    peps.Add(pep);

                initialDictionary.Add(cp, peps);
            }

            // apply parsimony to dictionary
            ProteinParsimonyEngine ae = new ProteinParsimonyEngine(dictionary, false, new CommonParameters(), new List<string>());
            var hah = (ProteinParsimonyResults)ae.Run();
            var proteinGroups = hah.ProteinGroups;

            var parsimonyProteinList = new List<Protein>();
            var parsimonyBaseSequences = new List<string>();

            foreach (var kvp in dictionary)
            {
                foreach (var virtualPeptide in kvp.Value)
                {
                    if (!parsimonyProteinList.Contains(virtualPeptide.Protein))
                    {
                        parsimonyProteinList.Add(virtualPeptide.Protein);
                        parsimonyBaseSequences.Add(virtualPeptide.Protein.BaseSequence);
                    }
                }
            }

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
                        case "A": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;
                        case "B": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 9, 0, scan, digestionParams)); break;
                        case "C": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 8, 0, scan, digestionParams)); break;
                        case "D": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 7, 0, scan, digestionParams)); break;
                        case "E": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 6, 0, scan, digestionParams)); break;
                        case "F": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 5, 0, scan, digestionParams)); break;
                        case "G": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 4, 0, scan, digestionParams)); break;
                        case "H": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 3, 0, scan, digestionParams)); break;
                        case "I": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 2, 0, scan, digestionParams)); break;
                    }
                }
            }

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            foreach (var hm in psms)
            {
                hm.MatchToProteinLinkedPeptides(initialDictionary);
                hm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            }

            ProteinScoringAndFdrEngine f = new ProteinScoringAndFdrEngine(proteinGroups, psms, true, false, true, new CommonParameters(), new List<string>());
            var ok = (ProteinScoringAndFdrResults)f.Run();
            proteinGroups = ok.SortedAndScoredProteinGroups;

            //prints initial dictionary
            List<Protein> proteinList = new List<Protein>();
            foreach (var kvp in initialDictionary)
            {
                proteinList = new List<Protein>();
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinList.Contains(peptide.Protein))
                    {
                        proteinList.Add(peptide.Protein);
                    }
                }
            }

            //prints parsimonious dictionary
            foreach (var kvp in dictionary)
            {
                proteinList = new List<Protein>();
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinList.Contains(peptide.Protein))
                    {
                        proteinList.Add(peptide.Protein);
                    }
                }
            }

            // check that correct proteins are in parsimony list
            Assert.Contains("AB--------", parsimonyBaseSequences);
            Assert.Contains("--C-------", parsimonyBaseSequences);
            Assert.Contains("-B-D---HHH--", parsimonyBaseSequences);
            Assert.Contains("-----F----*", parsimonyBaseSequences);
            Assert.Contains("----E----**", parsimonyBaseSequences);
            Assert.Contains("-B------I-", parsimonyBaseSequences);
            Assert.Contains("----EFG---", parsimonyBaseSequences);
            Assert.Contains("----EFG--J", parsimonyBaseSequences);
            Assert.AreEqual(8, parsimonyProteinList.Count);

            // sequence coverage test
            foreach (var proteinGroup in proteinGroups)
                foreach (var coverage in proteinGroup.SequenceCoveragePercent)
                    Assert.That(coverage <= 1.0);

            // protein group tests
            Assert.AreEqual(4, proteinGroups.Count);
            Assert.AreEqual(1, proteinGroups.First().Proteins.Count);
            Assert.AreEqual("AB--------", proteinGroups.First().Proteins.First().BaseSequence);
            Assert.AreEqual(4, proteinGroups.First().AllPsmsBelowOnePercentFDR.Count);
            Assert.AreEqual(19, proteinGroups.First().ProteinGroupScore);
        }

        [Test]
        public static void TestFragments()
        {
            // creates some test proteins, digest, and fragment
            string[] sequences = { "GLSDGEWQQVLNVWGK" }; // just one peptide

            var peptides = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString()));

            DigestionParams digestionParams = new DigestionParams();
            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                    peptides.Add(peptide);
            }

            var CfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var ZdotfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var BfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var YfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var BYfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();

            foreach (var peptide in peptides)
            {
                CfragmentMasses.Add(peptide, peptide.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.C }));
                ZdotfragmentMasses.Add(peptide, peptide.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.Zdot }));
                BfragmentMasses.Add(peptide, peptide.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B }));
                YfragmentMasses.Add(peptide, peptide.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.Y }));
                BYfragmentMasses.Add(peptide, peptide.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }));
            }
            Assert.That(BfragmentMasses.TryGetValue(peptides.First(), out double[] testB));

            Assert.That(YfragmentMasses.TryGetValue(peptides.First(), out double[] testY));

            Assert.That(CfragmentMasses.TryGetValue(peptides.First(), out double[] testC));

            Assert.That(ZdotfragmentMasses.TryGetValue(peptides.First(), out double[] testZ));
        }

        [Test]
        public static void TestNeutralLossFragments()
        {
            Protein p = new Protein("PR", "ac");

            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            ModificationWithMass nTermAmmoniaLoss = new ModificationWithMass("ntermammonialoss", "mt", motif, TerminusLocalization.NPep, 0, neutralLosses: new List<double> { 0, -17 });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);

            var cool = p.Digest(digestionParams, new List<ModificationWithMass> { nTermAmmoniaLoss }, new List<ModificationWithMass>()).First();
            var nice = cool.CompactPeptide(TerminusType.None);
            Assert.AreEqual(2, nice.NTerminalMasses.Length);
            Assert.AreEqual(1, nice.CTerminalMasses.Length);
        }

        [Test]
        public static void TestPTMOutput()
        {
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

            ModificationMotif.TryGetMotif("S", out ModificationMotif motif);
            variableModifications.Add(new ModificationWithMassAndCf("resMod", "HaHa", motif, TerminusLocalization.Any, ChemicalFormula.ParseFormula("H")));

            var proteinList = new List<Protein> { new Protein("MNNNSKQQQ", "accession") };
            var protease = new Protease("CustomProtease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>
            {
                {variableModifications.Last(), 1 }
            };

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1);

            var modPep = proteinList.First().Digest(digestionParams, fixedModifications, variableModifications).Last();
            HashSet<PeptideWithSetModifications> value = new HashSet<PeptideWithSetModifications> { modPep };
            CompactPeptide compactPeptide1 = new CompactPeptide(value.First(), TerminusType.None);
            Assert.AreEqual("QQQ", value.First().Sequence);

            var firstProtDigest = proteinList.First().Digest(digestionParams, fixedModifications, variableModifications).ToList();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { firstProtDigest[0] };
            CompactPeptide compactPeptide2 = new CompactPeptide(value2.First(), TerminusType.None);
            Assert.AreEqual("MNNNSK", value2.First().Sequence);

            HashSet<PeptideWithSetModifications> value2mod = new HashSet<PeptideWithSetModifications> { firstProtDigest[1] };
            CompactPeptide compactPeptide2mod = new CompactPeptide(value2mod.Last(), TerminusType.None);
            Assert.AreEqual("MNNNS[HaHa:resMod]K", value2mod.Last().Sequence);

            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { firstProtDigest[2] };
            CompactPeptide compactPeptide3 = new CompactPeptide(value3.First(), TerminusType.None);
            Assert.AreEqual("NNNSK", value3.First().Sequence);
            HashSet<PeptideWithSetModifications> value3mod = new HashSet<PeptideWithSetModifications> { firstProtDigest[3] };

            CompactPeptide compactPeptide3mod = new CompactPeptide(value3mod.Last(), TerminusType.None);
            Assert.AreEqual("NNNS[HaHa:resMod]K", value3mod.Last().Sequence);

            var peptideList = new HashSet<PeptideWithSetModifications>();
            foreach (var protein in proteinList)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications))
                {
                    peptideList.Add(peptide);
                }
            }

            compactPeptideToProteinPeptideMatching.Add(compactPeptide1, value);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide2, value2);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide3, value3);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide2mod, value2mod);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide3mod, value3mod);

            ProteinParsimonyEngine engine = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, true, new CommonParameters(), new List<string> { "ff" });
            var cool = (ProteinParsimonyResults)engine.Run();
            var proteinGroups = cool.ProteinGroups;

            MsDataScan jdfk = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass ms2scan = new Ms2ScanWithSpecificMass(jdfk, 2, 0, "File");

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            var match1 = new PeptideSpectralMatch(peptideList.ElementAt(0).CompactPeptide(TerminusType.None), 0, 10, 0, ms2scan, digestionParams)
            {
            };
            match1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            var match2 = new PeptideSpectralMatch(peptideList.ElementAt(1).CompactPeptide(TerminusType.None), 0, 10, 0, ms2scan, digestionParams)
            {
            };
            match2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            var match3 = new PeptideSpectralMatch(peptideList.ElementAt(1).CompactPeptide(TerminusType.None), 0, 10, 0, ms2scan, digestionParams)
            {
            };
            match3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            match1.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            match2.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            match3.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>
            {
                match1,
                match2,
                match3
            };
            ProteinScoringAndFdrEngine f = new ProteinScoringAndFdrEngine(proteinGroups, psms, false, false, true, new CommonParameters(), new List<string>());
            f.Run();

            Assert.AreEqual("#aa5[resMod,info:occupancy=0.67(2/3)];", proteinGroups.First().ModsInfo[0]);
        }

        [Test]
        public static void TestProteinGroupsAccessionOutputOrder()
        {
            var p = new HashSet<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();

            // make protein B
            p.Add(new Protein("-----F----*", "B", null, gn, new Dictionary<int, List<Modification>>(), isDecoy: true));

            // make protein A
            p.Add(new Protein("-----F----**", "A", null, gn, new Dictionary<int, List<Modification>>(), isDecoy: true));

            // add protein B and A to the protein group
            ProteinGroup testGroup = new ProteinGroup(p, null, null);

            // test order is AB and not BA
            Assert.That(testGroup.ProteinGroupName.Equals("A|B"));
            Assert.That(testGroup.Proteins.First().Accession.Equals("B"));
        }
    }
}