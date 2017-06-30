using Chemistry;
using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using FlashLFQ;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class RobTest
    {

        #region Public Methods

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

            IEnumerable<string> sequencesInducingCleavage = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "-" };
            var protease = new Protease("test", sequencesInducingCleavage, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString(), gn, new Dictionary<int, List<Modification>>()));
            p.Add(new Protein("-----F----*", "D1", gn, new Dictionary<int, List<Modification>>(), isDecoy: true));
            p.Add(new Protein("-----F----**", "C1", gn, new Dictionary<int, List<Modification>>(), isContaminant: true));
            p.Add(new Protein("----E----**", "C2", gn, new Dictionary<int, List<Modification>>(), isContaminant: true));

            IEnumerable<PeptideWithPossibleModifications> temp;
            IEnumerable<PeptideWithSetModifications> pepWithSetMods = null;
            foreach (var protein in p)
            {
                temp = protein.Digest(protease, 2, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>());

                foreach (var dbPeptide in temp)
                {
                    pepWithSetMods = dbPeptide.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4098, 3);
                    foreach (var peptide in pepWithSetMods)
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
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];
            HashSet<PeptideWithSetModifications>[] virtualPeptideSets = new HashSet<PeptideWithSetModifications>[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), modsDictionary);
            }

            // creates protein list
            for (int i = 0; i < virtualPeptideSets.Length; i++)
            {
                virtualPeptideSets[i] = new HashSet<PeptideWithSetModifications>();

                foreach (var virtualPeptide in peptideList)
                {
                    string peptideBaseSequence = string.Join("", peptides[i].BaseSequence.Select(b => char.ConvertFromUtf32(b)));

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
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> initialDictionary = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
            foreach (var kvp in dictionary)
            {
                CompactPeptide cp = kvp.Key;
                HashSet<PeptideWithSetModifications> peps = new HashSet<PeptideWithSetModifications>();
                foreach (var pep in kvp.Value)
                    peps.Add(pep);

                initialDictionary.Add(cp, peps);
            }

            // apply parsimony to dictionary
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            AnalysisEngine ae = new AnalysisEngine(new PsmParent[0][], dictionary, new List<Protein>(), null, null, null, null, null, null, null, null, null, null, true, true, true, 0, null, null, 0, false, new List<ProductType> { ProductType.B, ProductType.Y }, double.NaN, InitiatorMethionineBehavior.Variable, new List<string>(), modsDictionary, null, null);
            ae.ApplyProteinParsimony(out proteinGroups);

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
            List<NewPsmWithFdr> psms = new List<NewPsmWithFdr>();

            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> dfb = new MzmlScanWithPrecursor(0, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, new MzPeak(2, 2), 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "A": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 10, 0, scan))); break;
                        case "B": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 9, 0, scan))); break;
                        case "C": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 8, 0, scan))); break;
                        case "D": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 7, 0, scan))); break;
                        case "E": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 6, 0, scan))); break;
                        case "F": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 5, 0, scan))); break;
                        case "G": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 4, 0, scan))); break;
                        case "H": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 3, 0, scan))); break;
                        case "I": psms.Add(new NewPsmWithFdr(new PsmClassic(peptide, 0, 2, 0, scan))); break;
                    }
                }
            }

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);

            foreach (var hm in psms)
            {
                hm.thisPSM.GetProteinLinkedInfo(initialDictionary, fragmentTolerance, lp, modsDictionary);
            }

            //Console.WriteLine(psms.Count);
            //foreach (var ok in psms)
            //{
            //    Console.WriteLine(ok);
            //}

            ae.ScoreProteinGroups(proteinGroups, psms);
            proteinGroups = ae.DoProteinFdr(proteinGroups);

            //prints initial dictionary
            List<Protein> proteinList = new List<Protein>();
            System.Console.WriteLine("----Initial Dictionary----");
            System.Console.WriteLine("PEPTIDE\t\t\tPROTEIN");
            foreach (var kvp in initialDictionary)
            {
                proteinList = new List<Protein>();
                System.Console.Write(string.Join("", kvp.Key.BaseSequence.Select(b => char.ConvertFromUtf32(b))) + "  \t\t\t  ");
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinList.Contains(peptide.Protein))
                    {
                        Console.Write(peptide.Protein.BaseSequence + " ;; ");
                        proteinList.Add(peptide.Protein);
                    }
                }
                System.Console.WriteLine();
            }

            //prints parsimonious dictionary
            System.Console.WriteLine("----Parsimonious Dictionary----");
            System.Console.WriteLine("PEPTIDE\t\t\tPROTEIN");
            foreach (var kvp in dictionary)
            {
                proteinList = new List<Protein>();
                System.Console.Write(string.Join("", kvp.Key.BaseSequence.Select(b => char.ConvertFromUtf32(b))) + "  \t\t\t  ");
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinList.Contains(peptide.Protein))
                    {
                        System.Console.Write(peptide.Protein.BaseSequence + " ;; ");
                        proteinList.Add(peptide.Protein);
                    }
                }
                System.Console.WriteLine();
            }

            //prints protein groups after scoring /
            System.Console.WriteLine();
            System.Console.WriteLine("ProteinGroups:");
            foreach (var proteinGroup in proteinGroups)
            {
                System.Console.WriteLine(proteinGroup);
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

            var protease = new Protease("tryp", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var peptides = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString()));

            foreach (var protein in p)
            {
                var digestedProtein = protein.Digest(protease, 2, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>());

                foreach (var pepWithPossibleMods in digestedProtein)
                {
                    var pepWithSetMods = pepWithPossibleMods.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4098, 3);

                    foreach (var peptide in pepWithSetMods)
                        peptides.Add(peptide);
                }
            }

            var CfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var ZdotfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var BfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var YfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();
            var BYfragmentMasses = new Dictionary<PeptideWithSetModifications, double[]>();

            foreach (var peptide in peptides)
            {
                CfragmentMasses.Add(peptide, peptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.C }));
                ZdotfragmentMasses.Add(peptide, peptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.Zdot }));
                BfragmentMasses.Add(peptide, peptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B }));
                YfragmentMasses.Add(peptide, peptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.Y }));
                BYfragmentMasses.Add(peptide, peptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }));
            }
            double[] testB;
            Assert.That(BfragmentMasses.TryGetValue(peptides.First(), out testB));

            double[] testY;
            Assert.That(YfragmentMasses.TryGetValue(peptides.First(), out testY));

            double[] testC;
            Assert.That(CfragmentMasses.TryGetValue(peptides.First(), out testC));

            double[] testZ;
            Assert.That(ZdotfragmentMasses.TryGetValue(peptides.First(), out testZ));
        }

        [Test]
        public static void TestQuantification()
        {
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            int charge = 3;
            double intensity = 1000.0;
            double rt = 20.0;

            // creates some test proteins, digest, and fragment
            string sequence = "NVLIFDLGGGTFDVSILTIEDGIFEVK";
            var protease = new Protease("tryp", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var prot = (new Protein(sequence, "TestProtein"));

            var digestedProtein = prot.Digest(protease, 2, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>());
            var peptide = digestedProtein.First().GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4098, 3).First();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(peptide, charge, intensity, rt);

            var psms = new List<NewPsmWithFdr>();

            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> dfkj = new MzmlScanWithPrecursor(0, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfkj, new MzPeak(2, 2), 1, "TestDataFile");
            var psm = new PsmClassic(peptide, 0, 0, 0, scan);

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>
            {
                {psm.GetCompactPeptide(modsDictionary), new HashSet<PeptideWithSetModifications>{ peptide} }
            };

            psm.GetProteinLinkedInfo(compactPeptideToProteinPeptideMatching, fragmentTolerance, lp, modsDictionary);

            psms.Add(new NewPsmWithFdr(psm));

            FlashLFQEngine FlashLfqEngine = new FlashLFQEngine();
            FlashLfqEngine.PassFilePaths(new string[] { "TestDataFile" });

            AnalysisEngine ae = new AnalysisEngine(new PsmParent[0][], null, new List<Protein>(), null, null, null, null, null, null, null, null, null, null, false, false, false, 0, null, null, 0, false, new List<ProductType> { ProductType.B, ProductType.Y }, double.NaN, InitiatorMethionineBehavior.Variable, new List<string>(), modsDictionary, myMsDataFile, new List<string> { "TestMsFile" });
            //ae.RunQuantification(psms, 10);

            //var theIntensity = psms.First().thisPSM.QuantIntensity[0];
            //Assert.AreEqual(0, theIntensity);
        }

        [Test]
        public static void TestPTMOutput()
        {
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

            ModificationMotif motif;
            ModificationMotif.TryGetMotif("S", out motif);
            variableModifications.Add(new ModificationWithMassAndCf("resMod", null, motif, ModificationSites.Any, ChemicalFormula.ParseFormula("H"), PeriodicTable.GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, "HaHa"));

            var proteinList = new List<Protein> { new Protein("MNNNSKQQQ", "accession") };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>
            {
                {variableModifications.Last(), 1 }
            };

            PeptideWithPossibleModifications modPep = proteinList.First().Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).Last();
            HashSet<PeptideWithSetModifications> value = new HashSet<PeptideWithSetModifications> { modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide1 = new CompactPeptide(value.First(), modsDictionary);
            Assert.AreEqual("QQQ", value.First().Sequence);

            PeptideWithPossibleModifications modPep2 = proteinList.First().Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).First();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { modPep2.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide2 = new CompactPeptide(value2.First(), modsDictionary);
            Assert.AreEqual("MNNNSK", value2.First().Sequence);
            HashSet<PeptideWithSetModifications> value2mod = new HashSet<PeptideWithSetModifications> { modPep2.GetPeptidesWithSetModifications(variableModifications, 4096, 3).Last() };

            CompactPeptide compactPeptide2mod = new CompactPeptide(value2mod.Last(), modsDictionary);
            Assert.AreEqual("MNNNS[HaHa:resMod]K", value2mod.Last().Sequence);

            PeptideWithPossibleModifications modPep3 = proteinList.First().Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).ToList()[1];
            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide3 = new CompactPeptide(value3.First(), modsDictionary);
            Assert.AreEqual("NNNSK", value3.First().Sequence);
            HashSet<PeptideWithSetModifications> value3mod = new HashSet<PeptideWithSetModifications> { modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).Last() };

            CompactPeptide compactPeptide3mod = new CompactPeptide(value3mod.Last(), modsDictionary);
            Assert.AreEqual("NNNS[HaHa:resMod]K", value3mod.Last().Sequence);

            var peptideList = new HashSet<PeptideWithSetModifications>();
            foreach (var protein in proteinList)
            {
                var temp = protein.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>());
                foreach (var dbPeptide in temp)
                {
                    var pepWithSetMods = dbPeptide.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
                    foreach (var peptide in pepWithSetMods)
                    {
                        peptideList.Add(peptide);
                    }
                }
            }

            compactPeptideToProteinPeptideMatching.Add(compactPeptide1, value);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide2, value2);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide3, value3);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide2mod, value2mod);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide3mod, value3mod);

            AnalysisEngine engine = new AnalysisEngine(new PsmParent[0][], compactPeptideToProteinPeptideMatching, new List<Protein>(), null, null, null, null, null, null, null, null, null, null, true, true, true, 0, null, null, 0, false, new List<ProductType> { ProductType.B, ProductType.Y }, double.NaN, InitiatorMethionineBehavior.Variable, new List<string>(), modsDictionary, null, null);

            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            proteinGroups = engine.ConstructProteinGroups(new HashSet<PeptideWithSetModifications>(), peptideList);

            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> jdfk = new MzmlScanWithPrecursor(0, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null, null);
            Ms2ScanWithSpecificMass ms2scan = new Ms2ScanWithSpecificMass(jdfk, new MzPeak(2, 2), 0, "File");

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);

            var match1 = new PsmClassic(peptideList.ElementAt(0), 0, 10, 0, ms2scan);
            var match2 = new PsmClassic(peptideList.ElementAt(1), 0, 10, 0, ms2scan);

            match1.GetProteinLinkedInfo(compactPeptideToProteinPeptideMatching, fragmentTolerance, lp, modsDictionary);
            match2.GetProteinLinkedInfo(compactPeptideToProteinPeptideMatching, fragmentTolerance, lp, modsDictionary);

            List<NewPsmWithFdr> psms = new List<NewPsmWithFdr>
            {
                new NewPsmWithFdr(match1),
                new NewPsmWithFdr(match2),
                new NewPsmWithFdr(match2)
            };
            engine.ScoreProteinGroups(proteinGroups, psms);
            Assert.AreEqual("#aa5[resMod,info:occupancy=0.67(2/3)];", proteinGroups.First().ModsInfo[0]);
        }

        [Test]
        public static void TestFlashLFQ()
        {
            FlashLFQEngine e = new FlashLFQEngine();
            Assert.That(e != null);
            Assert.That(e.ReadPeriodicTable());
        }

        #endregion Public Methods

    }
}