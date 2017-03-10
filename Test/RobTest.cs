using Chemistry;
using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
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
                                   "--C-------",   // 2: contains unique
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
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], null, "", "", false, false, null));
            p.Add(new Protein("-----F----*", "D1", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], null, "", "", true, false, null));
            p.Add(new Protein("-----F----**", "C1", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], null, "", "", false, true, null));
            p.Add(new Protein("----E----**", "C2", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], null, "", "", false, true, null));

            IEnumerable<PeptideWithPossibleModifications> temp;
            IEnumerable<PeptideWithSetModifications> pepWithSetMods = null;
            foreach (var protein in p)
            {
                temp = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>());

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

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<ModificationWithMass>());
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
            AnalysisEngine ae = new AnalysisEngine(new PsmParent[0][], dictionary, new List<Protein>(), null, null, null, null, null, null, null, null, null, null, true, true, 0, 0, false, new List<ProductType> { ProductType.B, ProductType.Y }, double.NaN, InitiatorMethionineBehavior.Variable, new List<string>(), false, 0, 0);
            dictionary = ae.ApplyProteinParsimony(out proteinGroups);

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

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    HashSet<PeptideWithSetModifications> hashSet = new HashSet<PeptideWithSetModifications>();
                    hashSet.Add(peptide);

                    switch (peptide.BaseSequence)
                    {
                        case "A": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0), hashSet, null, null, null))); break;
                        case "B": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0), hashSet, null, null, null))); break;
                        case "C": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0), hashSet, null, null, null))); break;
                        case "D": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0), hashSet, null, null, null))); break;
                        case "E": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0), hashSet, null, null, null))); break;
                        case "F": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0), hashSet, null, null, null))); break;
                        case "G": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0), hashSet, null, null, null))); break;
                        case "H": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0), hashSet, null, null, null))); break;
                        case "I": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0), hashSet, null, null, null))); break;
                    }
                }
            }

            ae.ScoreProteinGroups(proteinGroups, psms);
            proteinGroups = ae.DoProteinFdr(proteinGroups);

            /*
            // prints initial dictionary
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
                        System.Console.Write(peptide.Protein.BaseSequence + " ;; ");
                        proteinList.Add(peptide.Protein);
                    }
                }
                System.Console.WriteLine();
            }

            // prints parsimonious dictionary
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

            // prints protein groups after scoring/fdr
            System.Console.WriteLine(ProteinGroup.TabSeparatedHeader);
            foreach (var proteinGroup in proteinGroups)
            {
                System.Console.WriteLine(proteinGroup);
            }
            */

            // check that correct proteins are in parsimony list
            Assert.That(parsimonyProteinList.Count == 8);
            Assert.That(parsimonyBaseSequences.Contains("AB--------"));
            Assert.That(parsimonyBaseSequences.Contains("--C-------"));
            Assert.That(parsimonyBaseSequences.Contains("-B-D---HHH--"));
            Assert.That(parsimonyBaseSequences.Contains("-----F----*"));  // decoy
            Assert.That(parsimonyBaseSequences.Contains("----E----**"));  // contaminant
            Assert.That(parsimonyBaseSequences.Contains("-B------I-"));
            Assert.That(parsimonyBaseSequences.Contains("----EFG---"));
            Assert.That(parsimonyBaseSequences.Contains("----EFG--J"));
            Assert.That(parsimonyBaseSequences.Count == 8);

            // protein group tests
            Assert.That(proteinGroups.Count == 5);
            Assert.That(proteinGroups.First().AllPsmsForStrictPeptideSequences.Count == 3);
            Assert.That(proteinGroups.First().proteinGroupScore == 19);

            // sequence coverage test
            foreach (var proteinGroup in proteinGroups)
                foreach (var coverage in proteinGroup.sequenceCoveragePercent)
                    Assert.That(coverage <= 1.0);
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
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], null, "", "", false, false, null));

            foreach (var protein in p)
            {
                var digestedProtein = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>());

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
                CfragmentMasses.Add(peptide, peptide.FastSortedProductMasses(new List<ProductType> { ProductType.C }));
                ZdotfragmentMasses.Add(peptide, peptide.FastSortedProductMasses(new List<ProductType> { ProductType.Zdot }));
                BfragmentMasses.Add(peptide, peptide.FastSortedProductMasses(new List<ProductType> { ProductType.B }));
                YfragmentMasses.Add(peptide, peptide.FastSortedProductMasses(new List<ProductType> { ProductType.Y }));
                BYfragmentMasses.Add(peptide, peptide.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }));
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
            int charge = 3;
            double mass = 2910.52664 + charge * Constants.protonMass;
            double mz = mass / charge;
            double intensity = 1000.0;
            double rt = 20.0;

            // creates some test proteins, digest, and fragment
            string sequence = "NVLIFDLGGGTFDVSILTIEDGIFEVK";
            var protease = new Protease("tryp", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var prot = (new Protein(sequence, "", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], null, "", "", false, false, null));
            var digestedProtein = prot.Digest(protease, 2, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>());
            var peptide = digestedProtein.First().GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4098, 3).First();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(peptide, charge, intensity, rt);

            var psms = new List<NewPsmWithFdr>();

            var psm = new PsmClassic(peptide, null, rt, intensity, mass, 2, 1, charge, 1, 0, mz, 0, 0);
            var t = new PsmWithMultiplePossiblePeptides(psm, new HashSet<PeptideWithSetModifications> { peptide }, null, null, null);
            psms.Add(new NewPsmWithFdr(t));

            AnalysisEngine ae = new AnalysisEngine(new PsmParent[0][], null, new List<Protein>(), null, null, null, null, null, myMsDataFile, null, null, null, null, true, true, 0, 0, false, new List<ProductType> { ProductType.B, ProductType.Y }, double.NaN, InitiatorMethionineBehavior.Variable, new List<string>(), false, 0, 0);
            ae.RunQuantification(psms, 0.2, 10);

            Assert.That(psms.First().thisPSM.newPsm.apexIntensity == 1000.0);
            Assert.That(double.IsNaN(psms.First().thisPSM.newPsm.fullWidthHalfMax));
        }

        #endregion Public Methods

    }
}