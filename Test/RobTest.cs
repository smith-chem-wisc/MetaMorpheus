using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using NUnit.Framework;

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
                                   "---D------",   // 3: subset
                                   "-B-D------",   // 4: D should go to 4, not 3
                                   "-B--E-----",   // 5: subsumable
                                   "----EFG---",   // 6: indistinguishable from 8 (J will not be a "detected" PSM)
                                   "-----F----",   // 7: only pep shared w/ decoy
                                   "--------I-",   // 8: HI should go to 9, not 8
                                   "-B------I-",   // 9: HI should go to 9, not 8
                                   "----EFG--J" }; // 10: indistinguishable from 6 (J will not be a "detected" PSM)

            IEnumerable<string> sequencesInducingCleavage = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "-" };
            var protease = new Protease("test", sequencesInducingCleavage, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString(), new Dictionary<int, List<MetaMorpheusModification>>(), new int[0], new int[0], null, "", "", 0, false, false));
            p.Add(new Protein("-----F----*", "D", new Dictionary<int, List<MetaMorpheusModification>>(), new int[0], new int[0], null, "", "", 0, true, false));

            IEnumerable<PeptideWithPossibleModifications> temp;
            IEnumerable<PeptideWithSetModifications> pepWithSetMods = null;
            foreach (var protein in p)
            {
                temp = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable);

                foreach (var dbPeptide in temp)
                {
                    pepWithSetMods = dbPeptide.GetPeptideWithSetModifications(new List<MetaMorpheusModification>(), 4098, 3, new List<MetaMorpheusModification>());
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
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count()];
            HashSet<PeptideWithSetModifications>[] virtualPeptideSets = new HashSet<PeptideWithSetModifications>[peptideList.Count()];

            // creates peptide list
            for (int i = 0; i < peptideList.Count(); i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), new List<MetaMorpheusModification>(), new List<MetaMorpheusModification>(), new List<MetaMorpheusModification>());
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
            AnalysisEngine ae = new AnalysisEngine(new PsmParent[0][], dictionary, new List<Protein>(), null, null, null, null, null, null, null, null, null, null, true, 0, 0, false, new List<ProductType> { ProductType.B, ProductType.Y }, double.NaN);
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
                        case "A": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 10), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "B": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 9), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "C": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 8), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "D": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 7), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "E": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 6), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "F": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 5), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "G": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 4), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "H": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 3), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "I": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 2), hashSet, null, null, null), 1, 0, 0.0)); break;
                    }
                }
            }

            ae.ScoreProteinGroups(proteinGroups, psms);
            ae.DoProteinFdr(proteinGroups);
            
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
            Assert.That(parsimonyProteinList.Count == 7);
            Assert.That(parsimonyBaseSequences.Contains("AB--------"));
            Assert.That(parsimonyBaseSequences.Contains("--C-------"));
            Assert.That(parsimonyBaseSequences.Contains("-B-D------"));
            Assert.That(parsimonyBaseSequences.Contains("----EFG---"));
            Assert.That(parsimonyBaseSequences.Contains("-----F----*"));
            Assert.That(parsimonyBaseSequences.Contains("-B------I-"));
            Assert.That(parsimonyBaseSequences.Contains("----EFG--J"));
            Assert.That(parsimonyBaseSequences.Count() == 7);

            // protein group tests
            Assert.That(proteinGroups.Count == 6);
            Assert.That(proteinGroups.First().AllPsmsForStrictPeptideSequences.Count() == 2);
            Assert.That(proteinGroups.First().proteinGroupScore == 19);

            // sequence coverage test
            foreach (var proteinGroup in proteinGroups)
                foreach (var coverage in proteinGroup.sequenceCoverage)
                    Assert.That(coverage <= 1.0);
                    
        }

        #endregion Public Methods
    }
}