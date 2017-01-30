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
            string[] sequences = { "AAKBBK", "BBKCCKDDK", "BBKCCKDDKEEK", "GGK", "GGKHHK", "HHKIIK", "IIKJJK", "LLK" };

            IEnumerable<string> sequencesInducingCleavage = new List<string> { "K", "R" };
            IEnumerable<string> sequencesPreventingCleavage = new List<string> { "KP", "RP" };
            var temp1 = new Dictionary<int, List<MetaMorpheusModification>>();
            var temp2 = new List<MetaMorpheusModification>();
            int[] temp3 = new int[0];
            var protease = new Protease("Trypsin", sequencesInducingCleavage, sequencesPreventingCleavage, TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], i.ToString(), temp1, temp3, temp3, null, "Test" + i.ToString(), "FullTest" + i.ToString(), 0, false, false));
            p.Add(new Protein("CCKEEK", "D", temp1, temp3, temp3, null, "Decoy ", "Decoy ", 0, true, false));

            // list of "detected" peptides
            IEnumerable<PeptideWithPossibleModifications> temp;
            IEnumerable<PeptideWithSetModifications> pepWithSetMods = null;
            foreach (var protein in p)
            {
                temp = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable);

                foreach (var dbPeptide in temp)
                {
                    pepWithSetMods = dbPeptide.GetPeptideWithSetModifications(temp2, 4098, 3, new List<MetaMorpheusModification>());
                    foreach (var peptide in pepWithSetMods)
                    {
                        switch (peptide.BaseSequence)
                        {
                            case "AAK": peptideList.Add(peptide); break;
                            case "BBK": peptideList.Add(peptide); break;
                            case "CCK": peptideList.Add(peptide); break;
                            case "DDK": peptideList.Add(peptide); break;
                            case "GGK": peptideList.Add(peptide); break;
                            case "HHK": peptideList.Add(peptide); break;
                            case "IIK": peptideList.Add(peptide); break;
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
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), temp2, temp2, new List<MetaMorpheusModification>());
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
                        case "AAK": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 10), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "BBK": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 9), hashSet, null, null, null), 1, 0, 0.1)); break;
                        case "CCK": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 8), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "DDK": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 7), hashSet, null, null, null), 1, 0, 0.2)); break;
                        case "GGK": psms.Add(new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new PsmClassic(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 6), hashSet, null, null, null), 1, 0, 0.3)); break;
                    }
                }
            }

            ae.ScoreProteinGroups(proteinGroups, psms);
            ae.DoProteinFdr(proteinGroups);

            /*
            // prints initial dictionary
            List<Protein> proteinList = new List<Protein>();
            System.Console.WriteLine("----Initial Dictionary----");
            System.Console.WriteLine("PEPTIDE\t\t\tPROTEIN\t\t\tPeptideWithSetModifications");
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
            System.Console.WriteLine("PEPTIDE\t\t\tPROTEIN\t\t\tPeptideWithSetModifications");
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
            Assert.That(parsimonyProteinList.Count == 6);
            Assert.That(parsimonyBaseSequences.Contains("AAKBBK"));
            Assert.That(parsimonyBaseSequences.Contains("CCKEEK"));
            Assert.That(parsimonyBaseSequences.Contains("GGKHHK"));
            Assert.That(parsimonyBaseSequences.Contains("BBKCCKDDK"));
            Assert.That(parsimonyBaseSequences.Contains("HHKIIK"));
            Assert.That(parsimonyBaseSequences.Contains("BBKCCKDDKEEK"));
            Assert.That(parsimonyBaseSequences.Count() == 6);

            // protein group tests
            Assert.That(proteinGroups.Count == 2);
            Assert.That(proteinGroups.First().proteinGroupScore == 10);

            // sequence coverage test
            foreach (var proteinGroup in proteinGroups)
                foreach (var coverage in proteinGroup.sequenceCoverage)
                    Assert.That(coverage <= 1.0);
        }

        #endregion Public Methods

    }
}