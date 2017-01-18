using InternalLogicEngineLayer;
using NUnit.Framework;
using OldInternalLogic;
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
            string sequence1 = "AKCKBK";
            string sequence2 = "DKCK";
            string sequence3 = "AAAAK";

            IEnumerable<string> sequencesInducingCleavage = new List<string> { "K", "R" };
            IEnumerable<string> sequencesPreventingCleavage = new List<string> { "KP", "RP" };
            var temp1 = new Dictionary<int, List<MorpheusModification>>();
            var temp2 = new List<MorpheusModification>();
            int[] temp3 = new int[0];
            var protease = new Protease("Trypsin", sequencesInducingCleavage, sequencesPreventingCleavage, OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);
            var totalVirtualPeptideList = new HashSet<PeptideWithSetModifications>();

            var p1 = new Protein(sequence1, "1", null, temp1, temp3, temp3, null, "Test1", "TestFullName1", 0, false);
            var p2 = new Protein(sequence2, "2", null, temp1, temp3, temp3, null, "Test2", "TestFullName2", 0, false);
            var p3 = new Protein(sequence3, "3", null, temp1, temp3, temp3, null, "Test3", "TestFullName3", 0, false);

            IEnumerable<PeptideWithPossibleModifications> digestedList1 = p1.Digest(protease, 2, InitiatorMethionineBehavior.Variable);
            IEnumerable<PeptideWithPossibleModifications> digestedList2 = p2.Digest(protease, 2, InitiatorMethionineBehavior.Variable);
            IEnumerable<PeptideWithPossibleModifications> digestedList3 = p3.Digest(protease, 2, InitiatorMethionineBehavior.Variable);

            foreach (var protein in digestedList1)
            {
                IEnumerable<PeptideWithSetModifications> peptides1 = protein.GetPeptideWithSetModifications(temp2, 4098, 3);

                foreach (var peptide in peptides1)
                    totalVirtualPeptideList.Add(peptide);
            }

            foreach (var protein in digestedList2)
            {
                IEnumerable<PeptideWithSetModifications> peptides2 = protein.GetPeptideWithSetModifications(temp2, 4098, 3);

                foreach (var peptide in peptides2)
                    totalVirtualPeptideList.Add(peptide);
            }

            foreach (var protein in digestedList3)
            {
                IEnumerable<PeptideWithSetModifications> peptides3 = protein.GetPeptideWithSetModifications(temp2, 4098, 3);

                foreach (var peptide in peptides3)
                    totalVirtualPeptideList.Add(peptide);
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var initialDictionary = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[totalVirtualPeptideList.Count()];
            HashSet<PeptideWithSetModifications>[] virtualPeptideSets = new HashSet<PeptideWithSetModifications>[totalVirtualPeptideList.Count()];

            // creates peptide list
            for (int i = 0; i < totalVirtualPeptideList.Count(); i++)
            {
                peptides[i] = new CompactPeptide(totalVirtualPeptideList.ElementAt(i), temp2, temp2);
            }

            // creates protein list
            for (int i = 0; i < virtualPeptideSets.Length; i++)
            {
                virtualPeptideSets[i] = new HashSet<PeptideWithSetModifications>();

                foreach (var virtualPeptide in totalVirtualPeptideList)
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
                if (!initialDictionary.ContainsKey(peptides[i]))
                {
                    initialDictionary.Add(peptides[i], virtualPeptideSets[i]);
                }
            }
            
            HashSet<CompactPeptide> pep = new HashSet<CompactPeptide>();

            // apply parsimony to initial dictionary
            var parsimonyTest = AnalysisEngine.ApplyProteinParsimony(initialDictionary, out pep);

            var parsimonyProteinList = new List<Protein>();
            string[] parsimonyBaseSequences = new string[3];
            int j = 0;

            foreach (var kvp in parsimonyTest)
            {
                foreach (var virtualPeptide in kvp.Value)
                {
                    if (!parsimonyProteinList.Contains(virtualPeptide.Protein))
                    {
                        if (j < 3)
                        {
                            parsimonyProteinList.Add(virtualPeptide.Protein);
                            parsimonyBaseSequences[j] = virtualPeptide.Protein.BaseSequence;
                            j++;
                        }
                    }
                }
            }

            /*
            // prints initial dictionary
            List<Protein> proteinList = new List<Protein>();

            Console.WriteLine("----Initial Dictionary----");
            Console.WriteLine("PEPTIDE\t\t\tPROTEIN\t\t\tPeptideWithSetModifications");
            foreach (var kvp in initialDictionary)
            {
                proteinList = new List<Protein>();
                Console.Write(string.Join("", kvp.Key.BaseSequence.Select(b => char.ConvertFromUtf32(b))) + "  \t\t\t  ");
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinList.Contains(peptide.protein))
                    {
                        Console.Write(peptide.protein.BaseSequence + " ;; ");
                        proteinList.Add(peptide.protein);
                    }
                }
                Console.WriteLine();
            }

            // prints parsimonious dictionary
            Console.WriteLine("----Parsimonious Dictionary----");
            Console.WriteLine("PEPTIDE\t\t\tPROTEIN\t\t\tPeptideWithSetModifications");
            foreach (var kvp in parsimonyTest)
            {
                proteinList = new List<Protein>();
                Console.Write(string.Join("", kvp.Key.BaseSequence.Select(b => char.ConvertFromUtf32(b))) + "  \t\t\t  ");
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinList.Contains(peptide.protein))
                    {
                        Console.Write(peptide.protein.BaseSequence + " ;; ");
                        proteinList.Add(peptide.protein);
                    }
                }
                Console.WriteLine();
            }
            */

            Assert.That(parsimonyProteinList.Count == 3);
            Assert.That(parsimonyBaseSequences.Contains(sequence1));
            Assert.That(parsimonyBaseSequences.Contains(sequence2));
            Assert.That(parsimonyBaseSequences.Contains(sequence3));
        }

        #endregion Public Methods

    }
}
