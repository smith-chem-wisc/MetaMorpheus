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

            var p1 = new Protein(sequence1, "1", temp1, temp3, temp3, null, "Test1", "TestFullName1", 0, false, false);
            var p2 = new Protein(sequence2, "2", temp1, temp3, temp3, null, "DECOY_Test2", "DECOY_TestFullName2", 0, true, false);
            var p3 = new Protein(sequence3, "3", temp1, temp3, temp3, null, "Test3", "TestFullName3", 0, false, false);

            IEnumerable<PeptideWithPossibleModifications> digestedList1 = p1.Digest(protease, 2, InitiatorMethionineBehavior.Variable);
            IEnumerable<PeptideWithPossibleModifications> digestedList2 = p2.Digest(protease, 2, InitiatorMethionineBehavior.Variable);
            IEnumerable<PeptideWithPossibleModifications> digestedList3 = p3.Digest(protease, 2, InitiatorMethionineBehavior.Variable);
            IEnumerable<PeptideWithSetModifications> peptides1 = null;
            IEnumerable<PeptideWithSetModifications> peptides2 = null;
            IEnumerable<PeptideWithSetModifications> peptides3 = null;

            foreach (var protein in digestedList1)
            {
                peptides1 = protein.GetPeptideWithSetModifications(temp2, 4098, 3);

                foreach (var peptide in peptides1)
                    totalVirtualPeptideList.Add(peptide);
            }

            foreach (var protein in digestedList2)
            {
                peptides2 = protein.GetPeptideWithSetModifications(temp2, 4098, 3);

                foreach (var peptide in peptides2)
                    totalVirtualPeptideList.Add(peptide);
            }

            foreach (var protein in digestedList3)
            {
                peptides3 = protein.GetPeptideWithSetModifications(temp2, 4098, 3);

                foreach (var peptide in peptides3)
                    totalVirtualPeptideList.Add(peptide);
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
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
                if (!dictionary.ContainsKey(peptides[i]))
                {
                    dictionary.Add(peptides[i], virtualPeptideSets[i]);
                }
            }
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            AnalysisEngine ae = new AnalysisEngine(null, dictionary, null, null, null, null, null, null, null, null, null, null, null, true, 0, 0, false, new List<ProductType> { ProductType.B, ProductType.Y });

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> initialDictionary = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
            foreach(var kvp in dictionary)
            {
                initialDictionary.Add(kvp.Key, kvp.Value);
            }


            // apply parsimony to initial dictionary
            ae.ApplyProteinParsimony(out proteinGroups);

            var parsimonyProteinList = new List<Protein>();
            string[] parsimonyBaseSequences = new string[3];
            int j = 0;

            foreach (var kvp in dictionary)
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

            

            List<NewPsmWithFdr> psms = new List<NewPsmWithFdr>();

            foreach(var kvp in dictionary)
            {
                foreach(var peptide in kvp.Value)
                {
                    HashSet<PeptideWithSetModifications> hashSet = new HashSet<PeptideWithSetModifications>();
                    hashSet.Add(peptide);

                    switch (peptide.BaseSequence)
                    {
                        case "AK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 1), hashSet, null, null, null), 1, 0, 0.29)); break;
                        case "CK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 9), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "BK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 8), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "AKCK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 7), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "CKBK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 6), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "AKCKBK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 5), hashSet, null, null, null), 1, 0, 0.0)); break;
                        case "DK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 4), hashSet, null, null, null), 1, 0, 0.2)); break;
                        case "DKCK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 3), hashSet, null, null, null), 1, 0, 0.4)); break;
                        case "AAAAK": psms.Add(new NewPsmWithFdr(new PSMwithProteinHashSet(new ClassicSpectrumMatch(peptide, null, 0, 0, 0, 0, 0, 0, 0, 0, 2), hashSet, null, null, null), 1, 0, 0.33)); break;
                    }
                }
            }

            ae.ScoreProteinGroups(proteinGroups, psms);
            ae.DoProteinFdr(proteinGroups);

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


            Assert.That(parsimonyProteinList.Count == 3);
            Assert.That(parsimonyBaseSequences.Contains(sequence1));
            Assert.That(parsimonyBaseSequences.Contains(sequence2));
            Assert.That(parsimonyBaseSequences.Contains(sequence3));
            Assert.That(proteinGroups.Count == 2);
            Assert.That(proteinGroups.First().proteinGroupScore > 10);
        }

        #endregion Public Methods

    }
}