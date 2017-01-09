using IndexSearchAndAnalyze;
using MassSpectrometry;
using MetaMorpheus;
using NUnit.Framework;
using Proteomics;
using Spectra;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class AnalysisEngineTest
    {
        [Test]
        public void TestAnalysis()
        {
            List<NewPsm>[] newPsms = null;
            Dictionary<CompactPeptide, ConcurrentDictionary<PeptideWithSetModifications, byte>> compactPeptideToProteinPeptideMatching = null;
            List<Protein> proteinList = null;
            List<MorpheusModification> variableModifications = null;
            List<MorpheusModification> fixedModifications = null;
            List<MorpheusModification> localizeableModifications = null;
            IMsDataFile<IMzSpectrum<MzPeak>> iMsDataFile = null;
            Protease protease = null;
            List<SearchMode> searchModes = null;
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            UsefulProteomicsDatabases.Generated.unimod unimodDeserialized = null;
            Dictionary<int, ChemicalFormulaModification> uniprotDeseralized = null;

            AllTasksParams po = null;
            AnalysisParams analysisParams = new AnalysisParams(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, iMsDataFile, fragmentTolerance, (BinTreeStructure myTreeStructure, string s) => { }, (List<NewPsmWithFDR> h, string s) => { }, po);
            AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);

            Assert.That(() => analysisEngine.Run(), Throws.TypeOf<ValidationException>()
                    .With.Property("Message").EqualTo("newPsms is null"));

            newPsms = new List<NewPsm>[1];
            newPsms[0] = new List<NewPsm>();

            analysisParams = new AnalysisParams(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, iMsDataFile, fragmentTolerance, (BinTreeStructure myTreeStructure, string s) => { }, (List<NewPsmWithFDR> h, string s) => { }, po);
            analysisEngine = new AnalysisEngine(analysisParams);
            Assert.That(() => analysisEngine.Run(), Throws.TypeOf<ValidationException>()
                    .With.Property("Message").EqualTo("proteinList is null"));
        }

        [Test]
        public void TestParsimony()
        {
            // creates some test proteins and digests them (simulating a protein database)
            string sequence1 = "AKCKBK";
            string sequence2 = "DKCK";
            string sequence3 = "AAKAAK";

            IEnumerable<string> sequencesInducingCleavage = new List<string>() { "K", "R" };
            IEnumerable<string> sequencesPreventingCleavage = new List<string>() { "KP", "RP" };
            Dictionary<int, List<MorpheusModification>> temp1 = new Dictionary<int, List<MorpheusModification>>();
            List<MorpheusModification> temp2 = new List<MorpheusModification>();
            int[] temp3 = new int[0];
            Protease protease = new Protease("Trypsin", sequencesInducingCleavage, sequencesPreventingCleavage, MetaMorpheus.Terminus.C, CleavageSpecificity.Full, null, null, null);
            HashSet<PeptideWithSetModifications> totalProteinList = new HashSet<PeptideWithSetModifications>();

            Protein p1 = new Protein(sequence1, "1", null, temp1, temp3, temp3, null, "Test1", "TestFullName1", 0, false);
            Protein p2 = new Protein(sequence2, "2", null, temp1, temp3, temp3, null, "Test2", "TestFullName2", 0, false);
            Protein p3 = new Protein(sequence3, "3", null, temp1, temp3, temp3, null, "Test3", "TestFullName3", 0, false);

            IEnumerable<PeptideWithPossibleModifications> digestedList1 = p1.Digest(protease, 2, InitiatorMethionineBehavior.Variable);
            IEnumerable<PeptideWithPossibleModifications> digestedList2 = p2.Digest(protease, 2, InitiatorMethionineBehavior.Variable);
            IEnumerable<PeptideWithPossibleModifications> digestedList3 = p3.Digest(protease, 2, InitiatorMethionineBehavior.Variable);

            foreach (var protein in digestedList1)
            {
                IEnumerable<PeptideWithSetModifications> peptides1 = protein.GetPeptideWithSetModifications(temp2, 4098, 3, temp2);

                foreach (var peptide in peptides1)
                    totalProteinList.Add(peptide);
            }

            foreach (var protein in digestedList2)
            {
                IEnumerable<PeptideWithSetModifications> peptides2 = protein.GetPeptideWithSetModifications(temp2, 4098, 3, temp2);

                foreach (var peptide in peptides2)
                    totalProteinList.Add(peptide);
            }

            foreach (var protein in digestedList3)
            {
                IEnumerable<PeptideWithSetModifications> peptides3 = protein.GetPeptideWithSetModifications(temp2, 4098, 3, temp2);

                foreach (var peptide in peptides3)
                    totalProteinList.Add(peptide);
            }

            // creates the initial dictionary of "peptide" and "protein" matches (protein must contain peptide sequence)
            Dictionary<CompactPeptide, ConcurrentDictionary<PeptideWithSetModifications, byte>> initialDictionary = new Dictionary<CompactPeptide, ConcurrentDictionary<PeptideWithSetModifications, byte>>();
            CompactPeptide[] peptides = new CompactPeptide[totalProteinList.Count()];
            ConcurrentDictionary<PeptideWithSetModifications, byte>[] proteinSets = new ConcurrentDictionary<PeptideWithSetModifications, byte>[totalProteinList.Count()];

            // creates peptide list
            for (int i = 0; i < totalProteinList.Count(); i++)
            {
                peptides[i] = new CompactPeptide(totalProteinList.ElementAt(i), temp2, temp2);
            }

            // creates protein list
            for (int i = 0; i < proteinSets.Length; i++)
            {
                proteinSets[i] = new ConcurrentDictionary<PeptideWithSetModifications, byte>();

                foreach (var protein in totalProteinList)
                {
                    string peptideBaseSequence = string.Join("", peptides[i].BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (protein.BaseSequence.Contains(peptideBaseSequence))
                    {
                        proteinSets[i].GetOrAdd(protein, 0);
                        //proteinSets[i].Add(protein);
                    }
                }
            }

            // populates initial peptide-protein dictionary
            for (int i = 0; i < peptides.Length; i++)
            {
                if (!initialDictionary.ContainsKey(peptides[i]))
                {
                    initialDictionary.Add(peptides[i], proteinSets[i]);
                }
            }

            // prints initial dictionary
            AnalysisParams analysisParams = new AnalysisParams(new List<NewPsm>(), null, null, null, null, null, null, null, null, null, null, null, null);
            var analysisEngine = new AnalysisEngine(analysisParams);

            // apply parsimony to initial dictionary
            var parsimonyTest = analysisEngine.ApplyProteinParsimony(initialDictionary);

            // apply the single pick version to parsimonious dictionary
            var singlePickTest = AnalysisEngine.GetSingleMatchDictionary(parsimonyTest);
        }
    }
}