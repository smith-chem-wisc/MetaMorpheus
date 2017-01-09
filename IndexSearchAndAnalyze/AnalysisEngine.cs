using MetaMorpheus;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace IndexSearchAndAnalyze
{
    public class AnalysisEngine : MyEngine
    {
        public AnalysisEngine(AnalysisParams analysisParams)
        {
            this.myParams = analysisParams;
        }

        protected override MyResults RunSpecific()
        {
            myParams.allTasksParams.status("Running analysis engine!");
            var analysisParams = (AnalysisParams)myParams;
            myParams.allTasksParams.status("Adding observed peptides to dictionary...");
            AddObservedPeptidesToDictionary(analysisParams.newPsms, analysisParams.compactPeptideToProteinPeptideMatching, analysisParams.proteinList, analysisParams.variableModifications, analysisParams.fixedModifications, analysisParams.localizeableModifications, analysisParams.protease);
            myParams.allTasksParams.status("Getting protein parsimony dictionary...");
            var parsimonyDictionary = ApplyProteinParsimony(analysisParams.compactPeptideToProteinPeptideMatching);
            myParams.allTasksParams.status("Getting single match dictionary...");
            var fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(parsimonyDictionary);
            //myParams.allTasksParams.status("Getting single match dictionary...");
            //var fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(analysisParams.compactPeptideToProteinPeptideMatching);

            List<NewPsmWithFDR>[] yeah = new List<NewPsmWithFDR>[analysisParams.searchModes.Count];
            for (int j = 0; j < analysisParams.searchModes.Count; j++)
            {
                PSMwithPeptide[] psmsWithPeptides = new PSMwithPeptide[analysisParams.newPsms[0].Count];

                Parallel.ForEach(Partitioner.Create(0, analysisParams.newPsms[0].Count), fff =>
                {
                    for (int i = fff.Item1; i < fff.Item2; i++)
                    {
                        if (analysisParams.newPsms[j] != null)
                        {
                            var huh = analysisParams.newPsms[j][i];
                            if (huh != null)
                                if (huh.ScoreFromSearch >= 1)
                                    psmsWithPeptides[i] = new PSMwithPeptide(huh, fullSequenceToProteinSingleMatch[huh.peptide], analysisParams.fragmentTolerance, analysisParams.myMsDataFile);
                        }
                    }
                });

                var orderedPsmsWithPeptides = psmsWithPeptides.Where(b => b != null).OrderByDescending(b => b.ScoreFromSearch);

                var orderedPsmsWithPeptidesAndFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides);
                var limitedpsms_with_fdr = orderedPsmsWithPeptidesAndFDR.Where(b => (b.QValue <= 0.01)).ToList();
                if (limitedpsms_with_fdr.Where(b => !b.isDecoy).Count() > 0)
                {
                    myParams.allTasksParams.status("Running histogram analysis...");
                    var hm = MyAnalysis(limitedpsms_with_fdr, analysisParams.unimodDeserialized, analysisParams.uniprotDeseralized);
                    analysisParams.action1(hm, analysisParams.searchModes[j].FileNameAddition);
                }

                analysisParams.action2(orderedPsmsWithPeptidesAndFDR, analysisParams.searchModes[j].FileNameAddition);

                yeah[j] = orderedPsmsWithPeptidesAndFDR;
            }

            return new AnalysisResults(analysisParams, yeah);
        }

        private static BinTreeStructure MyAnalysis(List<NewPsmWithFDR> limitedpsms_with_fdr, UsefulProteomicsDatabases.Generated.unimod unimodDeserialized, Dictionary<int, ChemicalFormulaModification> uniprotDeseralized)
        {
            BinTreeStructure myTreeStructure = new BinTreeStructure();
            myTreeStructure.GenerateBins(limitedpsms_with_fdr, 0.003);
            myTreeStructure.AddToBins(limitedpsms_with_fdr);

            MyAnalysisClass.IdentifyUnimodBins(myTreeStructure, 0.003, unimodDeserialized);
            MyAnalysisClass.IdentifyUniprotBins(myTreeStructure, 0.003, uniprotDeseralized);
            MyAnalysisClass.IdentifyAA(myTreeStructure, 0.003);

            MyAnalysisClass.IdentifyCombos(myTreeStructure, 0.003);

            MyAnalysisClass.IdentifyResidues(myTreeStructure);

            MyAnalysisClass.IdentifyMods(myTreeStructure);

            MyAnalysisClass.IdentifyAAsInCommon(myTreeStructure);

            MyAnalysisClass.IdentifyMine(myTreeStructure, 0.003);

            return myTreeStructure;
        }

        private static List<NewPsmWithFDR> DoFalseDiscoveryRateAnalysis(IEnumerable<PSMwithPeptide> items)
        {
            List<NewPsmWithFDR> ids = new List<NewPsmWithFDR>();

            int cumulative_target = 0;
            int cumulative_decoy = 0;
            foreach (PSMwithPeptide item in items)
            {
                var isDecoy = item.isDecoy;
                if (isDecoy)
                    cumulative_decoy++;
                else
                    cumulative_target++;
                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                ids.Add(new NewPsmWithFDR(item, cumulative_target, cumulative_decoy, temp_q_value));
            }

            double min_q_value = double.PositiveInfinity;
            for (int i = ids.Count - 1; i >= 0; i--)
            {
                NewPsmWithFDR id = ids[i];
                if (id.QValue > min_q_value)
                    id.QValue = min_q_value;
                else if (id.QValue < min_q_value)
                    min_q_value = id.QValue;
            }

            return ids;
        }

        private static void AddObservedPeptidesToDictionary(List<NewPsm>[] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease)
        {
            foreach (var ah in newPsms)
            {
                if (ah != null)
                    foreach (var fhh in ah)
                    {
                        if (fhh != null && !fullSequenceToProteinPeptideMatching.ContainsKey(fhh.peptide))
                            fullSequenceToProteinPeptideMatching.Add(fhh.peptide, new HashSet<PeptideWithSetModifications>());
                    }
            }

            foreach (var protein in proteinList)
                foreach (var peptide in protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable).ToList())
                {
                    if (peptide.Length == 1 || peptide.Length > 252)
                        continue;
                    peptide.SetFixedModifications(fixedModifications);
                    var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, 4098, 3, localizeableModifications).ToList();
                    foreach (var yyy in ListOfModifiedPeptides)
                    {
                        HashSet<PeptideWithSetModifications> v;
                        if (fullSequenceToProteinPeptideMatching.TryGetValue(new CompactPeptide(yyy, variableModifications, localizeableModifications), out v))
                        {
                            v.Add(yyy);
                        }
                    }
                }
        }

        public static Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> ApplyProteinParsimony(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // makes dictionary with proteins as keys and list of associated peptides as the value (swaps input parameter dictionary keys/values)
            Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>> newDict = new Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>>();
            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                foreach(var protein in kvp.Value)
                {
                    if (!newDict.ContainsKey(protein))
                    {
                        HashSet<CompactPeptide> peptides = new HashSet<CompactPeptide>();

                        foreach (var kvp1 in fullSequenceToProteinPeptideMatching)
                        {
                            if(kvp1.Value.Contains(protein))
                            {
                                peptides.Add(kvp1.Key);
                            }
                        }

                        newDict.Add(protein, peptides);
                    }
                }
            }
            
            // greedy algorithm adds the next protein that will account for the most unique peptides
            Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>> parsimonyDict = new Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>>();
            HashSet<CompactPeptide> usedPeptides = new HashSet<CompactPeptide>();
            HashSet<CompactPeptide> bestProteinPeptideList = new HashSet<CompactPeptide>();
            PeptideWithSetModifications bestProtein = null;
            int currentBestNumNewPeptides = -1;

            // as long as there are peptides that have not been accounted for, keep going
            while(currentBestNumNewPeptides != 0)
            {
                currentBestNumNewPeptides = 0;

                // attempt to find protein that best accounts for unaccounted-for peptides
                foreach (var kvp in newDict)
                {
                    int comparisonProteinNewPeptides = 0;

                    // determines number of unaccounted-for peptides for the current protein
                    foreach (CompactPeptide peptide in kvp.Value)
                    {
                        if(!usedPeptides.Contains(peptide))
                        {
                            comparisonProteinNewPeptides++;
                        }
                    }

                    // if the current protein is better than the best so far, current protein is the new best protein
                    if (comparisonProteinNewPeptides > currentBestNumNewPeptides)
                    {
                        bestProtein = kvp.Key;
                        bestProteinPeptideList = kvp.Value;
                        currentBestNumNewPeptides = comparisonProteinNewPeptides;
                    }
                }

                // adds the best protein if algo found unaccounted-for peptides
                if(currentBestNumNewPeptides > 0)
                {
                    parsimonyDict.Add(bestProtein, bestProteinPeptideList);

                    // the new best protein's peptides have been accounted for
                    foreach (CompactPeptide peptide in bestProteinPeptideList)
                    {
                        usedPeptides.Add(peptide);
                    }
                }
            }
            
            // swaps keys and values back for return
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> answer = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!answer.ContainsKey(peptide))
                    {
                        HashSet<PeptideWithSetModifications> proteins = new HashSet<PeptideWithSetModifications>();

                        foreach (var kvp1 in parsimonyDict)
                        {
                            if (kvp1.Value.Contains(peptide))
                            {
                                proteins.Add(kvp1.Key);
                            }
                        }

                        answer.Add(peptide, proteins);
                    }
                }
            }

            return answer;
        }

        public static Dictionary<CompactPeptide, PeptideWithSetModifications> GetSingleMatchDictionary(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // Right now very stupid, add the first decoy one, and if no decoy, add the first one
            Dictionary<CompactPeptide, PeptideWithSetModifications> outDict = new Dictionary<CompactPeptide, PeptideWithSetModifications>();
            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                var k = kvp.Key;
                var val = kvp.Value;
                bool sawDecoy = false;
                foreach (var entry in val)
                {
                    if (entry.protein.isDecoy)
                    {
                        outDict[kvp.Key] = entry;
                        sawDecoy = true;
                        break;
                    }
                }
                if (sawDecoy == false)
                    outDict[kvp.Key] = kvp.Value.First();
            }

            return outDict;
        }
    }
}