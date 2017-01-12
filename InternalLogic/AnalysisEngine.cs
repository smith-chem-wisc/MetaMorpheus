using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    public class AnalysisEngine : MyEngine
    {
        public ParentSpectrumMatch[][] newPsms { get; private set; }
        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching { get; private set; }
        public List<Protein> proteinList { get; private set; }
        public List<MorpheusModification> variableModifications { get; private set; }
        public List<MorpheusModification> fixedModifications { get; private set; }
        public List<MorpheusModification> localizeableModifications { get; private set; }
        public Protease protease { get; private set; }
        public List<SearchMode> searchModes { get; private set; }
        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile { get; private set; }
        public Tolerance fragmentTolerance { get; private set; }
        public Action<BinTreeStructure, string> action1 { get; private set; }
        public Action<List<NewPsmWithFDR>, string> action2 { get; private set; }
        public bool doParsimony { get; internal set; }

        public AnalysisEngine(ParentSpectrumMatch[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFDR>, string> action2, bool doParsimony)
        {
            this.doParsimony = doParsimony;
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.protease = protease;
            this.searchModes = searchModes;
            this.myMsDataFile = myMsDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.action1 = action1;
            this.action2 = action2;
        }

        protected override MyResults RunSpecific()
        {
            status("Running analysis engine!");
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            status("Adding observed peptides to dictionary...");
            AddObservedPeptidesToDictionary();

            status("Getting single match just for FDR purposes...");
            var fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(compactPeptideToProteinPeptideMatching);

            List<NewPsmWithFDR>[] yeah = new List<NewPsmWithFDR>[searchModes.Count];
            for (int j = 0; j < searchModes.Count; j++)
            {
                PSMwithTargetDecoyKnown[] psmsWithTargetDecoyKnown = new PSMwithTargetDecoyKnown[newPsms[0].Length];

                Parallel.ForEach(Partitioner.Create(0, newPsms[0].Length), fff =>
                {
                    for (int i = fff.Item1; i < fff.Item2; i++)
                    {
                        if (newPsms[j] != null)
                        {
                            var huh = newPsms[j][i];
                            if (huh != null)
                                if (huh.Score >= 1)
                                    psmsWithTargetDecoyKnown[i] = new PSMwithTargetDecoyKnown(huh, fullSequenceToProteinSingleMatch[huh.GetCompactPeptide(variableModifications, localizeableModifications)], fragmentTolerance, myMsDataFile);
                        }
                    }
                });

                var orderedPsmsWithPeptides = psmsWithTargetDecoyKnown.Where(b => b != null).OrderByDescending(b => b.Score);

                var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides);
                var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.QValue <= 0.01)).ToList();
                if (limitedpsms_with_fdr.Where(b => !b.isDecoy).Count() > 0)
                {
                    status("Running histogram analysis...");
                    var hm = MyAnalysis(limitedpsms_with_fdr);
                    action1(hm, searchModes[j].FileNameAddition);
                }

                action2(orderedPsmsWithFDR, searchModes[j].FileNameAddition);

                if (doParsimony)
                {
                    status("Getting protein parsimony dictionary...");
                    compactPeptideToProteinPeptideMatching = ApplyProteinParsimony(compactPeptideToProteinPeptideMatching);
                }
                yeah[j] = orderedPsmsWithFDR;
            }

            return new AnalysisResults(this, yeah, compactPeptideToProteinPeptideMatching);
        }

        private static BinTreeStructure MyAnalysis(List<NewPsmWithFDR> limitedpsms_with_fdr)
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

        private static List<NewPsmWithFDR> DoFalseDiscoveryRateAnalysis(IEnumerable<PSMwithTargetDecoyKnown> items)
        {
            List<NewPsmWithFDR> ids = new List<NewPsmWithFDR>();

            int cumulative_target = 0;
            int cumulative_decoy = 0;
            foreach (PSMwithTargetDecoyKnown item in items)
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

        private void AddObservedPeptidesToDictionary()
        {
            status("Adding new keys to peptide dictionary...");
            foreach (var ah in newPsms)
            {
                if (ah != null)
                    foreach (var fhh in ah)
                    {
                        if (fhh != null)
                        {
                            var cp = fhh.GetCompactPeptide(variableModifications, localizeableModifications);
                            if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                                compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                        }
                    }
            }

            int totalProteins = proteinList.Count;

            status("Adding possible sources to peptide dictionary...");

            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());

                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length == 1 || peptide.Length > 252)
                            continue;

                        peptide.SetFixedModifications(fixedModifications);

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, 4098, 3, localizeableModifications).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(yyy, variableModifications, localizeableModifications), out v))
                            {
                                v.Add(yyy);
                            }
                        }
                    }
                }
                lock (compactPeptideToProteinPeptideMatching)
                {
                    foreach (var ye in local)
                    {
                        HashSet<PeptideWithSetModifications> v;
                        if (compactPeptideToProteinPeptideMatching.TryGetValue(ye.Key, out v))
                        {
                            foreach (var huh in ye.Value)
                                v.Add(huh);
                        }
                    }
                }
            });
        }

        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> ApplyProteinParsimony(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // makes dictionary with proteins as keys and list of associated peptides as the value (swaps input parameter dictionary keys/values)
            Dictionary<Protein, HashSet<CompactPeptide>> newDict = new Dictionary<Protein, HashSet<CompactPeptide>>();

            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                foreach (var virtualPeptide in kvp.Value)
                {
                    // find peptides associated with new protein
                    if (!newDict.ContainsKey(virtualPeptide.protein))
                    {
                        HashSet<CompactPeptide> peptides = new HashSet<CompactPeptide>();

                        foreach (var kvp1 in fullSequenceToProteinPeptideMatching)
                        {
                            foreach(var virtualPeptide2 in kvp1.Value)
                            {
                                // if the virtual peptides come from the same protein, add the peptide to the protein's list
                                if (virtualPeptide.protein == virtualPeptide2.protein)
                                {
                                    peptides.Add(kvp1.Key);
                                }
                            }
                        }

                        newDict.Add(virtualPeptide.protein, peptides);
                    }
                }
            }



            // greedy algorithm adds the next protein that will account for the most unaccounted-for peptides
            Dictionary<Protein, HashSet<CompactPeptide>> parsimonyDict = new Dictionary<Protein, HashSet<CompactPeptide>>();
            HashSet<CompactPeptide> usedPeptides = new HashSet<CompactPeptide>();
            HashSet<CompactPeptide> bestProteinPeptideList = new HashSet<CompactPeptide>();
            Protein bestProtein = null;
            int currentBestNumNewPeptides = -1;
            int j = 0;

            // as long as there are peptides that have not been accounted for, keep going
            while (currentBestNumNewPeptides != 0)
            {
                currentBestNumNewPeptides = 0;

                // attempt to find protein that best accounts for unaccounted-for peptides
                foreach (var kvp in newDict)
                {
                    int comparisonProteinNewPeptides = 0;

                    // determines number of unaccounted-for peptides for the current protein
                    foreach (CompactPeptide peptide in kvp.Value)
                    {
                        if (!usedPeptides.Contains(peptide))
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
                if (currentBestNumNewPeptides > 0)
                {
                    if (j < 10)
                    {
                        Console.WriteLine(bestProtein.name + " is the best remaining protein associated with " + bestProteinPeptideList.Count + " peptides, adding it");
                        j++;
                    }

                    parsimonyDict.Add(bestProtein, bestProteinPeptideList);

                    // the new best protein's peptides have been accounted for
                    foreach (CompactPeptide peptide in bestProteinPeptideList)
                    {
                        usedPeptides.Add(peptide);
                    }
                }
            }

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> answer = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!answer.ContainsKey(peptide))
                    {
                        HashSet<PeptideWithSetModifications> virtualPeptides = new HashSet<PeptideWithSetModifications>();

                        foreach (var kvp1 in fullSequenceToProteinPeptideMatching)
                        {
                            foreach(var virtualPeptide in kvp1.Value)
                            {
                                if(virtualPeptide.protein == kvp.Key)
                                {
                                    virtualPeptides.Add(virtualPeptide);
                                }
                            }
                        }

                        answer.Add(peptide, virtualPeptides);
                    }
                }
            }

            return answer;
        }

        /*
        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> ApplyProteinParsimony(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // makes dictionary with proteins as keys and list of associated peptides as the value (swaps input parameter dictionary keys/values)
            Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>> newDict = new Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>>();
            

            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                foreach (var virtualPeptide in kvp.Value)
                {
                    if (!newDict.ContainsKey(virtualPeptide))
                    {
                        HashSet<CompactPeptide> peptides = new HashSet<CompactPeptide>();

                        foreach (var kvp1 in fullSequenceToProteinPeptideMatching)
                        {
                            if (kvp1.Value.Contains(virtualPeptide))
                            {
                                peptides.Add(kvp1.Key);
                            }
                        }
                        newDict.Add(virtualPeptide, peptides);
                    }
                }
            }

            // greedy algorithm adds the next protein that will account for the most unaccounted-for peptides
            Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>> parsimonyDict = new Dictionary<PeptideWithSetModifications, HashSet<CompactPeptide>>();
            HashSet<CompactPeptide> usedPeptides = new HashSet<CompactPeptide>();
            HashSet<CompactPeptide> bestProteinPeptideList = new HashSet<CompactPeptide>();
            PeptideWithSetModifications bestProtein = null;
            int currentBestNumNewPeptides = -1;
            int j = 0;
            // as long as there are peptides that have not been accounted for, keep going
            while (currentBestNumNewPeptides != 0)
            {
                currentBestNumNewPeptides = 0;

                // attempt to find protein that best accounts for unaccounted-for peptides
                foreach (var kvp in newDict)
                {
                    int comparisonProteinNewPeptides = 0;

                    // determines number of unaccounted-for peptides for the current protein
                    foreach (CompactPeptide peptide in kvp.Value)
                    {
                        if (!usedPeptides.Contains(peptide))
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
                if (currentBestNumNewPeptides > 0)
                {
                    if (j < 10)
                    {
                        Console.WriteLine(bestProtein.BaseSequence + " is the best remaining protein associated with " + bestProteinPeptideList.Count + " peptides, adding it");
                        j++;
                    }

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
        */

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

        public override void ValidateParams()
        {
            if (newPsms == null)
                throw new EngineValidationException("newPsms cannot be null");
        }
    }
}