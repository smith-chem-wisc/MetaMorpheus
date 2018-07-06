using EngineLayer.CrosslinkSearch;
using MassSpectrometry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkAnalysis
{
    public class CrosslinkAnalysisEngine : MetaMorpheusEngine
    {
        protected readonly TerminusType TerminusType;

        private readonly List<PsmCross> NewPsms;
        private readonly List<Protein> ProteinList;
        private readonly List<ModificationWithMass> VariableModifications;
        private readonly List<ModificationWithMass> FixedModifications;

        private readonly List<ProductType> ProductTypes;
        private readonly CrosslinkerTypeClass Crosslinker;

        private readonly Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> CompactPeptideToProteinPeptideMatching;
        private readonly string OutputFolder;

        public CrosslinkAnalysisEngine(List<PsmCross> newPsms, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ProductType> lp, string outputFolder, CrosslinkerTypeClass crosslinker, TerminusType terminusType, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            NewPsms = newPsms;
            CompactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            ProteinList = proteinList;
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            ProductTypes = lp;
            OutputFolder = outputFolder;
            Crosslinker = crosslinker;
            TerminusType = terminusType;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            CrosslinkAnalysisResults myAnalysisResults = new CrosslinkAnalysisResults(this);
            Status("Running analysis engine!");
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy

            // Match Seqeunces to PeptideWithSetModifications

            //myAnalysisResults.AddText("Starting compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            Status("Adding observed peptides to dictionary...");
            foreach (var psmpair in NewPsms)
            {
                if (psmpair != null)
                {
                    var cp = psmpair.compactPeptide;
                    if (!CompactPeptideToProteinPeptideMatching.ContainsKey(cp))
                        CompactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());

                    if (psmpair.BetaPsmCross != null)
                    {
                        var cp1 = psmpair.BetaPsmCross.compactPeptide;
                        if (!CompactPeptideToProteinPeptideMatching.ContainsKey(cp1))
                            CompactPeptideToProteinPeptideMatching.Add(cp1, new HashSet<PeptideWithSetModifications>());
                    }
                }
            }

            int proteinsSeen = 0;
            int old_progress = 0;
            Status("Adding possible sources to peptide dictionary...");
            Parallel.ForEach(Partitioner.Create(0, ProteinList.Count),
                new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (fff, loopState) =>
            {
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    foreach (var peptideWithSetModifications in ProteinList[i].Digest(commonParameters.DigestionParams, FixedModifications, VariableModifications))
                    {
                        var compactPeptide = peptideWithSetModifications.CompactPeptide(TerminusType);
                        if (CompactPeptideToProteinPeptideMatching.TryGetValue(compactPeptide, out HashSet<PeptideWithSetModifications> peptidesWithSetMods))
                        {
                            lock (peptidesWithSetMods)
                                peptidesWithSetMods.Add(peptideWithSetModifications);
                        }
                    }
                    proteinsSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)proteinsSeen / (ProteinList.Count) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In adding possible" +
                            " sources to peptide dictionary loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });

            // done with Match Seqeunces to PeptideWithSetModifications

            Status("Computing info about actual peptides with modifications...");
            for (int myScanWithMassIndex = 0; myScanWithMassIndex < NewPsms.Count; myScanWithMassIndex++)
            {
                if (NewPsms[myScanWithMassIndex] != null)
                {
                    NewPsms[myScanWithMassIndex].MatchToProteinLinkedPeptides(CompactPeptideToProteinPeptideMatching);
                    if (NewPsms[myScanWithMassIndex].BetaPsmCross != null)
                    {
                        NewPsms[myScanWithMassIndex].BetaPsmCross.MatchToProteinLinkedPeptides(CompactPeptideToProteinPeptideMatching);
                    }
                }
            }

            return myAnalysisResults;
        }
    }
}