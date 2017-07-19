﻿using Proteomics;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private const int max_mods_for_peptide = 3;
        private readonly List<PsmParent>[] allPsms;
        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;
        private readonly List<Protein> proteinList;
        private readonly List<MassDiffAcceptor> massDiffAcceptors;
        private readonly Protease protease;
        private readonly int maxMissedCleavages;
        private readonly int? minPeptideLength;
        private readonly int? maxPeptideLength;
        private readonly InitiatorMethionineBehavior initiatorMethionineBehavior;
        private readonly List<ModificationWithMass> fixedModifications;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly int maxModificationIsoforms;

        #endregion Private Fields

        #region Public Constructors

        public SequencesToActualProteinPeptidesEngine(List<PsmParent>[] allPsms, Dictionary<ModificationWithMass, ushort> modsDictionary, List<Protein> proteinList, List<MassDiffAcceptor> massDiffAcceptors, Protease protease, int maxMissedCleavages, int? minPeptideLength, int? maxPeptideLength, InitiatorMethionineBehavior initiatorMethionineBehavior, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> variableModifications, int maxModificationIsoforms, List<string> nestedIds) : base(nestedIds)
        {
            this.proteinList = proteinList;
            this.massDiffAcceptors = massDiffAcceptors;
            this.modsDictionary = modsDictionary;
            this.allPsms = allPsms;
            this.protease = protease;
            this.maxMissedCleavages = maxMissedCleavages;
            this.minPeptideLength = minPeptideLength;
            this.maxPeptideLength = maxPeptideLength;
            this.initiatorMethionineBehavior = initiatorMethionineBehavior;
            this.fixedModifications = fixedModifications;
            this.variableModifications = variableModifications;
            this.maxModificationIsoforms = maxModificationIsoforms;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            #region Match Seqeunces to PeptideWithSetModifications

            //myAnalysisResults.AddText("Starting compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            //Status("Adding observed peptides to dictionary...", new List<string> { taskId });
            foreach (var psmListForAspecificSerchMode in allPsms)
                if (psmListForAspecificSerchMode != null)
                    foreach (var psm in psmListForAspecificSerchMode)
                        if (psm != null)
                        {
                            var cp = psm.GetCompactPeptide(modsDictionary);
                            if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                                compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                        }
            //myAnalysisResults.AddText("Ending compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            int totalProteins = proteinList.Count;
            int proteinsSeen = 0;
            int old_progress = 0;
            var obj = new object();
            //Status("Adding possible sources to peptide dictionary...", new List<string> { taskId });
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                for (int i = fff.Item1; i < fff.Item2; i++)
                    foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maxMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications))
                    {
                        if (peptideWithPossibleModifications.Length <= 1)
                            continue;
                        foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModificationIsoforms, max_mods_for_peptide))
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications, modsDictionary), out v))
                                v.Add(peptideWithSetModifications);
                        }
                    }
                lock (obj)
                {
                    foreach (var ye in local)
                    {
                        HashSet<PeptideWithSetModifications> v;
                        if (compactPeptideToProteinPeptideMatching.TryGetValue(ye.Key, out v))
                            foreach (var huh in ye.Value)
                                v.Add(huh);
                    }
                    proteinsSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        //ReportProgress(new ProgressEventArgs(new_progress, "In adding possible sources to peptide dictionary loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });

            #endregion Match Seqeunces to PeptideWithSetModifications

            //Status("Computing info about actual peptides with modifications...", new List<string> { taskId });
            for (int j = 0; j < massDiffAcceptors.Count; j++)
            {
                if (allPsms[j] != null)
                {
                    foreach (var huh in allPsms[j])
                    {
                        if (huh != null && huh.Pli == null)
                            huh.SetProteinLinkedInfo(compactPeptideToProteinPeptideMatching, modsDictionary);
                    }
                }
            }
            return new SequencesToActualProteinPeptidesEngineResults(this, compactPeptideToProteinPeptideMatching);
        }

        #endregion Protected Methods

    }
}