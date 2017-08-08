using Proteomics;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngine : MetaMorpheusEngine
    {
        #region Protected Fields

        protected const int max_mods_for_peptide = 3;
        protected readonly int maxMissedCleavages;
        protected readonly int? minPeptideLength;
        protected readonly int? maxPeptideLength;
        protected readonly InitiatorMethionineBehavior initiatorMethionineBehavior;
        protected readonly List<ModificationWithMass> fixedModifications;
        protected readonly List<ModificationWithMass> variableModifications;
        protected readonly int maxModificationIsoforms;
        protected readonly List<Psm>[] allPsms;
        protected readonly List<Protein> proteinList;
        protected readonly List<MassDiffAcceptor> massDiffAcceptors;
        protected readonly Protease protease;
        protected readonly List<ProductType> lp;

        #endregion Protected Fields

        #region Public Constructors

        public SequencesToActualProteinPeptidesEngine(List<Psm>[] allPsms, List<Protein> proteinList, List<MassDiffAcceptor> massDiffAcceptors, Protease protease, int maxMissedCleavages, int? minPeptideLength, int? maxPeptideLength, InitiatorMethionineBehavior initiatorMethionineBehavior, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> variableModifications, int maxModificationIsoforms, List<string> nestedIds, List<ProductType> lp) : base(nestedIds)
        {
            this.proteinList = proteinList;
            this.massDiffAcceptors = massDiffAcceptors;
            this.allPsms = allPsms;
            this.protease = protease;
            this.maxMissedCleavages = maxMissedCleavages;
            this.minPeptideLength = minPeptideLength;
            this.maxPeptideLength = maxPeptideLength;
            this.initiatorMethionineBehavior = initiatorMethionineBehavior;
            this.fixedModifications = fixedModifications;
            this.variableModifications = variableModifications;
            this.maxModificationIsoforms = maxModificationIsoforms;
            this.lp = lp;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();

            #region Match Sequences to PeptideWithSetModifications

            //myAnalysisResults.AddText("Starting compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            //Status("Adding observed peptides to dictionary...", new List<string> { taskId });
            foreach (var psmListForAspecificSerchMode in allPsms)
                if (psmListForAspecificSerchMode != null)
                    foreach (var psm in psmListForAspecificSerchMode)
                        if (psm != null)
                        {
                            foreach (var cp in psm.CompactPeptides)
                                if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp.Key))
                                    compactPeptideToProteinPeptideMatching.Add(cp.Key, new HashSet<PeptideWithSetModifications>());
                        }
            //myAnalysisResults.AddText("Ending compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            int totalProteins = proteinList.Count;
            int proteinsSeen = 0;
            int old_progress = 0;
            var obj = new object();
            //Status("Adding possible sources to peptide dictionary...", new List<string> { taskId });
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                for (int i = fff.Item1; i < fff.Item2; i++)
                    foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maxMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications))
                    {
                        foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModificationIsoforms, max_mods_for_peptide))
                        {
                            if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications), out HashSet<PeptideWithSetModifications> v))

                                v.Add(peptideWithSetModifications);
                        }
                    }
                lock (obj)
                {
                    foreach (var ye in local)
                    {
                        if (compactPeptideToProteinPeptideMatching.TryGetValue(ye.Key, out HashSet<PeptideWithSetModifications> v))
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

            return new SequencesToActualProteinPeptidesEngineResults(this, compactPeptideToProteinPeptideMatching);
        }

        #endregion Protected Methods
    }
}