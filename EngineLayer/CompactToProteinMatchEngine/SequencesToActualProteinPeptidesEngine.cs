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

        protected readonly List<ModificationWithMass> fixedModifications;
        protected readonly List<ModificationWithMass> variableModifications;
        protected readonly List<Psm> allPsms;
        protected readonly List<Protein> proteinList;
        protected readonly TerminusType terminusType;
        protected readonly IEnumerable<DigestionParams> collectionOfDigestionParams;
        protected readonly bool reportAllAmbiguity;

        #endregion Protected Fields

        #region Public Constructors

        public SequencesToActualProteinPeptidesEngine(List<Psm> allPsms, List<Protein> proteinList, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> variableModifications, List<ProductType> ionTypes, IEnumerable<DigestionParams> collectionOfDigestionParams, bool reportAllAmbiguity, List<string> nestedIds) : base(nestedIds)
        {
            this.proteinList = proteinList;
            this.allPsms = allPsms;
            this.fixedModifications = fixedModifications;
            this.variableModifications = variableModifications;
            this.terminusType = ProductTypeMethod.IdentifyTerminusType(ionTypes);
            this.collectionOfDigestionParams = collectionOfDigestionParams;
            this.reportAllAmbiguity = reportAllAmbiguity;
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

            foreach (var psm in allPsms)
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
                    foreach (var digestionParam in collectionOfDigestionParams)
                        foreach (var peptideWithPossibleModifications in proteinList[i].Digest(digestionParam, fixedModifications))
                        {
                            foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(digestionParam, variableModifications))
                            {
                                if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications, terminusType), out HashSet<PeptideWithSetModifications> v))
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

            #endregion Match Sequences to PeptideWithSetModifications

            if (!reportAllAmbiguity)
                ResolveAmbiguities(compactPeptideToProteinPeptideMatching);

            return new SequencesToActualProteinPeptidesEngineResults(this, compactPeptideToProteinPeptideMatching);
        }

        protected static void ResolveAmbiguities(Dictionary<CompactPeptideBase,HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching)
        {
            //If ambiguities are not desired, a single compact peptide has survived to this point for each PSM.
            //This single peptide sequence can originate from multiple unique proteins
            //When there are multiple origins, both target and decoy peptides (even though they have the same sequence) can arise as of 9/15/17
            //There should not be any overlap between the target and decoy databases, as TDA assumes that the decoy sequence is incorrect
            //To prevent this error from dramatically overestimating the FDR, a violation of the 50-50 target-decoy ratio will be used for the time being, and decoys that share the same sequence as targets will be ignored
            List<CompactPeptideBase> keys = compactPeptideToProteinPeptideMatching.Keys.ToList();
            foreach (CompactPeptide key in keys)
            {
                HashSet<PeptideWithSetModifications> value = compactPeptideToProteinPeptideMatching[key];
                compactPeptideToProteinPeptideMatching[key] = new HashSet<PeptideWithSetModifications> { value.FirstOrDefault(b => !b.Protein.IsDecoy) ?? value.First() };
            }
        }

        #endregion Protected Methods
    }
}