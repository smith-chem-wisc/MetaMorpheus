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
        protected readonly List<PeptideSpectralMatch> allPsms;
        protected readonly List<Protein> proteins;
        protected readonly TerminusType terminusType;
        protected readonly IEnumerable<DigestionParams> collectionOfDigestionParams;
        protected readonly bool reportAllAmbiguity;

        #endregion Protected Fields

        #region Public Constructors

        //has input for collection of all digestion parameters--- can we split here for compact to protein matching
        public SequencesToActualProteinPeptidesEngine(List<PeptideSpectralMatch> allPsms, List<Protein> proteinList, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> variableModifications, List<ProductType> ionTypes, IEnumerable<DigestionParams> collectionOfDigestionParams, bool reportAllAmbiguity, List<string> nestedIds) : base(nestedIds)
        {
            this.proteins = proteinList;
            this.allPsms = allPsms;
            this.fixedModifications = fixedModifications;
            this.variableModifications = variableModifications;
            this.terminusType = ProductTypeMethod.IdentifyTerminusType(ionTypes);
            this.collectionOfDigestionParams = collectionOfDigestionParams;
            this.reportAllAmbiguity = reportAllAmbiguity;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected static void ResolveAmbiguities(Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching)
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

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // split by protease to group them together in a dictionary with key being the protease and the value being the Dictionary <CompactPeptideBase, HashSet<PeptideWithSetModification>>
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            HashSet<Protease> proteases = new HashSet<Protease>();
               
            foreach (DigestionParams dp in collectionOfDigestionParams)
            {
                if (!proteases.Contains(dp.Protease))
                {
                    proteases.Add(dp.Protease);
                }
                
            }

            Dictionary<Protease, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>> proteaseSortedCompactPeptideToProteinPeptideMatching = new Dictionary<Protease, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>>();

            foreach (Protease p in proteases)
            {
                Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
                List<PeptideSpectralMatch> proteasePsms = new List<PeptideSpectralMatch>();
                foreach (var psm in allPsms)
                {
                    if (psm.DigestionParams.Protease == p && psm != null)
                    {
                        proteasePsms.Add(psm);
                    }
                }
                foreach (var psm in proteasePsms)
                {
                    foreach (var compactPeptide in psm.CompactPeptides)
                    {
                        if (!compactPeptideToProteinPeptideMatching.ContainsKey(compactPeptide.Key))
                            compactPeptideToProteinPeptideMatching.Add(compactPeptide.Key, new HashSet<PeptideWithSetModifications>());
                    }
                }
               

                Parallel.ForEach(Partitioner.Create(0, proteins.Count), fff =>
                {
                    for (int i = fff.Item1; i < fff.Item2; i++)
                    {
                        foreach (var digestionParam in collectionOfDigestionParams)
                        {
                            foreach (var peptide in proteins[i].Digest(digestionParam, fixedModifications, variableModifications))
                            {
                                var compactPeptide = peptide.CompactPeptide(terminusType);

                                if (compactPeptideToProteinPeptideMatching.TryGetValue(compactPeptide, out var peptidesWithSetMods))
                                {
                                    lock (peptidesWithSetMods)
                                        peptidesWithSetMods.Add(peptide);
                                }
                            }
                        }
                    }
                });
                proteaseSortedCompactPeptideToProteinPeptideMatching.Add(p, compactPeptideToProteinPeptideMatching);
                if (!reportAllAmbiguity)
                    ResolveAmbiguities(compactPeptideToProteinPeptideMatching);
            }

            return new SequencesToActualProteinPeptidesEngineResults(this, proteaseSortedCompactPeptideToProteinPeptideMatching);
        }

        #endregion Protected Methods
    }
}