using MassSpectrometry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngine : MetaMorpheusEngine
    {
        protected readonly List<ModificationWithMass> FixedModifications;
        protected readonly List<ModificationWithMass> VariableModifications;
        protected readonly List<PeptideSpectralMatch> AllPsms;
        protected readonly List<Protein> Proteins;
        protected readonly TerminusType TerminusType;
        protected readonly IEnumerable<DigestionParams> CollectionOfDigestionParams;
        protected readonly bool ReportAllAmbiguity;

        public SequencesToActualProteinPeptidesEngine(List<PeptideSpectralMatch> allPsms, List<Protein> proteinList,
                List<ModificationWithMass> fixedModifications, List<ModificationWithMass> variableModifications,
                List<ProductType> ionTypes, IEnumerable<DigestionParams> collectionOfDigestionParams,
                bool reportAllAmbiguity, CommonParameters commonParameters, List<string> nestedIds)
            : base(commonParameters, nestedIds)
        {
            Proteins = proteinList;
            AllPsms = allPsms;
            FixedModifications = fixedModifications;
            VariableModifications = variableModifications;
            TerminusType = ProductTypeMethods.IdentifyTerminusType(ionTypes);
            CollectionOfDigestionParams = collectionOfDigestionParams;
            ReportAllAmbiguity = reportAllAmbiguity;
        }

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
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();

            // Match Sequences to PeptideWithSetModifications
            foreach (var psm in AllPsms)
            {
                if (psm != null)
                {
                    foreach (var compactPeptide in psm.CompactPeptides)
                    {
                        if (!compactPeptideToProteinPeptideMatching.ContainsKey(compactPeptide.Key))
                        {
                            compactPeptideToProteinPeptideMatching.Add(compactPeptide.Key, new HashSet<PeptideWithSetModifications>());
                        }
                    }
                }
            }

            double proteinsMatched = 0;
            int oldPercentProgress = 0;

            Parallel.ForEach(Partitioner.Create(0, Proteins.Count),
                new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
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

                    foreach (var digestionParam in CollectionOfDigestionParams.ToList())
                    {
                        // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                        foreach (var peptide in Proteins[i].Digest(digestionParam, FixedModifications, VariableModifications))
                        {
                            var compactPeptide = peptide.CompactPeptide(TerminusType);

                            if (compactPeptideToProteinPeptideMatching.TryGetValue(compactPeptide, out var peptidesWithSetMods))
                            {
                                lock (peptidesWithSetMods)
                                {
                                    peptidesWithSetMods.Add(peptide);
                                }
                            }
                        }
                    }
                }

                // report progress (proteins matched so far out of total proteins in database)
                proteinsMatched++;
                var percentProgress = (int)((proteinsMatched / Proteins.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Matching peptides to proteins... ", NestedIds));
                }
            });

            if (!ReportAllAmbiguity)
            {
                ResolveAmbiguities(compactPeptideToProteinPeptideMatching);
            }

            return new SequencesToActualProteinPeptidesEngineResults(this, compactPeptideToProteinPeptideMatching);
        }
    }
}