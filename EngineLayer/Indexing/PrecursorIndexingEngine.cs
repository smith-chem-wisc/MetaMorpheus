using MassSpectrometry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace EngineLayer.Indexing
{
    public class PrecursorIndexingEngine : IndexingEngine
    {
        public PrecursorIndexingEngine(List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications,
                List<ProductType> lp, int currentPartition, DecoyType decoyType, IEnumerable<DigestionParams> CollectionOfDigestionParams, CommonParameters commonParams,
                double maxFragmentSize, List<string> nestedIds)
            : base(proteinList, variableModifications, fixedModifications, lp, currentPartition, decoyType, CollectionOfDigestionParams, commonParams, maxFragmentSize, nestedIds)
        {
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("Precursor Mass Only");
            sb.AppendLine("Index partitions: " + CurrentPartition + "/" + CommonParameters.TotalPartitions);
            sb.AppendLine("Search Decoys: " + DecoyType);
            sb.AppendLine("Number of proteins: " + ProteinList.Count);
            sb.AppendLine("Number of fixed mods: " + FixedModifications.Count);
            sb.AppendLine("Number of variable mods: " + VariableModifications.Count);
            sb.AppendLine("lp: " + string.Join(",", ProductTypes));
            foreach (var digestionParams in CollectionOfDigestionParams)
            {
                sb.AppendLine("protease: " + digestionParams.Protease);
                sb.AppendLine("initiatorMethionineBehavior: " + digestionParams.InitiatorMethionineBehavior);
                sb.AppendLine("maximumMissedCleavages: " + digestionParams.MaxMissedCleavages);
                sb.AppendLine("minPeptideLength: " + digestionParams.MinPeptideLength);
                sb.AppendLine("maxPeptideLength: " + digestionParams.MaxPeptideLength);
                sb.AppendLine("maximumVariableModificationIsoforms: " + digestionParams.MaxModificationIsoforms);
            }
            sb.Append("Localizeable mods: " + ProteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Sum());
            return sb.ToString();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(ProductTypes);

            // digest database
            HashSet<CompactPeptide> peptideToId = new HashSet<CompactPeptide>();

            Parallel.ForEach(Partitioner.Create(0, ProteinList.Count),
                new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                (range, loopState) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    foreach (var digestionParams in CollectionOfDigestionParams)
                    {
                        foreach (var pepWithSetMods in ProteinList[i].Digest(digestionParams, FixedModifications, VariableModifications))
                        {
                            CompactPeptide compactPeptide = pepWithSetMods.CompactPeptide(terminusType);

                            var observed = peptideToId.Contains(compactPeptide);
                            if (observed)
                                continue;
                            lock (peptideToId)
                            {
                                observed = peptideToId.Contains(compactPeptide);
                                if (observed)
                                    continue;
                                peptideToId.Add(compactPeptide);
                            }
                        }
                    }

                    progress++;
                    var percentProgress = (int)((progress / ProteinList.Count) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins for precursor...", NestedIds));
                    }
                }
            });

            // sort peptides by mass
            var peptidesSortedByMass = peptideToId.AsParallel().WithDegreeOfParallelism(CommonParameters.MaxThreadsToUsePerFile).OrderBy(p => p.MonoisotopicMassIncludingFixedMods).ToList();
            peptideToId = null;

            // create fragment index
            int maxFragmentMass = 0;
            for (int i = peptidesSortedByMass.Count - 1; i >= 0; i--)
            {
                if (!Double.IsNaN(peptidesSortedByMass[i].MonoisotopicMassIncludingFixedMods))
                {
                    maxFragmentMass = (int)Math.Ceiling(Chemistry.ClassExtensions.ToMz(peptidesSortedByMass[i].MonoisotopicMassIncludingFixedMods, 1));
                    break;
                }
            }

            var fragmentIndex = new List<int>[maxFragmentMass * FragmentBinsPerDalton + 1];

            // populate fragment index
            progress = 0;
            oldPercentProgress = 0;
            for (int i = 0; i < peptidesSortedByMass.Count; i++)
            {
                double mz = Chemistry.ClassExtensions.ToMz(peptidesSortedByMass[i].MonoisotopicMassIncludingFixedMods, 1);
                if (!Double.IsNaN(mz))
                {
                    int fragmentBin = (int)Math.Round(mz * FragmentBinsPerDalton);

                    if (fragmentIndex[fragmentBin] == null)
                        fragmentIndex[fragmentBin] = new List<int> { i };
                    else
                        fragmentIndex[fragmentBin].Add(i);
                }
                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Creating fragment index for precursor...", NestedIds));
                }
            }

            return new IndexingResults(peptidesSortedByMass, fragmentIndex, this);
        }
    }
}