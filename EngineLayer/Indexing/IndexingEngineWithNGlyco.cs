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
    public class IndexingEngineWithNGlyco : IndexingEngine
    {
        public IndexingEngineWithNGlyco(List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ProductType> productTypes, int currentPartition, DecoyType decoyType, IEnumerable<DigestionParams> CollectionOfDigestionParams, CommonParameters commonParams, double maxFragmentSize, List<string> nestedIds) : base(proteinList, variableModifications, fixedModifications, productTypes, currentPartition, decoyType, CollectionOfDigestionParams, commonParams, maxFragmentSize, nestedIds)
        {

        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Partitions: " + CurrentPartition + "/" + commonParameters.TotalPartitions);
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

            Parallel.ForEach(Partitioner.Create(0, ProteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (range, loopState) =>
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
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins...", nestedIds));
                    }
                }
            });

            // sort peptides by mass
            var peptidesSortedByMass = peptideToId.AsParallel().WithDegreeOfParallelism(commonParameters.MaxThreadsToUsePerFile).OrderBy(p => p.MonoisotopicMassIncludingFixedMods).ToList();
            peptideToId = null;

            // create fragment index
            List<int>[] fragmentIndex;
            List<int>[] fragmentIndexNGly;

            try
            {
                fragmentIndex = new List<int>[(int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1];
                fragmentIndexNGly = new List<int>[(int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1];
            }
            catch (OutOfMemoryException)
            {
                throw new MetaMorpheusException("Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }

            // populate fragment index
            progress = 0;
            oldPercentProgress = 0;
            for (int peptideId = 0; peptideId < peptidesSortedByMass.Count; peptideId++)
            {
                var validFragments = peptidesSortedByMass[peptideId].ProductMassesMightHaveDuplicatesAndNaNs(ProductTypes).Distinct().Where(p => !Double.IsNaN(p));

                foreach (var theoreticalFragmentMass in validFragments)
                {
                    if (theoreticalFragmentMass < MaxFragmentSize && theoreticalFragmentMass > 0)
                    {
                        int fragmentBin = (int)Math.Round(theoreticalFragmentMass * FragmentBinsPerDalton);

                        if (fragmentIndex[fragmentBin] == null)
                            fragmentIndex[fragmentBin] = new List<int> { peptideId };
                        else
                            fragmentIndex[fragmentBin].Add(peptideId);
                    }
                }

                //TO DO: generate validFragmentsNGly
                var validFragmentsNGly = peptidesSortedByMass[peptideId].ProductMassesMightHaveDuplicatesAndNaNs(ProductTypes).Distinct().Where(p => !Double.IsNaN(p));
                foreach (var theoreticalFragmentMass in validFragmentsNGly)
                {
                    if (theoreticalFragmentMass < MaxFragmentSize && theoreticalFragmentMass > 0)
                    {
                        int fragmentBin = (int)Math.Round(theoreticalFragmentMass * FragmentBinsPerDalton);

                        if (fragmentIndexNGly[fragmentBin] == null)
                            fragmentIndexNGly[fragmentBin] = new List<int> { peptideId };
                        else
                            fragmentIndexNGly[fragmentBin].Add(peptideId);
                    }
                }

                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Creating fragment index...", nestedIds));
                }
            }

            return new IndexingResultsWithNGlyco(peptidesSortedByMass, fragmentIndex, fragmentIndexNGly, this);
        }
    }
}