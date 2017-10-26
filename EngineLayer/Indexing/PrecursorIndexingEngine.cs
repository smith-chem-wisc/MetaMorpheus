using Proteomics;
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
        #region Public Constructors

        public PrecursorIndexingEngine(List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ProductType> lp, int currentPartition, DecoyType decoyType, IEnumerable<DigestionParams> CollectionOfDigestionParams, CommonParameters commonParams, double maxFragmentSize, List<string> nestedIds) : base(proteinList, variableModifications, fixedModifications, lp, currentPartition, decoyType, CollectionOfDigestionParams, commonParams, maxFragmentSize, nestedIds)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("Precursor Mass Only");
            sb.AppendLine("Index partitions: " + currentPartition + "/" + commonParams.TotalPartitions);
            sb.AppendLine("Search Decoys: " + decoyType);
            sb.AppendLine("Number of proteins: " + proteinList.Count);
            sb.AppendLine("Number of fixed mods: " + fixedModifications.Count);
            sb.AppendLine("Number of variable mods: " + variableModifications.Count);
            sb.AppendLine("lp: " + string.Join(",", lp));
            foreach (var digestionParams in CollectionOfDigestionParams)
            {
                sb.AppendLine("protease: " + digestionParams.Protease);
                sb.AppendLine("initiatorMethionineBehavior: " + digestionParams.InitiatorMethionineBehavior);
                sb.AppendLine("maximumMissedCleavages: " + digestionParams.MaxMissedCleavages);
                sb.AppendLine("minPeptideLength: " + digestionParams.MinPeptideLength);
                sb.AppendLine("maxPeptideLength: " + digestionParams.MaxPeptideLength);
                sb.AppendLine("maximumVariableModificationIsoforms: " + digestionParams.MaxModificationIsoforms);
            }
            sb.Append("Localizeable mods: " + proteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Sum());
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(lp);

            // digest database
            HashSet<CompactPeptide> peptideToId = new HashSet<CompactPeptide>();

            Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = commonParams.MaxThreadsToUsePerFile }, range =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    foreach (var digestionParams in CollectionOfDigestionParams)
                    {
                        var digestedList = proteinList[i].Digest(digestionParams, fixedModifications, variableModifications).ToList();
                        foreach (var pepWithSetMods in digestedList)
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
                    var percentProgress = (int)((progress / proteinList.Count) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins for precursor...", nestedIds));
                    }
                }
            });

            // sort peptides by mass
            var peptidesSortedByMass = peptideToId.AsParallel().WithDegreeOfParallelism(commonParams.MaxThreadsToUsePerFile).OrderBy(p => p.MonoisotopicMassIncludingFixedMods).ToList();
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

            var fragmentIndex = new List<int>[maxFragmentMass * fragmentBinsPerDalton+1];

            // populate fragment index
            progress = 0;
            oldPercentProgress = 0;
            for (int i = 0; i < peptidesSortedByMass.Count; i++)
            {
                double mz = Chemistry.ClassExtensions.ToMz(peptidesSortedByMass[i].MonoisotopicMassIncludingFixedMods, 1);

                int fragmentBin = (int)Math.Round(mz * fragmentBinsPerDalton);

                if (fragmentIndex[fragmentBin] == null)
                    fragmentIndex[fragmentBin] = new List<int> { i };
                else
                    fragmentIndex[fragmentBin].Add(i);

                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Creating fragment index for precursor...", nestedIds));
                }
            }

            return new IndexingResults(peptidesSortedByMass, fragmentIndex, this);
        }

        #endregion Protected Methods
    }
}