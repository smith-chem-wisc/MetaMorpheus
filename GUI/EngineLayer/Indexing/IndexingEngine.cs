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
    public class IndexingEngine : MetaMorpheusEngine
    {
        #region Protected Fields

        protected const int fragmentBinsPerDalton = 1000;
        protected readonly List<Protein> proteinList;

        protected readonly List<ModificationWithMass> fixedModifications;
        protected readonly List<ModificationWithMass> variableModifications;
        protected readonly List<ProductType> lp;
        protected readonly int currentPartition;
        protected readonly DecoyType decoyType;
        protected readonly IEnumerable<DigestionParams> CollectionOfDigestionParams;
        protected readonly CommonParameters commonParams;
        protected readonly double maxFragmentSize;

        #endregion Protected Fields

        #region Public Constructors

        public IndexingEngine(List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ProductType> lp, int currentPartition, DecoyType decoyType, IEnumerable<DigestionParams> CollectionOfDigestionParams, CommonParameters commonParams, double maxFragmentSize, List<string> nestedIds) : base(nestedIds)
        {
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.lp = lp;
            this.currentPartition = currentPartition + 1;
            this.decoyType = decoyType;
            this.CollectionOfDigestionParams = CollectionOfDigestionParams;
            this.commonParams = commonParams;
            this.maxFragmentSize = maxFragmentSize;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Partitions: " + currentPartition + "/" + commonParams.TotalPartitions);
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
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins...", nestedIds));
                    }
                }
            });

            // sort peptides by mass
            var peptidesSortedByMass = peptideToId.AsParallel().WithDegreeOfParallelism(commonParams.MaxThreadsToUsePerFile).OrderBy(p => p.MonoisotopicMassIncludingFixedMods).ToList();
            peptideToId = null;

            // create fragment index
            List<int>[] fragmentIndex;

            try
            {
                fragmentIndex = new List<int>[(int)Math.Ceiling(maxFragmentSize) * fragmentBinsPerDalton + 1];
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
                var validFragments = peptidesSortedByMass[peptideId].ProductMassesMightHaveDuplicatesAndNaNs(lp).Distinct().Where(p => !Double.IsNaN(p));

                foreach (var theoreticalFragmentMass in validFragments)
                {
                    if (theoreticalFragmentMass < maxFragmentSize)
                    {
                        int fragmentBin = (int)Math.Round(theoreticalFragmentMass * fragmentBinsPerDalton);

                        if (fragmentIndex[fragmentBin] == null)
                            fragmentIndex[fragmentBin] = new List<int> { peptideId };
                        else
                            fragmentIndex[fragmentBin].Add(peptideId);
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

            return new IndexingResults(peptidesSortedByMass, fragmentIndex, this);
        }

        #endregion Protected Methods
    }
}