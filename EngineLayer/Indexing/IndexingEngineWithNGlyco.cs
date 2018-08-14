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
using Chemistry;
using EngineLayer.CrosslinkSearch;

namespace EngineLayer.Indexing
{
    public class IndexingEngineWithNGlyco : IndexingEngine
    {
        protected static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        protected static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        protected static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        protected static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        protected static readonly double hexNAcMass = 203.0793725330;
        protected static readonly double hexNAcCrossRingMass = 85.0527638520;

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
            List<int>[] fragmentIndexNGly = new List<int>[1];
            bool _searchGlyco = false;

            try
            {
                fragmentIndex = new List<int>[(int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1];
                if (_searchGlyco)
                {
                    fragmentIndexNGly = new List<int>[(int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1];
                }
                
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
                if (_searchGlyco)
                {
                    var validFragmentsNGly = GenerateBgYgFragments(peptidesSortedByMass[peptideId], ProductTypes).Distinct().Where(p => !Double.IsNaN(p));
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

        private List<double> GenerateBgYgFragments(CompactPeptide compactPeptide, List<ProductType> productTypes)
        {
            var len = compactPeptide.CTerminalMasses.Length;
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsBnoB1 = productTypes.Contains(ProductType.BnoB1ions);
            bool containsY = productTypes.Contains(ProductType.Y);

            var modPos = PsmCross.NGlyPosCal(compactPeptide);

            List<double> massesToReturn = new List<double>();

            foreach (var iPos in modPos)
            {
                if (compactPeptide.NTerminalMasses != null)
                {
                    for (int j = 0; j < compactPeptide.NTerminalMasses.Length; j++)
                    {
                        var hm = compactPeptide.NTerminalMasses[j];
                        if ((containsB || (containsBnoB1 && j > 0)) && j >= iPos)
                        {
                            //massesToReturn.Add( ClassExtensions.RoundedDouble(hm + 260).Value);
                            massesToReturn.Add(hm + hexNAcMass);
                            massesToReturn.Add(hm + hexNAcCrossRingMass);
                        }
                    }
                }
                if (compactPeptide.CTerminalMasses != null)
                {
                    for (int j = 0; j < compactPeptide.CTerminalMasses.Length; j++)
                    {
                        var hm = compactPeptide.CTerminalMasses[j];
                        if (containsY && j >= len - iPos + 2)
                        {
                            //massesToReturn.Add(ClassExtensions.RoundedDouble(hm + waterMonoisotopicMass + 260).Value);
                            massesToReturn.Add(hm + waterMonoisotopicMass + hexNAcMass);
                            massesToReturn.Add(hm + waterMonoisotopicMass + hexNAcCrossRingMass);
                        }
                    }
                }
            }

            return massesToReturn;
        }
    }
}