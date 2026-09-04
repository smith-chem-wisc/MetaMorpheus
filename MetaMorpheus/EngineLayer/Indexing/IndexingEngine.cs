using Chemistry;
using EngineLayer.NonSpecificEnzymeSearch;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using Omics;
using Omics.Fragmentation.Peptide;
using Omics.Modifications;
using UsefulProteomicsDatabases;

namespace EngineLayer.Indexing
{
    public class IndexingEngine : MetaMorpheusEngine
    {
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private const int FragmentBinsPerDalton = 1000;
        private readonly List<Protein> ProteinList;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly List<SilacLabel> SilacLabels;
        private readonly (SilacLabel StartLabel, SilacLabel EndLabel)? TurnoverLabels;
        private readonly int CurrentPartition;
        private readonly DecoyType DecoyType;
        private readonly double MaxFragmentSize;
        public readonly bool GeneratePrecursorIndex;
        public readonly List<FileInfo> ProteinDatabases;
        public readonly TargetContaminantAmbiguity TcAmbiguity;

        public IndexingEngine(List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications,
            List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel, int currentPartition, DecoyType decoyType,
            CommonParameters commonParams, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, 
            double maxFragmentSize, bool generatePrecursorIndex, List<FileInfo> proteinDatabases, TargetContaminantAmbiguity tcAmbiguity, List<string> nestedIds)
            : base(commonParams, fileSpecificParameters, nestedIds)
        {
            ProteinList = proteinList;
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            SilacLabels = silacLabels;
            if (startLabel != null || endLabel != null) //else it's null
            {
                TurnoverLabels = (startLabel, endLabel);
            }

            CurrentPartition = currentPartition + 1;
            DecoyType = decoyType;
            MaxFragmentSize = maxFragmentSize;
            GeneratePrecursorIndex = generatePrecursorIndex;
            ProteinDatabases = proteinDatabases;
            TcAmbiguity = tcAmbiguity;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Databases: " + string.Join(",", ProteinDatabases.OrderBy(p => p.Name).Select(p => p.Name + ":" + p.CreationTime)));
            sb.AppendLine("Partitions: " + CurrentPartition + "/" + CommonParameters.TotalPartitions);
            sb.AppendLine("Precursor Index: " + GeneratePrecursorIndex);
            sb.AppendLine("Search Decoys: " + DecoyType);
            sb.AppendLine("Number of proteins: " + ProteinList.Count);
            sb.AppendLine("Number of fixed mods: " + FixedModifications.Count);
            sb.AppendLine("Number of variable mods: " + VariableModifications.Count);
            sb.AppendLine("Dissociation Type: " + CommonParameters.DissociationType);
            sb.AppendLine("Contaminant Handling: " + TcAmbiguity);

            sb.AppendLine("protease: " + CommonParameters.DigestionParams.DigestionAgent);
            if (CommonParameters.DigestionParams is DigestionParams digestionParams)
                sb.AppendLine("initiatorMethionineBehavior: " + digestionParams.InitiatorMethionineBehavior);
            sb.AppendLine("maximumMissedCleavages: " + CommonParameters.DigestionParams.MaxMissedCleavages);
            sb.AppendLine("minPeptideLength: " + CommonParameters.DigestionParams.MinLength);
            sb.AppendLine("maxPeptideLength: " + CommonParameters.DigestionParams.MaxLength);
            sb.AppendLine("maximumVariableModificationIsoforms: " + CommonParameters.DigestionParams.MaxModificationIsoforms);
            sb.AppendLine("digestionTerminus: " + CommonParameters.DigestionParams.FragmentationTerminus);
            sb.AppendLine("maxModsForEachPeptide: " + CommonParameters.DigestionParams.MaxMods);
            sb.AppendLine("cleavageSpecificity: " + CommonParameters.DigestionParams.SearchModeType);
            if (CommonParameters.DigestionParams is DigestionParams digestionParam)
                sb.AppendLine("specificProtease: " + digestionParam.SpecificProtease);
            sb.AppendLine("maximumFragmentSize" + (int)Math.Round(MaxFragmentSize));

            sb.Append("Localizeable mods: " + ProteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Sum());
            return sb.ToString();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;

            // digest database
            List<PeptideWithSetModifications> peptides = new List<PeptideWithSetModifications>();

            if (CommonParameters.DigestionParams is not DigestionParams digestionParams)
                throw new MetaMorpheusException("Digestion parameters must be of type DigestionParams. Not yet implemented for Rna Digestion");
            
            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            // Each thread's peptides land in its own slot and are concatenated in thread order, so the
            // peptide index does not depend on which thread finishes first. Appending under a lock made
            // the order vary run to run, which changed how equal-mass peptides were ordered by the sort
            // below and so changed reported Delta Scores between otherwise identical runs.
            var peptidesPerThread = new List<PeptideWithSetModifications>[maxThreadsPerFile];
            Parallel.ForEach(threads, (i) =>
            {
                int threadIndex = i;
                List<PeptideWithSetModifications> localPeptides = new List<PeptideWithSetModifications>();

                for (; i < ProteinList.Count; i += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    localPeptides.AddRange(ProteinList[i].Digest(digestionParams, FixedModifications, VariableModifications, SilacLabels, TurnoverLabels));

                    progress++;
                    var percentProgress = (int)((progress / ProteinList.Count) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins...", NestedIds));
                    }
                }

                // left null on cancellation, matching the previous behaviour of discarding partial results
                peptidesPerThread[threadIndex] = localPeptides;
            });

            foreach (var threadPeptides in peptidesPerThread)
            {
                if (threadPeptides != null)
                {
                    peptides.AddRange(threadPeptides);
                }
            }

            // sort peptides by mass
            peptides = SortByMonoisotopicMass(peptides);

            //create precursor index (if specified)
            List<int>[] precursorIndex = null;
            if (GeneratePrecursorIndex)
            {
                precursorIndex = CreateNewPrecursorIndex(peptides);
            }
            bool addInteriorTerminalModsToPrecursorIndex = GeneratePrecursorIndex && CommonParameters.DigestionParams.DigestionAgent.Name.Contains("single");
            List<Modification> terminalModifications = addInteriorTerminalModsToPrecursorIndex ?
                NonSpecificEnzymeSearchEngine.GetVariableTerminalMods(CommonParameters.DigestionParams.FragmentationTerminus, VariableModifications) :
                null;

            // Create the fragment index in compressed sparse row form. Peptide ids are split into
            // contiguous blocks, one per thread; each block fragments its own peptides into its own
            // emission run, and the runs are concatenated in block order and counting-sorted into the
            // two flat arrays. Concatenating in block order reproduces the sequence a single sequential
            // pass produced, so the index does not depend on the thread count and peptide ids still
            // ascend within every bin.
            //
            // Fragmenting once and recording where each fragment landed replaces counting in one pass
            // and filling in a second: fragmentation is the expensive part of this loop, and the
            // count-then-fill shape paid for it twice. Digestion above was already parallel; this stage
            // was not, so a whole phase of indexing ran on one core.
            int binCount = (int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1;

            List<PeptideBlock> blocks = PartitionPeptidesIntoBlocks(peptides.Count);
            var fragmentEmissions = new FragmentEmissionRun[blocks.Count];
            var interiorTerminalModEmissions = addInteriorTerminalModsToPrecursorIndex ? new List<(int Bin, int PeptideId)>[blocks.Count] : null;
            var fragmentationProgress = new FragmentationProgress();

            progress = 0;
            oldPercentProgress = 0;

            Parallel.For(0, blocks.Count, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, CommonParameters.MaxThreadsToUsePerFile) }, blockIndex =>
            {
                FragmentPeptideBlock(blocks[blockIndex], blockIndex, peptides, addInteriorTerminalModsToPrecursorIndex,
                    terminalModifications, fragmentEmissions, interiorTerminalModEmissions, fragmentationProgress);
            });

            FragmentIndex fragmentIndex = FragmentIndexBuilder.Build(binCount, fragmentEmissions);

            // Replayed in block order, and within a block in ascending peptide id, which is the order the
            // sequential loop appended them in. The precursor index is still a List<int>[] here, so the
            // append order is the bin contents.
            if (addInteriorTerminalModsToPrecursorIndex)
            {
                foreach (var blockEmissions in interiorTerminalModEmissions)
                {
                    foreach ((int precursorBin, int peptideId) in blockEmissions)
                    {
                        if (precursorIndex[precursorBin] == null)
                        {
                            precursorIndex[precursorBin] = new List<int> { peptideId };
                        }
                        else
                        {
                            precursorIndex[precursorBin].Add(peptideId);
                        }
                    }
                }
            }


            return new IndexingResults(peptides, fragmentIndex, precursorIndex, this);
        }

        /// <summary>
        /// Counts fragmented peptides across the parallel blocks. Each whole percent is announced by
        /// exactly one thread -- whichever wins the compare-exchange that advances the counter -- so the
        /// progress bar neither stalls nor reports the same percent several times.
        /// </summary>
        private sealed class FragmentationProgress
        {
            private long _peptidesFragmented;
            private int _lastPercentReported;

            internal bool AdvanceOnePeptide(int totalPeptides, out int percentProgress)
            {
                long done = Interlocked.Increment(ref _peptidesFragmented);
                percentProgress = (int)(done * 100 / totalPeptides);

                int last = Volatile.Read(ref _lastPercentReported);
                return percentProgress > last && Interlocked.CompareExchange(ref _lastPercentReported, percentProgress, last) == last;
            }
        }

        private readonly struct PeptideBlock
        {
            public PeptideBlock(int start, int end)
            {
                Start = start;
                End = end;
            }

            /// <summary>First peptide id in the block.</summary>
            public int Start { get; }

            /// <summary>One past the last peptide id in the block.</summary>
            public int End { get; }
        }

        /// <summary>
        /// Splits the peptide ids into contiguous blocks, one per thread. Emissions are concatenated in
        /// block order, so the index that comes out does not depend on how many blocks there are.
        /// </summary>
        private List<PeptideBlock> PartitionPeptidesIntoBlocks(int peptideCount)
        {
            var blocks = new List<PeptideBlock>();
            if (peptideCount == 0)
            {
                return blocks;
            }

            int blockCount = Math.Max(1, Math.Min(CommonParameters.MaxThreadsToUsePerFile, peptideCount));
            int blockSize = (peptideCount + blockCount - 1) / blockCount;

            for (int start = 0; start < peptideCount; start += blockSize)
            {
                blocks.Add(new PeptideBlock(start, Math.Min(start + blockSize, peptideCount)));
            }

            return blocks;
        }

        /// <summary>
        /// Fragments one contiguous block of peptides, recording the fragment bins each peptide lands in
        /// and, when the search needs them, the interior terminal mod precursor bins -- which come out of
        /// the same fragmentation, so peptides are still fragmented exactly once.
        /// </summary>
        private void FragmentPeptideBlock(PeptideBlock block, int blockIndex, List<PeptideWithSetModifications> peptides,
            bool addInteriorTerminalModsToPrecursorIndex, List<Modification> terminalModifications,
            FragmentEmissionRun[] fragmentEmissions, List<(int Bin, int PeptideId)>[] interiorTerminalModEmissions,
            FragmentationProgress fragmentationProgress)
        {
            int blockLength = block.End - block.Start;
            var fragmentRun = new FragmentEmissionRun(block.Start, blockLength);
            var interiorEmissions = addInteriorTerminalModsToPrecursorIndex ? new List<(int Bin, int PeptideId)>() : null;

            List<Product> fragments = new List<Product>();

            for (int peptideId = block.Start; peptideId < block.End; peptideId++)
            {
                fragmentRun.BeginPeptide();

                peptides[peptideId].Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, fragments, CommonParameters.FragmentationParameters);

                foreach (var theoreticalFragment in fragments)
                {
                    if (TryGetFragmentBin(theoreticalFragment, out int fragmentBin))
                    {
                        fragmentRun.Add(fragmentBin);
                    }
                }

                //Add terminal mods if needed (do it here rather than earlier so that we don't have to fragment twice)
                if (addInteriorTerminalModsToPrecursorIndex)
                {
                    CollectInteriorTerminalModPrecursorBins(interiorEmissions, fragments, peptides[peptideId], peptideId, terminalModifications);
                }

                if (fragmentationProgress.AdvanceOnePeptide(peptides.Count, out int percentProgress))
                {
                    ReportProgress(new ProgressEventArgs(percentProgress, "Fragmenting peptides...", NestedIds));
                }
            }

            fragmentRun.Complete();

            fragmentEmissions[blockIndex] = fragmentRun;
            if (interiorTerminalModEmissions != null)
            {
                interiorTerminalModEmissions[blockIndex] = interiorEmissions;
            }
        }

        /// <summary>
        /// The bin a theoretical fragment belongs in, or false if it falls outside the index. Shared by the
        /// counting and filling passes so the two cannot disagree about which fragments are indexed — in
        /// particular about the low-resolution rounding, which only applies to LowCID and which a second
        /// hand-written copy of this arithmetic would be easy to omit.
        /// </summary>
        private bool TryGetFragmentBin(Product theoreticalFragment, out int fragmentBin)
        {
            double theoreticalFragmentMass = theoreticalFragment.NeutralMass;

            //if low res round
            if (CommonParameters.DissociationType == MassSpectrometry.DissociationType.LowCID)
            {
                theoreticalFragmentMass = Math.Round(theoreticalFragmentMass / 1.0005079, 0) * 1.0005079;
            }

            if (theoreticalFragmentMass < MaxFragmentSize && theoreticalFragmentMass > 0)
            {
                fragmentBin = (int)Math.Round(theoreticalFragmentMass * FragmentBinsPerDalton);
                return true;
            }

            fragmentBin = -1;
            return false;
        }

        /// <summary>
        /// Sorts by monoisotopic mass, reading each mass once instead of once per comparison —
        /// PeptideWithSetModifications.MonoisotopicMass rounds to 9 places on every read even though the
        /// unrounded value is cached, so a comparison delegate pays for it ~2n*log(n) times.
        /// Sorting (mass, position) pairs through the same comparison-based introsort reproduces the exact
        /// permutation the previous in-place sort produced, which matters because equal masses are
        /// pervasive: a reversed decoy has its target's residue composition and therefore its exact mass.
        /// </summary>
        internal static List<PeptideWithSetModifications> SortByMonoisotopicMass(List<PeptideWithSetModifications> peptides)
        {
            var order = new KeyValuePair<double, int>[peptides.Count];
            for (int i = 0; i < order.Length; i++)
            {
                order[i] = new KeyValuePair<double, int>(peptides[i].MonoisotopicMass, i);
            }

            Array.Sort(order, (x, y) => x.Key.CompareTo(y.Key));

            var sorted = new List<PeptideWithSetModifications>(order.Length);
            for (int i = 0; i < order.Length; i++)
            {
                sorted.Add(peptides[order[i].Value]);
            }
            return sorted;
        }

        private List<int>[] CreateNewPrecursorIndex(List<PeptideWithSetModifications> peptidesSortedByMass)
        {
            // create precursor index
            List<int>[] precursorIndex = null;
            try
            {
                precursorIndex = new List<int>[(int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1];
            }
            catch (OutOfMemoryException)
            {
                throw new MetaMorpheusException("Max precursor mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }

            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(0, "Creating precursor index...", NestedIds));

            //Add all the precursors
            for (int i = 0; i < peptidesSortedByMass.Count; i++)
            {
                double mass = peptidesSortedByMass[i].MonoisotopicMass;
                if (!double.IsNaN(mass))
                {
                    if (mass > MaxFragmentSize) //if the precursor is larger than the index allows, then stop adding precursors
                    {
                        break;
                    }

                    int precursorBin = (int)Math.Round(mass * FragmentBinsPerDalton);

                    if (precursorIndex[precursorBin] == null)
                    {
                        precursorIndex[precursorBin] = new List<int> { i };
                    }
                    else
                    {
                        precursorIndex[precursorBin].Add(i);
                    }
                }
                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Creating precursor index...", NestedIds));
                }
            }
            return precursorIndex;
        }

        //add possible protein/peptide terminal modifications that aren't on the terminal amino acids
        //The purpose is for terminal mods that are contained WITHIN the Single peptide
        /// <summary>
        /// Records the interior terminal mod precursor bins for one peptide instead of writing them into
        /// the precursor index directly. The blocks fragment in parallel and the precursor index is a
        /// shared List&lt;int&gt;[], so the writes are replayed in block order afterwards -- which is the
        /// order the sequential loop appended them in, and the bin contents are that order.
        /// </summary>
        private void CollectInteriorTerminalModPrecursorBins(List<(int Bin, int PeptideId)> emissions, List<Product> fragmentMasses, PeptideWithSetModifications peptide, int peptideId, List<Modification> variableModifications)
        {
            //Get database annotated mods
            Dictionary<int, List<Modification>> databaseAnnotatedMods = NonSpecificEnzymeSearchEngine.GetTerminalModPositions(peptide, CommonParameters.DigestionParams, variableModifications);
            foreach (KeyValuePair<int, List<Modification>> relevantDatabaseMod in databaseAnnotatedMods)
            {
                int fragmentNumber = relevantDatabaseMod.Key;
                Product fragmentAtIndex = fragmentMasses.FirstOrDefault(x => x.FragmentNumber == fragmentNumber);

                double basePrecursorMass;
                if (fragmentAtIndex is null)
                {
                    basePrecursorMass = peptide.MonoisotopicMass;
                }
                else
                {
                    basePrecursorMass = fragmentAtIndex.NeutralMass -
                                        DissociationTypeCollection.GetMassShiftFromProductType(fragmentAtIndex.ProductType) +
                                        WaterMonoisotopicMass;
                }

                foreach (Modification mod in relevantDatabaseMod.Value)
                {
                    double modifiedMass = basePrecursorMass + mod.MonoisotopicMass.Value;
                    if (modifiedMass <= MaxFragmentSize) //if the precursor is larger than the index allows, then don't add it
                    {
                        int precursorBin = (int)Math.Round(modifiedMass * FragmentBinsPerDalton);
                        emissions.Add((precursorBin, peptideId));
                    }
                }
            }
        }
    }
}