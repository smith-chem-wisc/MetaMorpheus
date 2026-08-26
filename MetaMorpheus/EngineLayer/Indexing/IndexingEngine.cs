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
            Parallel.ForEach(threads, (i) =>
            {
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

                lock (peptides)
                {
                    peptides.AddRange(localPeptides);
                }
            });

            // sort peptides by mass
            peptides.Sort((x, y) => x.MonoisotopicMass.CompareTo(y.MonoisotopicMass));

            int binCount = (int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1;

            //create precursor index (if specified)
            MassBinEmissionRun precursorBaseEmissions = null;
            if (GeneratePrecursorIndex)
            {
                precursorBaseEmissions = CreatePrecursorBaseEmissions(peptides);
            }
            bool addInteriorTerminalModsToPrecursorIndex = GeneratePrecursorIndex && CommonParameters.DigestionParams.DigestionAgent.Name.Contains("single");
            List<Modification> terminalModifications = addInteriorTerminalModsToPrecursorIndex ?
                NonSpecificEnzymeSearchEngine.GetVariableTerminalMods(CommonParameters.DigestionParams.FragmentationTerminus, VariableModifications) :
                null;

            // Fragment every peptide once, emitting into per-block runs. The runs are concatenated in
            // block order below, which reproduces the emission sequence a single sequential pass would
            // have produced -- and therefore the exact per-bin ordering the old List<int>[] build gave.
            // That ordering is load-bearing: BinarySearchBinForPrecursorIndex binary-searches a bin by
            // peptide mass, which only works while peptide ids ascend within the bin.
            // Blocks are independent -- each one fragments its own peptides into its own run -- so they
            // run in parallel. Digestion above this point was already parallel; fragmenting was not.
            List<PeptideBlock> blocks = PartitionPeptidesIntoBlocks(peptides.Count);
            var fragmentEmissions = new MassBinEmissionRun[blocks.Count];
            var interiorTerminalModEmissions = addInteriorTerminalModsToPrecursorIndex ? new MassBinEmissionRun[blocks.Count] : null;
            var fragmentationProgress = new FragmentationProgress();

            Parallel.For(0, blocks.Count, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, CommonParameters.MaxThreadsToUsePerFile) }, blockIndex =>
            {
                FragmentPeptideBlock(blocks[blockIndex], blockIndex, peptides, addInteriorTerminalModsToPrecursorIndex,
                    terminalModifications, fragmentEmissions, interiorTerminalModEmissions, fragmentationProgress);
            });

            MassBinIndex fragmentIndex;
            MassBinIndex precursorIndex = null;

            try
            {
                fragmentIndex = MassBinIndex.Build(binCount, fragmentEmissions);

                if (GeneratePrecursorIndex)
                {
                    // Base precursor masses first, then the interior terminal mods, matching the order
                    // in which the old two-pass build appended them.
                    var precursorEmissions = new List<MassBinEmissionRun> { precursorBaseEmissions };
                    if (addInteriorTerminalModsToPrecursorIndex)
                    {
                        precursorEmissions.AddRange(interiorTerminalModEmissions);
                    }
                    precursorIndex = MassBinIndex.Build(binCount, precursorEmissions);
                }
            }
            catch (OutOfMemoryException)
            {
                throw new MetaMorpheusException("Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }

            return new IndexingResults(peptides, fragmentIndex, precursorIndex, this);
        }

        /// <summary>
        /// Counts fragmented peptides across the parallel blocks. Each whole percent is announced by
        /// exactly one thread -- whichever wins the compare-exchange that advances the counter.
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
        /// Splits the peptide ids into contiguous blocks. Emissions are concatenated in block order, so
        /// the index that comes out does not depend on how many blocks there are.
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
        /// the same fragmentation so that peptides are not fragmented twice.
        /// </summary>
        private void FragmentPeptideBlock(PeptideBlock block, int blockIndex, List<PeptideWithSetModifications> peptides,
            bool addInteriorTerminalModsToPrecursorIndex, List<Modification> terminalModifications,
            MassBinEmissionRun[] fragmentEmissions, MassBinEmissionRun[] interiorTerminalModEmissions,
            FragmentationProgress fragmentationProgress)
        {
            int blockLength = block.End - block.Start;
            var fragmentRun = new MassBinEmissionRun(block.Start, blockLength);
            MassBinEmissionRun interiorRun = addInteriorTerminalModsToPrecursorIndex ? new MassBinEmissionRun(block.Start, blockLength) : null;

            List<Product> fragments = new List<Product>();

            for (int peptideId = block.Start; peptideId < block.End; peptideId++)
            {
                fragmentRun.BeginPeptide();
                interiorRun?.BeginPeptide();

                peptides[peptideId].Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, fragments, CommonParameters.FragmentationParameters);

                foreach (var theoreticalFragment in fragments)
                {
                    double theoreticalFragmentMass = theoreticalFragment.NeutralMass;

                    //if low res round
                    if (CommonParameters.DissociationType == MassSpectrometry.DissociationType.LowCID)
                    {
                        theoreticalFragmentMass = Math.Round(theoreticalFragmentMass / 1.0005079, 0) * 1.0005079;
                    }

                    if (theoreticalFragmentMass < MaxFragmentSize && theoreticalFragmentMass > 0)
                    {
                        fragmentRun.Add((int)Math.Round(theoreticalFragmentMass * FragmentBinsPerDalton));
                    }
                }

                //Add terminal mods if needed (do it here rather than earlier so that we don't have to fragment twice)
                if (addInteriorTerminalModsToPrecursorIndex)
                {
                    AddInteriorTerminalModsToPrecursorIndex(interiorRun, fragments, peptides[peptideId], terminalModifications);
                }

                if (fragmentationProgress.AdvanceOnePeptide(peptides.Count, out int percentProgress))
                {
                    ReportProgress(new ProgressEventArgs(percentProgress, "Fragmenting peptides...", NestedIds));
                }
            }

            fragmentRun.Complete();
            interiorRun?.Complete();

            fragmentEmissions[blockIndex] = fragmentRun;
            if (interiorTerminalModEmissions != null)
            {
                interiorTerminalModEmissions[blockIndex] = interiorRun;
            }
        }

        /// <summary>
        /// Emits the base precursor bin for every peptide, in peptide id order. These land in the
        /// precursor index ahead of any interior terminal mod entries.
        /// </summary>
        private MassBinEmissionRun CreatePrecursorBaseEmissions(List<PeptideWithSetModifications> peptidesSortedByMass)
        {
            var run = new MassBinEmissionRun(0, peptidesSortedByMass.Count);

            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(0, "Creating precursor index...", NestedIds));

            //Add all the precursors
            for (int i = 0; i < peptidesSortedByMass.Count; i++)
            {
                run.BeginPeptide();

                double mass = peptidesSortedByMass[i].MonoisotopicMass;
                if (!double.IsNaN(mass))
                {
                    if (mass > MaxFragmentSize) //if the precursor is larger than the index allows, then stop adding precursors
                    {
                        break;
                    }

                    run.Add((int)Math.Round(mass * FragmentBinsPerDalton));
                }
                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Creating precursor index...", NestedIds));
                }
            }

            run.Complete();
            return run;
        }

        //add possible protein/peptide terminal modifications that aren't on the terminal amino acids
        //The purpose is for terminal mods that are contained WITHIN the Single peptide
        private void AddInteriorTerminalModsToPrecursorIndex(MassBinEmissionRun precursorRun, List<Product> fragmentMasses, PeptideWithSetModifications peptide, List<Modification> variableModifications)
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
                        precursorRun.Add((int)Math.Round(modifiedMass * FragmentBinsPerDalton));
                    }
                }
            }
        }
    }
}