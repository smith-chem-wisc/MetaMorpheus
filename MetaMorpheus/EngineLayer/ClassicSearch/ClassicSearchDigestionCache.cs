#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading;
using Omics;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;

namespace EngineLayer.ClassicSearch
{
    /// <summary>
    /// The digestion products of a protein list, kept from one Classic search to the next within a task, so the list is
    /// digested once rather than once per spectra file (and, in calibration, once per round). Each protein's products are
    /// kept in their own slot, in the order digestion yields them, and a search reads them back in that order - exactly
    /// the sequence a fresh digestion would yield, since digestion depends only on the protein and these inputs.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Memory: the products are kept only while the estimated total stays within <see cref="BudgetBytes"/>, a quarter of
    /// the memory free when the cache is made (the machine's memory, a container limit or DOTNET_GCHeapHardLimit, whichever
    /// applies, less what is in use). A protein whose products would exceed it is digested as before, every time, and so
    /// is every protein after it; the products already kept are still used. On a 32 GB workstation a reviewed mouse or a
    /// human+yeast+E. coli database (9-13 million products, about 420 bytes each with mzLib's lazy full sequences) fits.
    /// </para>
    /// <para>
    /// Keeping the products costs garbage collection: the first search's products survive, so the collector copies them
    /// into the older generations instead of discarding them - about 7.5 s of collector pauses for the 13.5 million
    /// products of a human+yeast+E. coli database, against about 1.5 s saved on each later search (measured on Astral
    /// files: a 3-file search lost 4.5 s, a 6-file search broke even). So a task makes a cache only if it will run at
    /// least <see cref="MinimumSearches"/> searches with the same inputs.
    /// </para>
    /// <para>
    /// Results cannot change. The products of protein i are the same objects every search reads, in the same order the
    /// first search read them; which proteins are kept may vary with timing when the budget runs out, but a kept protein
    /// and a re-digested one yield equal products in the same order. Products are shared between files: digestion sets
    /// nothing on a product after making it except lazily computed values (mass, full sequence), and the one thing a
    /// search sets on one - PairedTargetDecoySequence, for a spectral-library search - turns the cache off. SILAC and
    /// turnover labels make digestion create new labelled proteins each time, so they turn it off as well.
    /// </para>
    /// <para>
    /// The cache serves only a search whose inputs are its own: the same protein list (same object, same count), equal
    /// digestion parameters, and the same fixed and variable modification lists. File-specific parameters that change
    /// digestion therefore get a fresh digestion, as before. The task drops the cache once its searches are done, before
    /// post-search analysis and PEP, so it does not outlive the searches it serves.
    /// </para>
    /// </remarks>
    public sealed class ClassicSearchDigestionCache
    {
        /// <summary>The share of the memory free at creation the products may take.</summary>
        public const double ShareOfFreeMemory = 0.25;

        /// <summary>
        /// The fewest Classic searches with the same inputs that make a cache worth its cost: it broke even at six (one
        /// filling it, five reading it), so eight leaves a margin.
        /// </summary>
        public const int MinimumSearches = 8;

        /// <summary>
        /// Estimated bytes per product besides its base sequence: the object, its modification dictionary and the slot
        /// holding it. Kept products measure 416-427 bytes each in all, base sequence included, with the mzLib that shares
        /// a product's base sequence between its forms and builds a digested form's full sequence only when it is read
        /// (bottom-up mouse and human+yeast+E. coli databases, top-down human). The estimate assumes that mzLib: a full
        /// sequence is not counted, since the cache cannot ask whether one was built without building it.
        /// </summary>
        private const long BytesPerProduct = 450;

        /// <summary>For tests: when set, the budget in bytes, instead of the share of free memory.</summary>
        internal static long? BudgetOverrideBytes { get; set; }

        /// <summary>For tests: when set, used instead of <see cref="MinimumSearches"/>.</summary>
        internal static int? MinimumSearchesOverride { get; set; }

        /// <summary>For tests: how many times, in this process, a search read a protein's kept products.</summary>
        internal static long ProteinsServed => Interlocked.Read(ref _proteinsServed);
        private static long _proteinsServed;

        // every cache made in this process, held weakly, so that tests can see which are still reachable
        private static readonly ConditionalWeakTable<ClassicSearchDigestionCache, object?> Instances = new();

        /// <summary>For tests: the caches made in this process that are still reachable.</summary>
        internal static List<ClassicSearchDigestionCache> ReachableInstances() => Instances.Select(kv => kv.Key).ToList();

        private readonly List<IBioPolymer> _proteins;
        private readonly int _proteinCount;
        private readonly IDigestionParams _digestionParams;
        private readonly List<Modification> _fixedModifications;
        private readonly List<Modification> _variableModifications;
        private readonly IBioPolymerWithSetMods[]?[] _products;
        private long _estimatedBytes;
        private int _full;

        private ClassicSearchDigestionCache(List<IBioPolymer> proteins, IDigestionParams digestionParams,
            List<Modification> fixedModifications, List<Modification> variableModifications, long budgetBytes)
        {
            _proteins = proteins;
            _proteinCount = proteins.Count;
            _digestionParams = digestionParams;
            _fixedModifications = fixedModifications;
            _variableModifications = variableModifications;
            _products = new IBioPolymerWithSetMods[]?[proteins.Count];
            BudgetBytes = budgetBytes;
            Instances.Add(this, null);
        }

        /// <summary>The estimated bytes the kept products may take.</summary>
        public long BudgetBytes { get; }

        /// <summary>The estimated bytes the kept products take.</summary>
        public long EstimatedBytes => Interlocked.Read(ref _estimatedBytes);

        /// <summary>How many proteins' products are kept.</summary>
        public int ProteinsKept => _products.Count(p => p != null);

        /// <summary>
        /// Lets go of every kept product; from then on the cache keeps nothing and every search digests afresh. A task calls
        /// this when its searches are done. Dropping its own reference is not enough: the search task runs post-search
        /// analysis (FDR, PEP, quantification) inside the same long method, which runs as tier-0 code, and tier-0 code
        /// reports every local and every spilled temporary of a method as live until the method returns - so the cache
        /// object itself can stay reachable from that frame after the variable holding it is set to null. Its products must
        /// not: they are most of the live heap, and PEP training collects that heap again and again.
        /// </summary>
        public void Release()
        {
            Volatile.Write(ref _full, 1);
            Array.Clear(_products);
            Interlocked.Exchange(ref _estimatedBytes, 0);
        }

        /// <summary>
        /// A cache for the Classic searches of a task that will run <paramref name="searches"/> of them with these inputs,
        /// or null when that is fewer than <see cref="MinimumSearches"/>, the inputs cannot be cached (SILAC or turnover
        /// labels), or no memory is free for it.
        /// </summary>
        public static ClassicSearchDigestionCache? Create(List<IBioPolymer> proteins, IDigestionParams digestionParams,
            List<Modification> fixedModifications, List<Modification> variableModifications, int searches,
            List<SilacLabel>? silacLabels = null, SilacLabel? startTurnoverLabel = null, SilacLabel? endTurnoverLabel = null)
        {
            if (proteins == null || digestionParams == null || fixedModifications == null || variableModifications == null
                || silacLabels != null || startTurnoverLabel != null || endTurnoverLabel != null
                || searches < (MinimumSearchesOverride ?? MinimumSearches))
            {
                return null;
            }

            long budget = BudgetOverrideBytes ?? FreeMemoryBudget();
            return budget > 0 ? new ClassicSearchDigestionCache(proteins, digestionParams, fixedModifications, variableModifications, budget) : null;
        }

        /// <summary>
        /// A quarter of the memory free now: the memory this process may use (physical memory, a container limit or the
        /// GC heap hard limit) less the larger of the machine's memory load and this process's managed heap.
        /// </summary>
        private static long FreeMemoryBudget()
        {
            GCMemoryInfo info = GC.GetGCMemoryInfo();
            long inUse = Math.Max(info.MemoryLoadBytes, GC.GetTotalMemory(false));
            long free = info.TotalAvailableMemoryBytes - inUse;
            return free > 0 ? (long)(free * ShareOfFreeMemory) : 0;
        }

        /// <summary>
        /// Whether a Classic search with these inputs may read its products from this cache.
        /// </summary>
        public bool Serves(List<IBioPolymer> proteins, IDigestionParams digestionParams, List<Modification> fixedModifications,
            List<Modification> variableModifications, List<SilacLabel>? silacLabels, (SilacLabel StartLabel, SilacLabel EndLabel)? turnoverLabels)
        {
            return ReferenceEquals(proteins, _proteins) && proteins.Count == _proteinCount
                && ReferenceEquals(fixedModifications, _fixedModifications) && ReferenceEquals(variableModifications, _variableModifications)
                && silacLabels == null && turnoverLabels == null
                && digestionParams != null
                && (ReferenceEquals(digestionParams, _digestionParams)
                    || (digestionParams.GetType() == _digestionParams.GetType() && digestionParams.Equals(_digestionParams)));
        }

        /// <summary>
        /// The products of the protein at <paramref name="proteinIndex"/>: those kept from an earlier search, or a fresh
        /// digestion with the searching engine's <paramref name="digestionParams"/> (equal to this cache's, see
        /// <see cref="Serves"/>), kept if the budget allows. Only one thread at a time may ask for a given protein.
        /// </summary>
        public IEnumerable<IBioPolymerWithSetMods> Digest(int proteinIndex, IDigestionParams digestionParams)
        {
            IBioPolymerWithSetMods[]? kept = Volatile.Read(ref _products[proteinIndex]);
            if (kept != null)
            {
                Interlocked.Increment(ref _proteinsServed);
                return kept;
            }

            IEnumerable<IBioPolymerWithSetMods> digestion = _proteins[proteinIndex].Digest(digestionParams, _fixedModifications, _variableModifications);
            if (Volatile.Read(ref _full) != 0)
            {
                return digestion;
            }

            IBioPolymerWithSetMods[] products = digestion.ToArray();
            if (Interlocked.Add(ref _estimatedBytes, EstimateBytes(products)) > BudgetBytes)
            {
                Volatile.Write(ref _full, 1);
                return products;
            }

            Volatile.Write(ref _products[proteinIndex], products);
            return products;
        }

        private long EstimateBytes(IBioPolymerWithSetMods[] products)
        {
            long bytes = 0;
            string? previousBaseSequence = null;
            foreach (IBioPolymerWithSetMods product in products)
            {
                bytes += BytesPerProduct;

                // the forms of one digestion product follow one another, and may share its base sequence string
                string baseSequence = product.BaseSequence;
                if (!ReferenceEquals(baseSequence, previousBaseSequence))
                {
                    bytes += 2L * baseSequence.Length + 24;
                    previousBaseSequence = baseSequence;
                }
            }
            return bytes;
        }
    }
}
