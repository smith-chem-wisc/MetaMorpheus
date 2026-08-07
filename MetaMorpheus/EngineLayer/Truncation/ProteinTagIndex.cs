using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Proteomics;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// A k-mer index over a protein database for sequence-tag candidate filtering (doc §11.2.4): maps each
    /// I/L-normalized k-mer (k = tag length) to the proteins whose sequence contains it, so a scan's de-novo
    /// tags (<see cref="SequenceTagExtractor"/>) select a small candidate-protein set without materializing
    /// or scoring the whole database. Built in parallel over proteins (per the perf-parallelism convention),
    /// merging thread-local maps so each k-mer's protein-id list stays ascending.
    /// </summary>
    public class ProteinTagIndex
    {
        private readonly Dictionary<string, List<int>> _kmerToProteinIds;

        /// <summary>The proteins this index was built over; candidate ids returned by <see cref="GetCandidateProteinIds"/> index into this list.</summary>
        public IReadOnlyList<Protein> Proteins { get; }

        public int TagLength { get; }

        public ProteinTagIndex(IReadOnlyList<Protein> proteins, int tagLength, int maxThreads)
        {
            if (tagLength < 1) throw new ArgumentOutOfRangeException(nameof(tagLength));
            Proteins = proteins;
            TagLength = tagLength;
            _kmerToProteinIds = Build(proteins, tagLength, Math.Max(1, maxThreads));
        }

        private static Dictionary<string, List<int>> Build(IReadOnlyList<Protein> proteins, int k, int threads)
        {
            int chunks = Math.Max(1, Math.Min(threads, proteins.Count));
            var partials = new Dictionary<string, List<int>>[chunks];

            Parallel.For(0, chunks, chunk =>
            {
                int lo = (int)((long)proteins.Count * chunk / chunks);
                int hi = (int)((long)proteins.Count * (chunk + 1) / chunks);
                var local = new Dictionary<string, List<int>>();
                var seenInProtein = new HashSet<string>();
                var buffer = new char[k];

                for (int id = lo; id < hi; id++)
                {
                    string seq = proteins[id].BaseSequence;
                    if (seq == null || seq.Length < k) continue;
                    seenInProtein.Clear();

                    for (int start = 0; start + k <= seq.Length; start++)
                    {
                        for (int t = 0; t < k; t++) buffer[t] = SequenceTagExtractor.Normalize(seq[start + t]);
                        string kmer = new string(buffer);
                        if (!seenInProtein.Add(kmer)) continue; // a protein contributes each k-mer at most once
                        if (!local.TryGetValue(kmer, out List<int> ids)) local[kmer] = ids = new List<int>();
                        ids.Add(id);
                    }
                }

                partials[chunk] = local;
            });

            // Merge in chunk (ascending id-range) order so each posting list stays ascending.
            var global = new Dictionary<string, List<int>>();
            foreach (var part in partials)
            {
                foreach (KeyValuePair<string, List<int>> kv in part)
                {
                    if (!global.TryGetValue(kv.Key, out List<int> ids)) global[kv.Key] = kv.Value;
                    else ids.AddRange(kv.Value);
                }
            }
            return global;
        }

        /// <summary>Number of distinct k-mers indexed (diagnostic).</summary>
        public int DistinctKmerCount => _kmerToProteinIds.Count;

        /// <summary>
        /// Returns the ids of proteins supported by at least <paramref name="minTagHits"/> distinct query
        /// tags (tags whose length differs from <see cref="TagLength"/> are ignored). Higher
        /// <paramref name="minTagHits"/> trades recall for a smaller, more specific candidate set.
        /// </summary>
        public List<int> GetCandidateProteinIds(IEnumerable<string> tags, int minTagHits = 1)
        {
            var hitsPerProtein = new Dictionary<int, int>();
            // Count each DISTINCT query tag at most once; the contract is "distinct query tags",
            // so a tag repeated in the input must not inflate a protein toward minTagHits.
            foreach (string tag in new HashSet<string>(tags))
            {
                if (tag.Length != TagLength) continue;
                if (!_kmerToProteinIds.TryGetValue(tag, out List<int> ids)) continue;
                foreach (int id in ids)
                {
                    hitsPerProtein[id] = (hitsPerProtein.TryGetValue(id, out int v) ? v : 0) + 1;
                }
            }

            var candidates = new List<int>();
            foreach (KeyValuePair<int, int> kv in hitsPerProtein)
            {
                if (kv.Value >= minTagHits) candidates.Add(kv.Key);
            }
            return candidates;
        }
    }
}
