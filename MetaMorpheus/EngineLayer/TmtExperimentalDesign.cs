using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    // Mirrors ExperimentalDesign.cs style for TMT designs.
    // Reads a fixed "TmtDesign.txt" file and returns per-file info plus per-plex channel annotations.
    public static class TmtExperimentalDesign
    {
        public const string Header = "File\tPlex\tSample Name\tTMT Channel\tCondition\tBiological Replicate\tFraction\tTechnical Replicate";

        public static (List<TmtFileInfo> Files, Dictionary<string, List<TmtPlexAnnotation>> PlexAnnotations)
            Read(string tmtDesignPath, List<string> fullFilePathsWithExtension, out List<string> errors)
        {
            errors = new List<string>();
            var files = new List<TmtFileInfo>();
            var plexToAnnotations = new Dictionary<string, Dictionary<string, TmtPlexAnnotation>>(StringComparer.OrdinalIgnoreCase);

            if (!File.Exists(tmtDesignPath))
            {
                errors.Add("TMT design file not found!");
                return (files, ConvertAnnotations(plexToAnnotations));
            }

            string[] lines;
            try
            {
                lines = File.ReadAllLines(tmtDesignPath);
            }
            catch (Exception ex)
            {
                errors.Add("Could not read TMT design file: " + ex.Message);
                return (files, ConvertAnnotations(plexToAnnotations));
            }

            if (lines.Length == 0 || !IsHeaderValid(lines[0]))
            {
                errors.Add("TMT design file header is invalid.");
                return (files, ConvertAnnotations(plexToAnnotations));
            }

            // Parse header indices
            var headers = lines[0].Split('\t');
            int idxFile = Array.FindIndex(headers, h => string.Equals(h.Trim(), "File", StringComparison.OrdinalIgnoreCase));
            int idxPlex = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Plex", StringComparison.OrdinalIgnoreCase));
            int idxSample = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Sample Name", StringComparison.OrdinalIgnoreCase));
            int idxChannel = Array.FindIndex(headers, h => string.Equals(h.Trim(), "TMT Channel", StringComparison.OrdinalIgnoreCase));
            int idxCond = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Condition", StringComparison.OrdinalIgnoreCase));
            int idxBio = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Biological Replicate", StringComparison.OrdinalIgnoreCase));
            int idxFrac = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Fraction", StringComparison.OrdinalIgnoreCase));
            int idxTech = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Technical Replicate", StringComparison.OrdinalIgnoreCase));

            var fileState = new Dictionary<string, (string Plex, int Fraction, int TechRep)>(StringComparer.OrdinalIgnoreCase);
            var fileStateConflicts = new Dictionary<string, HashSet<string>>(StringComparer.OrdinalIgnoreCase);

            for (int i = 1; i < lines.Length; i++)
            {
                var line = lines[i];
                if (string.IsNullOrWhiteSpace(line)) continue;

                var cols = line.Split('\t');
                if (cols.Length < headers.Length)
                {
                    errors.Add($"Line {i + 1} has fewer columns than expected.");
                    continue;
                }

                var file = cols[idxFile].Trim();
                if (string.IsNullOrEmpty(file)) continue;

                string full;
                try { full = Path.GetFullPath(file); }
                catch { full = file; }

                // Only seed for files we currently care about (like ExperimentalDesign)
                if (!fullFilePathsWithExtension.Any() || fullFilePathsWithExtension.Contains(full, StringComparer.OrdinalIgnoreCase))
                {
                    var plex = cols[idxPlex].Trim();
                    var sample = cols[idxSample].Trim();
                    var tag = cols[idxChannel].Trim();
                    var cond = cols[idxCond].Trim();

                    if (!int.TryParse(cols[idxBio].Trim(), out var bio))
                    {
                        errors.Add($"Line {i + 1}: invalid Biological Replicate.");
                        continue;
                    }
                    if (!int.TryParse(cols[idxFrac].Trim(), out var frac) || frac < 1)
                    {
                        errors.Add($"Line {i + 1}: Fraction must be >= 1.");
                        continue;
                    }
                    if (!int.TryParse(cols[idxTech].Trim(), out var tech) || tech < 1)
                    {
                        errors.Add($"Line {i + 1}: Technical Replicate must be >= 1.");
                        continue;
                    }

                    // Per-file: ensure consistent (Plex, Fraction, TechRep) for all rows referencing this file
                    if (fileState.TryGetValue(full, out var state))
                    {
                        if (!string.Equals(state.Plex, plex, StringComparison.OrdinalIgnoreCase) ||
                            state.Fraction != frac ||
                            state.TechRep != tech)
                        {
                            if (!fileStateConflicts.TryGetValue(full, out var set))
                            {
                                set = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
                                fileStateConflicts[full] = set;
                            }
                            set.Add($"{plex}|{frac}|{tech}");
                        }
                    }
                    else
                    {
                        fileState[full] = (plex, frac, tech);
                    }

                    // Per-plex annotations, dedupe by tag; first occurrence wins
                    if (!plexToAnnotations.TryGetValue(plex, out var byTag))
                    {
                        byTag = new Dictionary<string, TmtPlexAnnotation>(StringComparer.OrdinalIgnoreCase);
                        plexToAnnotations[plex] = byTag;
                    }
                    if (!byTag.ContainsKey(tag))
                    {
                        byTag[tag] = new TmtPlexAnnotation
                        {
                            Tag = tag,
                            SampleName = sample,
                            Condition = cond,
                            BiologicalReplicate = bio
                        };
                    }
                }
            }

            // Report conflicts (if any)
            foreach (var kvp in fileStateConflicts)
            {
                errors.Add($"File '{kvp.Key}' has inconsistent Plex/Fraction/TechRep assignments in TmtDesign.txt: {string.Join(", ", kvp.Value)}");
            }

            // Build return list for files
            foreach (var kv in fileState)
            {
                files.Add(new TmtFileInfo(kv.Key, kv.Value.Plex, kv.Value.Fraction, kv.Value.TechRep));
            }

            return (files, ConvertAnnotations(plexToAnnotations));
        }

        private static Dictionary<string, List<TmtPlexAnnotation>> ConvertAnnotations(Dictionary<string, Dictionary<string, TmtPlexAnnotation>> src)
        {
            var result = new Dictionary<string, List<TmtPlexAnnotation>>(StringComparer.OrdinalIgnoreCase);
            foreach (var (plex, tags) in src)
            {
                // Keep insertion order by Tag alphabetically to be stable (could also leave as-is)
                var list = tags.Values.OrderBy(a => a.Tag, StringComparer.OrdinalIgnoreCase).ToList();
                result[plex] = list;
            }
            return result;
        }

        private static bool IsHeaderValid(string headerLine)
        {
            // exact match or case-insensitive contains of all required headers
            if (string.Equals(headerLine, Header, StringComparison.Ordinal))
                return true;

            var parts = headerLine.Split('\t').Select(s => s.Trim().ToLowerInvariant()).ToHashSet();
            var required = Header.Split('\t').Select(s => s.Trim().ToLowerInvariant());
            return required.All(parts.Contains);
        }
    }

    public sealed class TmtFileInfo
    {
        public TmtFileInfo(string fullFilePathWithExtension, string plex, int fraction, int technicalReplicate)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            Plex = plex ?? string.Empty;
            Fraction = fraction;
            TechnicalReplicate = technicalReplicate;
        }

        public string FullFilePathWithExtension { get; }
        public string Plex { get; }
        public int Fraction { get; }             // 1-based
        public int TechnicalReplicate { get; }   // 1-based
    }

    public sealed class TmtPlexAnnotation
    {
        public string Tag { get; set; } = "";
        public string SampleName { get; set; } = "";
        public string Condition { get; set; } = "";
        public int BiologicalReplicate { get; set; }
    }
}