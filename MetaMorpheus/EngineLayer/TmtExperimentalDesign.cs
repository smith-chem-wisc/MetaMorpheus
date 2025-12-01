using FlashLFQ;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    public static class TmtExperimentalDesign
    {
        public const string Header = "File\tPlex\tSample Name\tTMT Channel\tCondition\tBiological Replicate\tFraction\tTechnical Replicate";

        // New validation helper: every (SampleName, BioRep, Fraction, TechRep) must be unique
        private static string ValidateUniqueSampleBioFracTech(IEnumerable<(string Sample, int Bio, int Fraction, int Tech)> tuples)
        {
            var duplicates = tuples
                .GroupBy(t => t)
                .Where(g => g.Count() > 1)
                .Select(g => g.Key)
                .ToList();

            if (duplicates.Count == 0)
                return null;

            var msgs = duplicates.Select(d =>
                $"Duplicate combination detected: Sample \"{d.Sample}\" Biorep {d.Bio} Fraction {d.Fraction} Techrep {d.Tech}");
            return string.Join(Environment.NewLine, msgs);
        }

        // RETURN: files where each file carries its plex annotations
        public static List<TmtFileInfo> Read(string tmtDesignPath, List<string> fullFilePathsWithExtension, out List<string> errors)
        {
            errors = new List<string>();
            var files = new List<TmtFileInfo>();

            // Collect per-plex annotations (unique by tag) and per-file (plex, frac, tech)
            var plexToAnnotations = new Dictionary<string, Dictionary<string, TmtPlexAnnotation>>(StringComparer.OrdinalIgnoreCase);
            var fileState = new Dictionary<string, (string Plex, int Fraction, int TechRep)>(StringComparer.OrdinalIgnoreCase);
            var fileStateConflicts = new Dictionary<string, HashSet<string>>(StringComparer.OrdinalIgnoreCase);

            // Collect tuples for uniqueness validation
            var uniqueTuples = new List<(string Sample, int Bio, int Fraction, int Tech)>();

            if (!File.Exists(tmtDesignPath))
            {
                errors.Add("TMT design file not found!");
                return files;
            }

            string[] lines;
            try
            {
                lines = File.ReadAllLines(tmtDesignPath);
            }
            catch (Exception ex)
            {
                errors.Add("Could not read TMT design file: " + ex.Message);
                return files;
            }

            if (lines.Length == 0 || !IsHeaderValid(lines[0]))
            {
                errors.Add("TMT design file header is invalid.");
                return files;
            }

            var headers = lines[0].Split('\t');
            int idxFile   = Array.FindIndex(headers, h => string.Equals(h.Trim(), "File", StringComparison.OrdinalIgnoreCase));
            int idxPlex   = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Plex", StringComparison.OrdinalIgnoreCase));
            int idxSample = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Sample Name", StringComparison.OrdinalIgnoreCase));
            int idxChannel= Array.FindIndex(headers, h => string.Equals(h.Trim(), "TMT Channel", StringComparison.OrdinalIgnoreCase));
            int idxCond   = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Condition", StringComparison.OrdinalIgnoreCase));
            int idxBio    = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Biological Replicate", StringComparison.OrdinalIgnoreCase));
            int idxFrac   = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Fraction", StringComparison.OrdinalIgnoreCase));
            int idxTech   = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Technical Replicate", StringComparison.OrdinalIgnoreCase));

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

                // Only consider lines for files actually in this run (mirrors ExperimentalDesign behavior)
                if (!fullFilePathsWithExtension.Any() ||
                    fullFilePathsWithExtension.Contains(full, StringComparer.OrdinalIgnoreCase))
                {
                    var plex    = cols[idxPlex].Trim();
                    var sample  = cols[idxSample].Trim();
                    var tag     = cols[idxChannel].Trim();
                    var cond    = cols[idxCond].Trim();

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

                    // Record tuple for uniqueness validation (ignore completely blank sample name)
                    if (!string.IsNullOrWhiteSpace(sample))
                        uniqueTuples.Add((sample, bio, frac, tech));

                    // Per-file consistency
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

                    // Plex annotations
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

            // Consistency errors
            foreach (var kvp in fileStateConflicts)
                errors.Add($"File '{kvp.Key}' has inconsistent Plex/Fraction/TechRep assignments: {string.Join(", ", kvp.Value)}");

            // Uniqueness validation of (Sample,Bio,Fraction,Tech)
            var uniquenessError = ValidateUniqueSampleBioFracTech(uniqueTuples);
            if (uniquenessError != null)
                errors.Add(uniquenessError);

            // Build file list with annotations embedded
            foreach (var kv in fileState)
            {
                var plex = kv.Value.Plex ?? string.Empty;
                var anns = plexToAnnotations.TryGetValue(plex, out var tags)
                    ? tags.Values.OrderBy(a => a.Tag, StringComparer.OrdinalIgnoreCase).ToList()
                    : new List<TmtPlexAnnotation>();

                files.Add(new TmtFileInfo(kv.Key, plex, kv.Value.Fraction, kv.Value.TechRep, anns));
            }

            // Ensure all provided files are defined in the TMT design (mirrors ExperimentalDesign behavior)
            if (fullFilePathsWithExtension != null && fullFilePathsWithExtension.Count > 0)
            {
                // normalize and compare case-insensitively using full paths
                var provided = new HashSet<string>(
                    fullFilePathsWithExtension.Select(p => { try { return Path.GetFullPath(p); } catch { return p; } }),
                    StringComparer.OrdinalIgnoreCase);

                var defined = new HashSet<string>(
                    fileState.Keys.Select(p => { try { return Path.GetFullPath(p); } catch { return p; } }),
                    StringComparer.OrdinalIgnoreCase);

                var notDefined = provided.Where(p => !defined.Contains(p)).ToList();
                if (notDefined.Any())
                {
                    errors.Add("Error: The TMT design did not contain the file(s): " + string.Join(", ", notDefined));
                }
            }

            return files;
        }

        public static string Write(List<TmtFileInfo> files)
        {
            if (files == null || files.Count == 0)
                throw new InvalidOperationException("No TMT files to write.");

            var dir = Directory.GetParent(files.First().FullFilePathWithExtension)!.FullName;
            var path = Path.Combine(dir, GlobalVariables.TmtExperimentalDesignFileName);

            using var sw = new StreamWriter(path);
            sw.WriteLine(Header);

            foreach (var file in files.OrderBy(f => f.Fraction).ThenBy(f => f.TechnicalReplicate))
            {
                var plex = file.Plex ?? "";
                var anns = file.Annotations ?? Array.Empty<TmtPlexAnnotation>();

                if (anns.Count > 0)
                {
                    foreach (var a in anns)
                    {
                        sw.WriteLine($"{file.FullFilePathWithExtension}\t{plex}\t{a.SampleName}\t{a.Tag}\t{a.Condition}\t{a.BiologicalReplicate}\t{file.Fraction}\t{file.TechnicalReplicate}");
                    }
                }
                else
                {
                    sw.WriteLine($"{file.FullFilePathWithExtension}\t{plex}\t\t\t\t\t{file.Fraction}\t{file.TechnicalReplicate}");
                }
            }

            return path;
        }

        private static bool IsHeaderValid(string headerLine)
        {
            if (string.Equals(headerLine, Header, StringComparison.Ordinal))
                return true;

            var parts = headerLine.Split('\t').Select(s => s.Trim().ToLowerInvariant()).ToHashSet();
            var required = Header.Split('\t').Select(s => s.Trim().ToLowerInvariant());
            return required.All(parts.Contains);
        }
    }

    public sealed class TmtFileInfo
    {
        public TmtFileInfo(string fullFilePathWithExtension, string plex, int fraction, int technicalReplicate, IReadOnlyList<TmtPlexAnnotation> annotations)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            Plex = plex ?? string.Empty;
            Fraction = fraction;
            TechnicalReplicate = technicalReplicate;
            Annotations = annotations ?? Array.Empty<TmtPlexAnnotation>();
        }

        public string FullFilePathWithExtension { get; }
        public string Plex { get; }
        public int Fraction { get; }             // 1-based
        public int TechnicalReplicate { get; }   // 1-based
        public IReadOnlyList<TmtPlexAnnotation> Annotations { get; } // All tags for this file's plex
        public override bool Equals(object obj)
        {
            if (base.Equals(obj))
            {
                return ((TmtFileInfo)obj).FullFilePathWithExtension.Equals(FullFilePathWithExtension);
            }

            return false;
        }

        public override int GetHashCode()
        {
            return FullFilePathWithExtension.GetHashCode();
        }

        public override string ToString()
        {
            return Path.GetFileName(FullFilePathWithExtension);
        }
    }

    public sealed class TmtPlexAnnotation
    {
        public string Tag { get; set; } = "";
        public string SampleName { get; set; } = "";
        public string Condition { get; set; } = "";
        public int BiologicalReplicate { get; set; }
    }
}