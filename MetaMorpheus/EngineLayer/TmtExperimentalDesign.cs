using FlashLFQ;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    public static class TmtExperimentalDesign
    {
        /// <summary>
        /// The header this class writes. "Sample Type" is the newest column and is OPTIONAL on read
        /// (see <see cref="RequiredColumns"/>), so a design file written before it existed still loads.
        /// </summary>
        public const string Header = "File\tPlex\tSample Name\tTMT Channel\tCondition\tBiological Replicate\tFraction\tTechnical Replicate\tSample Type";

        /// <summary>
        /// The columns a design file must carry to be readable. "Sample Type" is deliberately absent:
        /// requiring it would invalidate every TmtDesign.txt written before it was added, and a file
        /// without it is unambiguous - every channel is a study sample.
        /// </summary>
        public static readonly IReadOnlyList<string> RequiredColumns = new[]
        {
            "File", "Plex", "Sample Name", "TMT Channel", "Condition",
            "Biological Replicate", "Fraction", "Technical Replicate"
        };

        /// <summary>The column name carrying <see cref="TmtSampleType"/>.</summary>
        public const string SampleTypeColumn = "Sample Type";

        /// <summary>
        /// Every (Sample, BioRep, Fraction, TechRep) must be unique WITHIN A PLEX.
        /// </summary>
        /// <remarks>
        /// Per plex rather than per file, because a bridge channel is by definition the same material
        /// carried in more than one plex under one name -- which is what <see cref="TmtSampleType.Bridge"/>
        /// exists to describe. Collecting these across the whole file rejected exactly that design, and
        /// the only way through was to rename the bridges apart, which throws away the fact that they
        /// are the same material and so defeats the point of having them.
        ///
        /// Keyed per plex still catches the case worth catching: one plex describing the same sample
        /// twice.
        /// </remarks>
        private static string ValidateUniqueSampleBioFracTech(IEnumerable<(string Plex, string Sample, int Bio, int Fraction, int Tech)> tuples)
        {
            var duplicates = tuples
                .GroupBy(t => t)
                .Where(g => g.Count() > 1)
                .Select(g => g.Key)
                .ToList();

            if (duplicates.Count == 0)
                return null;

            var msgs = duplicates.Select(d =>
                $"Duplicate combination detected in plex \"{d.Plex}\": Sample \"{d.Sample}\" " +
                $"Biorep {d.Bio} Fraction {d.Fraction} Techrep {d.Tech}");
            return string.Join(Environment.NewLine, msgs);
        }

        // RETURN: files where each file carries its plex annotations
        /// <param name="fullFilePathsWithExtension">
        /// The files this run actually searched, used to ignore design rows naming anything else.
        /// Null or empty means "do not filter" -- the same meaning the guard in ToMzLibDesign gives it,
        /// and the reason this no longer throws on null.
        /// </param>
        public static List<TmtFileInfo> Read(string tmtDesignPath, List<string> fullFilePathsWithExtension, out List<string> errors)
        {
            errors = new List<string>();
            fullFilePathsWithExtension ??= new List<string>();
            var files = new List<TmtFileInfo>();

            // How many data rows the file had, and how many named a file in this run. A design whose
            // rows all fail to resolve is not an empty design -- it is a design pointed at the wrong
            // place -- and it used to return no errors at all.
            int dataRowsSeen = 0;
            int dataRowsMatched = 0;

            // Collect per-plex annotations (unique by tag) and per-file (plex, frac, tech)
            var plexToAnnotations = new Dictionary<string, Dictionary<string, TmtPlexAnnotation>>(StringComparer.OrdinalIgnoreCase);
            var fileState = new Dictionary<string, (string Plex, int Fraction, int TechRep)>(StringComparer.OrdinalIgnoreCase);
            var fileStateConflicts = new Dictionary<string, HashSet<string>>(StringComparer.OrdinalIgnoreCase);

            // Collect tuples for uniqueness validation
            var uniqueTuples = new List<(string Plex, string Sample, int Bio, int Fraction, int Tech)>();

            if (!File.Exists(tmtDesignPath))
            {
                errors.Add("TMT design file not found!");
                return files;
            }

            string designDirectory;
            try { designDirectory = Path.GetDirectoryName(Path.GetFullPath(tmtDesignPath)); }
            catch { designDirectory = null; }

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
            // -1 when the file predates the column; every channel then reads as a study sample.
            int idxType   = Array.FindIndex(headers, h => string.Equals(h.Trim(), SampleTypeColumn, StringComparison.OrdinalIgnoreCase));

            // How many cells a row needs to carry every required column.
            int minimumCells = new[] { idxFile, idxPlex, idxSample, idxChannel, idxCond, idxBio, idxFrac, idxTech }.Max() + 1;

            for (int i = 1; i < lines.Length; i++)
            {
                var line = lines[i];
                if (string.IsNullOrWhiteSpace(line)) continue;

                var cols = line.Split('\t');
                // A row must reach every REQUIRED column, not every column in the header. Comparing
                // against the full header width would reject a row whose only omission is the
                // optional trailing Sample Type cell — which is exactly what a design file written
                // before that column existed looks like.
                if (cols.Length < minimumCells)
                {
                    errors.Add($"Line {i + 1} has fewer columns than expected.");
                    continue;
                }

                var file = cols[idxFile].Trim();
                if (string.IsNullOrEmpty(file)) continue;

                // A row naming a file and no channel is the placeholder Write emits for a file with no
                // annotations yet, so that a part-authored design does not lose the file. Read used to
                // reject it -- the blank Biological Replicate failed int.TryParse, the row was skipped,
                // and the file was then reported as missing from the design. The placeholder therefore
                // lost exactly what it existed to preserve, and a part-authored design is the NORMAL
                // state while someone is filling one in through the GUI.
                bool channelless = string.IsNullOrWhiteSpace(cols[idxChannel]);

                // Resolve a relative name against the design file's own directory, not the process
                // working directory. A TmtDesign.txt written with bare file names sits beside the raw
                // files it names, and Path.GetFullPath alone would only match when the process happened
                // to be running in that folder -- so the same design file worked or silently matched
                // nothing depending on where MetaMorpheus was launched from.
                string full = ResolveAgainstDesignFile(file, designDirectory);

                // Only consider lines for files actually in this run (mirrors ExperimentalDesign behavior)
                bool inThisRun = !fullFilePathsWithExtension.Any() ||
                    fullFilePathsWithExtension.Contains(full, StringComparer.OrdinalIgnoreCase);

                dataRowsSeen++;
                if (inThisRun)
                {
                    dataRowsMatched++;

                    if (channelless)
                    {
                        // Register the file with no annotations and move on. Fraction and technical
                        // replicate still parse when present, so a placeholder carrying them keeps them.
                        int.TryParse(cols[idxFrac].Trim(), out int placeholderFrac);
                        int.TryParse(cols[idxTech].Trim(), out int placeholderTech);
                        if (!fileState.ContainsKey(full))
                        {
                            fileState[full] = (cols[idxPlex].Trim(),
                                placeholderFrac < 1 ? 1 : placeholderFrac,
                                placeholderTech < 1 ? 1 : placeholderTech);
                        }
                        continue;
                    }

                    var plex    = cols[idxPlex].Trim();
                    var sample  = cols[idxSample].Trim();
                    var tag     = cols[idxChannel].Trim();
                    var cond    = cols[idxCond].Trim();

                    var rawType = idxType >= 0 && idxType < cols.Length ? cols[idxType].Trim() : string.Empty;
                    if (!TryParseSampleType(rawType, out var sampleType))
                    {
                        errors.Add($"Line {i + 1}: '{rawType}' is not a recognised Sample Type. " +
                                   $"Use one of: {string.Join(", ", SampleTypeNames)}.");
                        continue;
                    }

                    // >= 1, matching Fraction and Technical Replicate. A bare TryParse accepted 0 and
                    // -3 silently, and ToMzLibDesign passes the value straight through to
                    // ISampleInfo.BiologicalReplicate. If these are ever meant to be 0-based, all three
                    // columns have to change together rather than one drifting from the other two.
                    if (!int.TryParse(cols[idxBio].Trim(), out var bio) || bio < 1)
                    {
                        errors.Add($"Line {i + 1}: Biological Replicate must be >= 1.");
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
                        uniqueTuples.Add((plex, sample, bio, frac, tech));

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
                    if (!byTag.TryGetValue(tag, out var existing))
                    {
                        byTag[tag] = new TmtPlexAnnotation
                        {
                            Tag = tag,
                            SampleName = sample,
                            Condition = cond,
                            BiologicalReplicate = bio,
                            SampleType = sampleType
                        };
                    }
                    // Fractions of one plex legitimately repeat the same channel-to-sample map, which is
                    // why this is keyed per plex and why a repeat is collapsed rather than refused. But
                    // collapsing rows that DISAGREE silently discards one of them: two rows naming
                    // channel 126 as different samples kept the first and said nothing. The uniqueness
                    // rule above cannot catch it either -- two different sample names are two different
                    // tuples, so it only fires when the duplicates agree, which is the harmless case.
                    else if (!string.Equals(existing.SampleName, sample, StringComparison.OrdinalIgnoreCase)
                             || !string.Equals(existing.Condition, cond, StringComparison.OrdinalIgnoreCase)
                             || existing.BiologicalReplicate != bio
                             || existing.SampleType != sampleType)
                    {
                        errors.Add(
                            $"Line {i + 1}: plex \"{plex}\" channel \"{tag}\" is described twice and the " +
                            $"descriptions disagree — \"{existing.SampleName}\"/\"{existing.Condition}\" " +
                            $"versus \"{sample}\"/\"{cond}\".");
                    }
                }
            }

            // A design that named files, none of which belong to this run, is a mistake rather than an
            // empty design -- most often a design written for different data, or one whose relative
            // paths do not resolve. Reporting nothing here let the command line print success over an
            // annotation set that was silently empty.
            if (dataRowsSeen > 0 && dataRowsMatched == 0 && fullFilePathsWithExtension.Any())
            {
                errors.Add($"None of the {dataRowsSeen} row(s) in the TMT design name a file in this run. " +
                           "Check that the File column matches the spectra files being searched.");
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

        /// <summary>
        /// Resolves a File cell to a full path. A rooted path is taken as written; a bare or relative
        /// name is resolved against the design file's own directory first, since that is where a
        /// TmtDesign.txt written alongside its raw files expects to be read from, and only then
        /// against the process working directory.
        /// </summary>
        private static string ResolveAgainstDesignFile(string file, string designDirectory)
        {
            try
            {
                if (Path.IsPathRooted(file))
                {
                    return Path.GetFullPath(file);
                }

                if (!string.IsNullOrEmpty(designDirectory))
                {
                    string besideDesign = Path.GetFullPath(Path.Combine(designDirectory, file));
                    if (File.Exists(besideDesign))
                    {
                        return besideDesign;
                    }
                }

                return Path.GetFullPath(file);
            }
            catch
            {
                return file;
            }
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
                        sw.WriteLine($"{file.FullFilePathWithExtension}\t{plex}\t{a.SampleName}\t{a.Tag}\t{a.Condition}\t{a.BiologicalReplicate}\t{file.Fraction}\t{file.TechnicalReplicate}\t{ToDesignFileValue(a.SampleType)}");
                    }
                }
                else
                {
                    sw.WriteLine($"{file.FullFilePathWithExtension}\t{plex}\t\t\t\t\t{file.Fraction}\t{file.TechnicalReplicate}\t");
                }
            }

            return path;
        }

        #region Projection to mzLib

        /// <summary>
        /// Projects a design onto mzLib's <see cref="IExperimentalDesign"/>, the contract
        /// <c>QuantificationEngine</c> consumes.
        ///
        /// The design file is only a carrier; this is the thing quantification actually runs on, so
        /// a second carrier (an SDRF-Proteomics file, say) means a second method like this one rather
        /// than a different engine.
        /// </summary>
        /// <param name="files">The parsed design, from <see cref="Read"/>.</param>
        /// <param name="tag">
        /// The isobaric tag the search used. It supplies the channel order and the reporter m/z
        /// values, neither of which the design file states.
        /// </param>
        /// <param name="errors">Everything that made the projection impossible or lossy.</param>
        /// <returns>
        /// A design keyed by file name WITH extension, matching what
        /// <c>QuantificationEngine.PivotByFile</c> looks up. Null when <paramref name="errors"/> is
        /// non-empty, so a caller cannot half-use a broken design.
        /// </returns>
        /// <remarks>
        /// Two contracts are load-bearing here, and getting either wrong mislabels every channel
        /// silently rather than failing:
        ///
        /// 1. ORDER. mzLib requires the ISampleInfo array to line up positionally with
        ///    <c>ISpectralMatch.Intensities</c>. MetaMorpheus fills that array in
        ///    <see cref="IsobaricMassTag.GetReporterIonIntensities"/>, which walks
        ///    <see cref="IsobaricMassTag.ReporterIonMzs"/> — sorted ascending. So the array emitted
        ///    here is in tag order, NOT the order rows happened to appear in the design file.
        ///
        /// 2. LENGTH. That same method returns one intensity per channel the tag defines, whether or
        ///    not the design annotated it. So every tag channel gets an entry; an unannotated one
        ///    becomes an <see cref="TmtSampleType.Empty"/> channel rather than being skipped, which
        ///    would shift every later channel by one.
        /// </remarks>
        public static IExperimentalDesign ToMzLibDesign(
            IEnumerable<TmtFileInfo> files,
            IsobaricMassTag tag,
            out List<string> errors)
        {
            errors = new List<string>();

            var fileList = files?.ToList();
            if (fileList == null || fileList.Count == 0)
            {
                errors.Add("No TMT design entries to convert.");
                return null;
            }
            if (tag == null)
            {
                errors.Add("No isobaric mass tag was supplied, so reporter ion m/z values are unknown. " +
                           "Check that the search's multiplex label is set.");
                return null;
            }

            var channelLabels = IsobaricMassTag.GetReporterIonLabels(tag.TagType);
            if (channelLabels == null || channelLabels.Count == 0)
            {
                errors.Add($"No channel labels are known for isobaric tag {tag.TagType}.");
                return null;
            }
            if (tag.ReporterIonMzs == null || tag.ReporterIonMzs.Length != channelLabels.Count)
            {
                errors.Add($"Isobaric tag {tag.TagType} defines {channelLabels.Count} channels but " +
                           $"{tag.ReporterIonMzs?.Length ?? 0} reporter ion m/z values. The modification " +
                           "definition and the channel list disagree.");
                return null;
            }

            // Plex is a free-text label in the design file but an int on IsobaricQuantSampleInfo.
            // Assign by sorted order so the same design always yields the same ids.
            var plexIds = fileList
                .Select(f => f.Plex ?? string.Empty)
                .Distinct(StringComparer.OrdinalIgnoreCase)
                .OrderBy(p => p, StringComparer.OrdinalIgnoreCase)
                .Select((plex, index) => (plex, id: index + 1))
                .ToDictionary(x => x.plex, x => x.id, StringComparer.OrdinalIgnoreCase);

            var design = new Dictionary<string, ISampleInfo[]>(StringComparer.OrdinalIgnoreCase);

            foreach (var file in fileList)
            {
                string fileName = Path.GetFileName(file.FullFilePathWithExtension);
                if (string.IsNullOrEmpty(fileName))
                {
                    errors.Add($"A design entry has no file name: '{file.FullFilePathWithExtension}'.");
                    continue;
                }
                if (design.ContainsKey(fileName))
                {
                    errors.Add($"Two design entries share the file name '{fileName}'. " +
                               "Quantification keys channels by file name, so the names must be distinct.");
                    continue;
                }

                var byTag = new Dictionary<string, TmtPlexAnnotation>(StringComparer.OrdinalIgnoreCase);
                foreach (var annotation in file.Annotations ?? Array.Empty<TmtPlexAnnotation>())
                {
                    if (!channelLabels.Contains(annotation.Tag, StringComparer.OrdinalIgnoreCase))
                    {
                        errors.Add($"File '{fileName}' annotates channel '{annotation.Tag}', which is not a " +
                                   $"channel of {tag.TagType}. Expected one of: {string.Join(", ", channelLabels)}.");
                        continue;
                    }
                    byTag[annotation.Tag] = annotation;
                }

                var samples = new ISampleInfo[channelLabels.Count];
                for (int i = 0; i < channelLabels.Count; i++)
                {
                    string label = channelLabels[i];
                    byTag.TryGetValue(label, out var annotation);

                    var sampleType = annotation?.SampleType ?? TmtSampleType.Empty;

                    samples[i] = new IsobaricQuantSampleInfo(
                        fullFilePathWithExtension: file.FullFilePathWithExtension,
                        condition: annotation?.Condition ?? string.Empty,
                        biologicalReplicate: annotation?.BiologicalReplicate ?? 0,
                        technicalReplicate: file.TechnicalReplicate,
                        fraction: file.Fraction,
                        plexId: plexIds[file.Plex ?? string.Empty],
                        channelLabel: label,
                        reporterIonMz: tag.ReporterIonMzs[i],
                        isReferenceChannel: IsReferenceChannel(sampleType));
                }

                design[fileName] = samples;
            }

            return errors.Count > 0 ? null : new TmtMzLibExperimentalDesign(design);
        }

        /// <summary>
        /// The <see cref="IExperimentalDesign"/> returned by <see cref="ToMzLibDesign"/>.
        /// </summary>
        private sealed class TmtMzLibExperimentalDesign : IExperimentalDesign
        {
            public TmtMzLibExperimentalDesign(Dictionary<string, ISampleInfo[]> design)
            {
                FileNameSampleInfoDictionary = design;
            }

            public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }
        }

        #endregion

        #region Sample Type

        /// <summary>
        /// The accepted spellings, for error messages and for the GUI's drop-down.
        /// </summary>
        public static readonly IReadOnlyList<string> SampleTypeNames = new[]
        {
            "study sample", "reference", "bridge", "carrier", "empty"
        };

        /// <summary>
        /// Reads a Sample Type cell. An empty cell is a study sample, which is what a design file
        /// written before the column existed means. Comparison is case-insensitive, and the
        /// SDRF-Proteomics spellings are accepted as-is so the two formats agree.
        /// </summary>
        public static bool TryParseSampleType(string value, out TmtSampleType sampleType)
        {
            sampleType = TmtSampleType.StudySample;

            if (string.IsNullOrWhiteSpace(value))
                return true;

            switch (value.Trim().ToLowerInvariant())
            {
                case "study sample":
                case "studysample":
                case "sample":
                    sampleType = TmtSampleType.StudySample; return true;
                case "reference":
                case "reference channel":
                case "pooled":
                    sampleType = TmtSampleType.Reference; return true;
                case "bridge":
                case "bridge channel":
                    sampleType = TmtSampleType.Bridge; return true;
                case "carrier":
                case "carrier channel":
                    sampleType = TmtSampleType.Carrier; return true;
                case "empty":
                case "unused":
                    sampleType = TmtSampleType.Empty; return true;
                default:
                    return false;
            }
        }

        /// <summary>The spelling written back out, so a read/write round trip is stable.</summary>
        public static string ToDesignFileValue(TmtSampleType sampleType) => sampleType switch
        {
            TmtSampleType.StudySample => "study sample",
            TmtSampleType.Reference => "reference",
            TmtSampleType.Bridge => "bridge",
            TmtSampleType.Carrier => "carrier",
            TmtSampleType.Empty => "empty",
            _ => "study sample"
        };

        /// <summary>
        /// True when this channel normalizes the others rather than being quantified against them.
        /// This is the single bit mzLib's <see cref="IsobaricQuantSampleInfo.IsReferenceChannel"/>
        /// carries, and the only thing <c>ReferenceChannelNormalization</c> asks about.
        /// </summary>
        public static bool IsReferenceChannel(TmtSampleType sampleType) =>
            sampleType is TmtSampleType.Reference or TmtSampleType.Bridge;

        #endregion

        private static bool IsHeaderValid(string headerLine)
        {
            if (string.Equals(headerLine, Header, StringComparison.Ordinal))
                return true;

            var parts = headerLine.Split('\t').Select(s => s.Trim().ToLowerInvariant()).ToHashSet();
            return RequiredColumns.Select(c => c.ToLowerInvariant()).All(parts.Contains);
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
        /// <summary>
        /// Identity is the file path, compared the way <see cref="TmtExperimentalDesign.Read"/>
        /// compares it -- case-insensitively, because that is how the design file keys files.
        /// </summary>
        public override bool Equals(object obj)
        {
            return obj is TmtFileInfo other
                && string.Equals(FullFilePathWithExtension, other.FullFilePathWithExtension,
                    StringComparison.OrdinalIgnoreCase);
        }

        public override int GetHashCode()
        {
            return StringComparer.OrdinalIgnoreCase.GetHashCode(FullFilePathWithExtension ?? string.Empty);
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

        /// <summary>
        /// What this channel is for. Defaults to <see cref="TmtSampleType.StudySample"/>, which is
        /// what a design file written before the column existed means.
        /// </summary>
        public TmtSampleType SampleType { get; set; } = TmtSampleType.StudySample;
    }

    /// <summary>
    /// What a channel is for within its plex.
    ///
    /// The names and spellings are taken from SDRF-Proteomics'
    /// <c>characteristics[sample type]</c> rather than invented here, so that a design read from an
    /// SDRF file later maps through the same switch as one read from a TmtDesign.txt. That is also
    /// why this is an enum rather than the single boolean mzLib consumes: mzLib's
    /// <c>IsobaricQuantSampleInfo.IsReferenceChannel</c> only needs to know whether a channel
    /// normalizes the others, but the design file is worth keeping more faithful than that.
    /// </summary>
    public enum TmtSampleType
    {
        /// <summary>An ordinary experimental channel. The default.</summary>
        StudySample,

        /// <summary>A pooled reference channel, used to normalize the other channels in its plex.</summary>
        Reference,

        /// <summary>A bridge channel, shared across plexes to make them comparable.</summary>
        Bridge,

        /// <summary>A carrier channel: high-load material added to boost signal, never quantified.</summary>
        Carrier,

        /// <summary>An empty channel, present in the plex but carrying no sample.</summary>
        Empty
    }
}