using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    /// <summary>What <see cref="ExperimentalDesignCheck.Check"/> decides about a run.</summary>
    public enum ExperimentalDesignVerdict
    {
        /// <summary>The run may start: a design was read, or none is needed.</summary>
        Proceed,

        /// <summary>The run must not start.</summary>
        Refuse,

        /// <summary>A design file is there but cannot be read, and nothing needs it: ask whether to run without it.</summary>
        AskToContinueWithoutDesign,
    }

    /// <summary>Why a run is refused or asked about.</summary>
    public enum ExperimentalDesignProblem
    {
        None,

        /// <summary>ExperimentalDesign.tsv and TmtDesign.txt are both in the folder.</summary>
        BothDesignFiles,

        /// <summary>No design file, and a search task normalizes.</summary>
        MissingForNormalization,

        /// <summary>The design file could not be read.</summary>
        Unreadable,
    }

    /// <summary>The outcome of <see cref="ExperimentalDesignCheck.Check"/>, for a front end to show.</summary>
    public sealed class ExperimentalDesignCheckResult
    {
        internal ExperimentalDesignCheckResult(ExperimentalDesignVerdict verdict, ExperimentalDesignProblem problem,
            string designFileName, List<string> messages, List<string> errors)
        {
            Verdict = verdict;
            Problem = problem;
            DesignFileName = designFileName;
            Messages = messages;
            Errors = errors;
        }

        public ExperimentalDesignVerdict Verdict { get; }

        public ExperimentalDesignProblem Problem { get; }

        /// <summary>The design file the folder holds (ExperimentalDesign.tsv or TmtDesign.txt); null when it holds neither.</summary>
        public string DesignFileName { get; }

        /// <summary>
        /// What to tell the user. Proceeding: "Read ... successfully" and any warnings about the design, or nothing when
        /// there is no design. Refused: the reason, or every error reading the design. Asking: the question, then the
        /// first error.
        /// </summary>
        public IReadOnlyList<string> Messages { get; }

        /// <summary>Every error reading the design file; empty unless <see cref="Problem"/> is <see cref="ExperimentalDesignProblem.Unreadable"/>.</summary>
        public IReadOnlyList<string> Errors { get; }

        /// <summary>The question to ask when <see cref="Verdict"/> is <see cref="ExperimentalDesignVerdict.AskToContinueWithoutDesign"/>.</summary>
        public string Question { get; internal init; }
    }

    /// <summary>
    /// Decides, before a run starts, what its design files allow: the command line (CMD Program.ResolveExperimentalDesign)
    /// and the GUI's Run button both call this, so a folder one refuses the other refuses too. Before, the GUI read only
    /// ExperimentalDesign.tsv: it never looked at TmtDesign.txt, did not refuse both files in one folder, and refused a
    /// normalizing TMT run whose only design was TmtDesign.txt. Each front end only shows the result (the command line
    /// prints and asks y/n; the GUI shows a dialog), and neither ever deletes a design file.
    /// </summary>
    public static class ExperimentalDesignCheck
    {
        /// <summary>
        /// What the design files in <paramref name="designDirectory"/> allow a run of <paramref name="spectraFilePaths"/>.
        /// </summary>
        /// <param name="designDirectory">The first spectra file's folder, where a run reads its design.</param>
        /// <param name="spectraFilePaths">Full paths of the spectra files the run will search.</param>
        /// <param name="normalizationRequested">True when any search task normalizes, which makes a design mandatory.</param>
        public static ExperimentalDesignCheckResult Check(string designDirectory, List<string> spectraFilePaths,
            bool normalizationRequested)
        {
            string pathToExperDesign = Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName);
            string pathToTmtDesign = Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName);

            bool hasClassicDesign = File.Exists(pathToExperDesign);
            bool hasTmtDesign = File.Exists(pathToTmtDesign);

            // Only one may be present. Two files carrying replicate structure can drift apart, and nothing here could say
            // which of them the user meant.
            if (hasClassicDesign && hasTmtDesign)
            {
                return new ExperimentalDesignCheckResult(ExperimentalDesignVerdict.Refuse, ExperimentalDesignProblem.BothDesignFiles,
                    null, new List<string> { "Both ExperimentalDesign and TmtDesign.txt were found. Only one design file type may be present." },
                    new List<string>());
            }

            if (!hasClassicDesign && !hasTmtDesign)
            {
                // A design is optional until something needs it. Normalization is that something.
                return normalizationRequested
                    ? new ExperimentalDesignCheckResult(ExperimentalDesignVerdict.Refuse, ExperimentalDesignProblem.MissingForNormalization,
                        null, new List<string> { "No experimental design file present. Normalization requires a design (ExperimentalDesign.tsv or TmtDesign.txt)." },
                        new List<string>())
                    : new ExperimentalDesignCheckResult(ExperimentalDesignVerdict.Proceed, ExperimentalDesignProblem.None,
                        null, new List<string>(), new List<string>());
            }

            string designPath = hasClassicDesign ? pathToExperDesign : pathToTmtDesign;
            string designFileName = Path.GetFileName(designPath);
            List<string> errors;
            var warnings = new List<string>();

            if (hasClassicDesign)
            {
                var design = ExperimentalDesign.ReadExperimentalDesign(designPath, spectraFilePaths, out errors);
                if (!errors.Any())
                {
                    warnings = ExperimentalDesign.GetWarningsInExperimentalDesign(design);
                }
            }
            else
            {
                var tmtFiles = TmtExperimentalDesign.Read(designPath, spectraFilePaths, out errors);
                if (!errors.Any())
                {
                    warnings = TmtExperimentalDesign.GetWarningsInDesign(tmtFiles);
                }
            }

            if (!errors.Any())
            {
                return new ExperimentalDesignCheckResult(ExperimentalDesignVerdict.Proceed, ExperimentalDesignProblem.None,
                    designFileName, new[] { "Read " + designFileName + " successfully" }.Concat(warnings).ToList(), new List<string>());
            }

            // Normalizing against a design that could not be read would silently produce different numbers, so this is the
            // one case that cannot be recovered from.
            if (normalizationRequested)
            {
                return new ExperimentalDesignCheckResult(ExperimentalDesignVerdict.Refuse, ExperimentalDesignProblem.Unreadable,
                    designFileName, new List<string>(errors), errors);
            }

            // Otherwise the run may continue without reading the design, and the file is left alone (#2256, #2780: a design
            // that fails to parse is the user's file, not ours to discard).
            string question = (hasClassicDesign ? "An experimental design file" : "A TMT design file")
                + " was found, but an error occurred reading it. Continue without an experimental design?";
            return new ExperimentalDesignCheckResult(ExperimentalDesignVerdict.AskToContinueWithoutDesign, ExperimentalDesignProblem.Unreadable,
                designFileName, new List<string> { question, "First error: " + errors.First() }, errors)
            {
                Question = question,
            };
        }
    }
}
