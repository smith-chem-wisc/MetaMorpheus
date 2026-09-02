using EngineLayer;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    internal static class TmtExperimentalDesignTests
    {
        [Test]
        public static void WriteTmtExperimentalDesignTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtExperimentalDesignTest");
            Directory.CreateDirectory(outputFolder);

            var file1 = Path.Combine(outputFolder, "file.1.raw");
            var file2 = Path.Combine(outputFolder, "file.2.raw");
            var file3 = Path.Combine(outputFolder, "file.3.raw");

            var anns = new List<TmtPlexAnnotation>
            {
                new TmtPlexAnnotation { Tag="126", SampleName="S1", Condition="C1", BiologicalReplicate=1 }
            };

            // define simple per-file state (same plex and channel set for all)
            var files = new List<TmtFileInfo>
            {
                new TmtFileInfo(file1, "PlexA", 1, 1, anns),
                new TmtFileInfo(file2, "PlexA", 2, 1, anns),
                new TmtFileInfo(file3, "PlexA", 3, 1, anns)
            };

            // write and read back
            _ = TmtExperimentalDesign.Write(files);

            var readFiles = TmtExperimentalDesign.Read(
                Path.Combine(outputFolder, GlobalVariables.TmtExperimentalDesignFileName),
                new List<string> { file1, file2, file3 },
                out var errors);

            Assert.That(!errors.Any(), "No errors expected");
            Assert.That(readFiles.Count == 3);
            Assert.That(readFiles.All(f => f.Plex.Equals("PlexA", StringComparison.OrdinalIgnoreCase)));
            Assert.That(readFiles.All(f => f.Annotations.Count == 1));
            Assert.That(readFiles.All(f => f.Annotations.First().Tag == "126"));

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestTmtExperimentalDesignErrors()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtExperimentalDesignErrors");
            Directory.CreateDirectory(outputFolder);

            var file1 = Path.Combine(outputFolder, "file.1.raw");
            var file2 = Path.Combine(outputFolder, "file.2.raw");
            var designPath = Path.Combine(outputFolder, GlobalVariables.TmtExperimentalDesignFileName);

            // common annotations
            var anns = new List<TmtPlexAnnotation>
            {
                new TmtPlexAnnotation { Tag="126", SampleName="SampleX", Condition="CondX", BiologicalReplicate=1 }
            };

            // 1) Duplicate (Sample, Bio, Fraction, Tech) across files -> error
            {
                var files = new List<TmtFileInfo>
                {
                    new TmtFileInfo(file1, "PlexA", 1, 1, anns),
                    new TmtFileInfo(file2, "PlexA", 1, 1, anns) // same Sample/Bio/Fraction/Tech
                };

                _ = TmtExperimentalDesign.Write(files);
                _ = TmtExperimentalDesign.Read(designPath, new List<string> { file1, file2 }, out var errors);
                Assert.That(errors.Any(), "Duplicate Sample/Bio/Fraction/Tech should error");
            }

            // 2) Conflicting per-file assignments -> error (same file appears with different fraction/tech)
            {
                // write manually to simulate conflict
                var lines = new List<string>
                {
                    TmtExperimentalDesign.Header,
                    $"{file1}\tPlexA\tS1\t126\tC1\t1\t1\t1",
                    $"{file1}\tPlexA\tS1\t127\tC1\t1\t2\t1"
                };
                File.WriteAllLines(designPath, lines);

                _ = TmtExperimentalDesign.Read(designPath, new List<string> { file1 }, out var errors);
                Assert.That(errors.Any(), "Conflicting per-file state should error");
            }

            // 3) Missing files from design -> error
            {
                var files = new List<TmtFileInfo>
                {
                    new TmtFileInfo(file1, "PlexA", 1, 1, anns)
                };

                _ = TmtExperimentalDesign.Write(files);

                _ = TmtExperimentalDesign.Read(designPath, new List<string> { file1, file2 }, out var errors);
                Assert.That(errors.Any(), "Missing file(s) in design should error");
            }

            // 4) Non-integer fields -> error (manually write invalid file)
            {
                var lines = new List<string>
                {
                    TmtExperimentalDesign.Header,
                    // bad Biological Replicate
                    $"{file1}\tPlexA\tS1\t126\tC1\ta\t1\t1"
                };
                File.WriteAllLines(designPath, lines);

                _ = TmtExperimentalDesign.Read(designPath, new List<string> { file1 }, out var errors1);
                Assert.That(errors1.Any(), "Non-integer biological replicate should error");

                // bad Fraction
                lines[1] = $"{file1}\tPlexA\tS1\t126\tC1\t1\ta\t1";
                File.WriteAllLines(designPath, lines);
                _ = TmtExperimentalDesign.Read(designPath, new List<string> { file1 }, out var errors2);
                Assert.That(errors2.Any(), "Non-integer fraction should error");

                // bad Technical Replicate
                lines[1] = $"{file1}\tPlexA\tS1\t126\tC1\t1\t1\ta";
                File.WriteAllLines(designPath, lines);
                _ = TmtExperimentalDesign.Read(designPath, new List<string> { file1 }, out var errors3);
                Assert.That(errors3.Any(), "Non-integer technical replicate should error");
            }

            Directory.Delete(outputFolder, true);
        }

        /// <summary>
        /// A directory of this test class's own, created fresh so a previous run cannot pollute it.
        /// </summary>
        private static string NewFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        /// <summary>
        /// A design whose rows name no file in the run is a design pointed at the wrong data, not an
        /// empty one. It previously returned no errors and an empty annotation set, so the command line
        /// reported that it had read the design successfully while quantification had nothing to use.
        /// </summary>
        [Test]
        public static void DesignThatMatchesNoFileInTheRunIsAnError()
        {
            string folder = NewFolder("TmtDesignNoMatch");
            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);

            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                "someone_elses.raw	Plex1	S1	126	Control	1	1	1	study sample",
                "someone_elses.raw	Plex1	S2	127N	Treated	1	1	1	study sample",
            });

            var files = TmtExperimentalDesign.Read(
                designPath,
                new List<string> { Path.Combine(folder, "actually_searched.raw") },
                out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(files, Is.Empty);
                Assert.That(errors.Any(e => e.Contains("None of the 2 row(s)")), Is.True,
                    "a design that resolves to nothing must say so rather than reading as success");
            });
        }

        /// <summary>
        /// A design written with bare file names sits beside the raw files it names. Resolving those
        /// against the process working directory made the same design file work or match nothing
        /// depending on where MetaMorpheus was launched from; they resolve against the design file's
        /// own directory now.
        /// </summary>
        [Test]
        public static void BareFileNamesResolveAgainstTheDesignFileNotTheWorkingDirectory()
        {
            string folder = NewFolder("TmtDesignRelativePaths");
            string rawPath = Path.Combine(folder, "sample.raw");
            File.WriteAllText(rawPath, string.Empty);   // only its existence matters here

            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                "sample.raw	Plex1	S1	126	Control	1	1	1	study sample",
                "sample.raw	Plex1	S2	127N	Treated	1	1	1	reference",
            });

            // Deliberately not the design's directory -- this is the condition that used to break it.
            string originalWorkingDirectory = Directory.GetCurrentDirectory();
            Directory.SetCurrentDirectory(TestContext.CurrentContext.TestDirectory);
            try
            {
                var files = TmtExperimentalDesign.Read(designPath, new List<string> { rawPath }, out var errors);

                Assert.Multiple(() =>
                {
                    Assert.That(errors, Is.Empty);
                    Assert.That(files, Has.Count.EqualTo(1));
                    Assert.That(files.Single().FullFilePathWithExtension, Is.EqualTo(rawPath));
                    Assert.That(files.Single().Annotations, Has.Count.EqualTo(2));
                });
            }
            finally
            {
                Directory.SetCurrentDirectory(originalWorkingDirectory);
            }
        }

        /// <summary>
        /// A rooted path in the File column is taken as written, so a design that names absolute paths
        /// keeps working regardless of where it or the process happens to sit.
        /// </summary>
        [Test]
        public static void RootedPathsAreTakenAsWritten()
        {
            string folder = NewFolder("TmtDesignRootedPaths");
            string rawPath = Path.Combine(folder, "rooted.raw");
            File.WriteAllText(rawPath, string.Empty);

            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                $"{rawPath}	Plex1	S1	126	Control	1	1	1	study sample",
            });

            var files = TmtExperimentalDesign.Read(designPath, new List<string> { rawPath }, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors, Is.Empty);
                Assert.That(files.Single().FullFilePathWithExtension, Is.EqualTo(rawPath));
            });
        }

        /// <summary>
        /// A provided path that is not already normalized -- here one routed through a '..' segment --
        /// still matches the design row that names the same file. Both the row-level in-this-run test
        /// and the "did not contain the file(s)" check normalize their inputs, so they cannot disagree
        /// and report a file as missing from a design that names it.
        /// </summary>
        [Test]
        public static void AnUnnormalizedProvidedPathStillMatchesItsDesignRow()
        {
            string folder = NewFolder("TmtDesignUnnormalizedPath");
            string rawPath = Path.Combine(folder, "unnormalized.raw");
            File.WriteAllText(rawPath, string.Empty);

            // Same file, reached through a redundant directory hop.
            string detoured = Path.Combine(folder, "sub", "..", "unnormalized.raw");
            Assert.That(detoured, Is.Not.EqualTo(rawPath));

            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                $"{rawPath}	Plex1	S1	126	Control	1	1	1	study sample",
            });

            var files = TmtExperimentalDesign.Read(designPath, new List<string> { detoured }, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors, Is.Empty);
                Assert.That(files.Single().FullFilePathWithExtension, Is.EqualTo(rawPath));
                Assert.That(files.Single().Annotations, Has.Count.EqualTo(1));
            });
        }

        /// <summary>
        /// A row that stops short of the required columns is reported by line number rather than
        /// silently dropped. The count is against the required width, not the full header, so a design
        /// written before the optional Sample Type column existed still loads -- that case is covered
        /// separately; this one is a genuinely truncated row.
        /// </summary>
        [Test]
        public static void TruncatedRowIsReportedByLineNumber()
        {
            string folder = NewFolder("TmtDesignTruncatedRow");
            string rawPath = Path.Combine(folder, "sample.raw");
            File.WriteAllText(rawPath, string.Empty);

            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                "sample.raw	Plex1	S1	126	Control	1	1	1	study sample",
                "sample.raw	Plex1	S2",                                   // stops after three cells
            });

            TmtExperimentalDesign.Read(designPath, new List<string> { rawPath }, out var errors);

            Assert.That(errors.Any(e => e.Contains("Line 3") && e.Contains("fewer columns")), Is.True,
                "the offending line number is the only thing that makes this actionable");
        }

        /// <summary>
        /// A non-numeric Biological Replicate is rejected rather than silently becoming zero, which
        /// would collapse two samples onto one another during roll-up.
        /// </summary>
        [Test]
        public static void NonNumericBiologicalReplicateIsRejected()
        {
            string folder = NewFolder("TmtDesignBadBioRep");
            string rawPath = Path.Combine(folder, "sample.raw");
            File.WriteAllText(rawPath, string.Empty);

            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                "sample.raw	Plex1	S1	126	Control	first	1	1	study sample",
            });

            TmtExperimentalDesign.Read(designPath, new List<string> { rawPath }, out var errors);

            Assert.That(errors.Any(e => e.Contains("Biological Replicate")), Is.True);
        }

        /// <summary>
        /// One file cannot be in two plexes, or at two fractions, at once. The design is per-file for
        /// those three fields, so a row that disagrees with an earlier row about them is a mistake that
        /// would otherwise be resolved by whichever row happened to be read first.
        /// </summary>
        [Test]
        public static void ConflictingPlexOrFractionForOneFileIsReported()
        {
            string folder = NewFolder("TmtDesignConflict");
            string rawPath = Path.Combine(folder, "sample.raw");
            File.WriteAllText(rawPath, string.Empty);

            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                "sample.raw	Plex1	S1	126	Control	1	1	1	study sample",
                "sample.raw	Plex2	S2	127N	Treated	1	2	1	study sample",   // different plex AND fraction
            });

            TmtExperimentalDesign.Read(designPath, new List<string> { rawPath }, out var errors);

            Assert.That(errors.Any(e => e.Contains("inconsistent Plex")), Is.True,
                "the same file cannot carry two plex or fraction assignments");
        }

        [Test]
        public static void MissingDesignFileIsAnError()
        {
            var files = TmtExperimentalDesign.Read(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "no-such-tmt-design.txt"),
                new List<string>(), out var errors);

            Assert.That(files.Count == 0);
            Assert.That(errors.Count == 1);
            Assert.That(errors.Single().Contains("not found"));
        }

        [Test]
        public static void EmptyFileAndBadHeaderAreBothRejected()
        {
            string folder = NewFolder("TmtDesignHeader");
            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);

            // An empty file has no header at all.
            File.WriteAllText(designPath, string.Empty);
            _ = TmtExperimentalDesign.Read(designPath, new List<string>(), out var emptyErrors);
            Assert.That(emptyErrors.Any(e => e.Contains("header is invalid")), "Empty file should be rejected");

            // A header missing required columns.
            File.WriteAllLines(designPath, new[] { "File\tPlex\tSample Name" });
            _ = TmtExperimentalDesign.Read(designPath, new List<string>(), out var shortErrors);
            Assert.That(shortErrors.Any(e => e.Contains("header is invalid")), "Incomplete header should be rejected");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The header is matched by column NAME, not by position, so a design written by hand with its
        /// columns in a different order still parses. That is what the per-column index lookup buys.
        /// </summary>
        [Test]
        public static void HeaderColumnsMayBeReordered()
        {
            string folder = NewFolder("TmtDesignReordered");
            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            string file1 = Path.Combine(folder, "file.1.raw");

            File.WriteAllLines(designPath, new[]
            {
                "Technical Replicate\tFraction\tBiological Replicate\tCondition\tTMT Channel\tSample Name\tPlex\tFile",
                $"1\t2\t3\tCond\t127N\tSampleA\tPlexB\t{file1}"
            });

            var files = TmtExperimentalDesign.Read(designPath, new List<string> { file1 }, out var errors);

            Assert.That(!errors.Any(), "A reordered but complete header should parse");
            Assert.That(files.Count == 1);
            Assert.That(files[0].Plex == "PlexB");
            Assert.That(files[0].Fraction == 2);
            Assert.That(files[0].TechnicalReplicate == 1);
            Assert.That(files[0].Annotations.Single().Tag == "127N");
            Assert.That(files[0].Annotations.Single().SampleName == "SampleA");
            Assert.That(files[0].Annotations.Single().BiologicalReplicate == 3);

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Blank lines are skipped, a row with too few columns is reported by line number, and a row
        /// naming no file is ignored rather than reported -- the last because a tab-only trailing row
        /// is what a spreadsheet export produces, and it is not worth stopping a run for.
        /// </summary>
        [Test]
        public static void MalformedRowsAreSkippedOrReportedByLineNumber()
        {
            string folder = NewFolder("TmtDesignMalformedRows");
            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            string file1 = Path.Combine(folder, "file.1.raw");

            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                $"{file1}\tPlexA\tS1\t126\tC1\t1\t1\t1",
                "",                                 // blank -> skipped silently
                $"{file1}\tPlexA\tS1",              // too few columns -> reported
                "\tPlexA\tS2\t127\tC1\t1\t1\t1"     // names no file -> ignored
            });

            var files = TmtExperimentalDesign.Read(designPath, new List<string> { file1 }, out var errors);

            Assert.That(files.Count == 1, "Only the one well-formed file should be returned");
            Assert.That(errors.Count == 1, "Exactly one error: the short row");
            Assert.That(errors.Single().Contains("Line 4"), "Errors are reported by 1-based line number");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Fraction and technical replicate are 1-based in the file. A zero is what a caller gets by
        /// forwarding a 0-based index, so it is rejected rather than accepted as "unset".
        /// </summary>
        [Test]
        public static void FractionAndTechnicalReplicateMustBeOneBased()
        {
            string folder = NewFolder("TmtDesignOneBased");
            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            string file1 = Path.Combine(folder, "file.1.raw");

            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                $"{file1}\tPlexA\tS1\t126\tC1\t1\t0\t1"
            });
            _ = TmtExperimentalDesign.Read(designPath, new List<string>(), out var fractionErrors);
            Assert.That(fractionErrors.Any(e => e.Contains("Fraction must be >= 1")));

            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                $"{file1}\tPlexA\tS1\t126\tC1\t1\t1\t0"
            });
            _ = TmtExperimentalDesign.Read(designPath, new List<string>(), out var techErrors);
            Assert.That(techErrors.Any(e => e.Contains("Technical Replicate must be >= 1")));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// An empty spectra-file list means "describe whatever the design names", so every row is kept
        /// and the file-not-in-design check does not run. This is the path the GUI seeds from.
        /// </summary>
        [Test]
        public static void AnEmptyFileListAcceptsEveryRowInTheDesign()
        {
            string folder = NewFolder("TmtDesignNoFileFilter");
            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            string file1 = Path.Combine(folder, "file.1.raw");
            string file2 = Path.Combine(folder, "file.2.raw");

            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                $"{file1}\tPlexA\tS1\t126\tC1\t1\t1\t1",
                $"{file2}\tPlexA\tS2\t127\tC2\t2\t1\t1"
            });

            var files = TmtExperimentalDesign.Read(designPath, new List<string>(), out var errors);

            Assert.That(!errors.Any(), "With no file list, nothing can be missing from the design");
            Assert.That(files.Count == 2);
            // Both files share a plex, so both carry that plex's whole channel set.
            Assert.That(files.All(f => f.Annotations.Count == 2));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The (sample, biorep, fraction, techrep) uniqueness rule applies only to rows that name a
        /// sample. Two unnamed channels in one plex are not duplicates of each other.
        /// </summary>
        [Test]
        public static void RowsWithNoSampleNameAreExemptFromTheUniquenessRule()
        {
            string folder = NewFolder("TmtDesignBlankSample");
            string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            string file1 = Path.Combine(folder, "file.1.raw");

            File.WriteAllLines(designPath, new[]
            {
                TmtExperimentalDesign.Header,
                $"{file1}\tPlexA\t\t126\tC1\t1\t1\t1",
                $"{file1}\tPlexA\t\t127\tC1\t1\t1\t1"
            });

            var files = TmtExperimentalDesign.Read(designPath, new List<string> { file1 }, out var errors);

            Assert.That(!errors.Any(), "Blank sample names are not duplicates");
            Assert.That(files.Single().Annotations.Count == 2);

            Directory.Delete(folder, true);
        }

        [Test]
        public static void WriteRefusesAnEmptyDesign()
        {
            Assert.Throws<InvalidOperationException>(() => TmtExperimentalDesign.Write(null));
            Assert.Throws<InvalidOperationException>(() => TmtExperimentalDesign.Write(new List<TmtFileInfo>()));
        }

        /// <summary>
        /// A file with no channel annotations still gets a row, so a design the user is part-way
        /// through authoring does not silently lose the file. Every per-channel cell on that row is
        /// blank — including the trailing Sample Type cell, hence the tab at the end.
        /// </summary>
        [Test]
        public static void WriteEmitsAPlaceholderRowForAFileWithNoAnnotations()
        {
            string folder = NewFolder("TmtDesignNoAnnotations");
            string file1 = Path.Combine(folder, "file.1.raw");

            string written = TmtExperimentalDesign.Write(
                new List<TmtFileInfo> { new TmtFileInfo(file1, "PlexA", 2, 3, null) });

            var lines = File.ReadAllLines(written);
            Assert.That(lines.Length == 2);
            Assert.That(lines[0] == TmtExperimentalDesign.Header);
            Assert.That(lines[1] == $"{file1}\tPlexA\t\t\t\t\t2\t3\t");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Rows are written fraction-major, so a design round-trips in a stable order rather than in
        /// whatever order the GUI grid happened to hold.
        /// </summary>
        [Test]
        public static void WriteOrdersRowsByFractionThenTechnicalReplicate()
        {
            string folder = NewFolder("TmtDesignWriteOrder");
            string file1 = Path.Combine(folder, "file.1.raw");
            string file2 = Path.Combine(folder, "file.2.raw");
            string file3 = Path.Combine(folder, "file.3.raw");
            var anns = new List<TmtPlexAnnotation>
            {
                new TmtPlexAnnotation { Tag = "126", SampleName = "S1", Condition = "C1", BiologicalReplicate = 1 }
            };

            string written = TmtExperimentalDesign.Write(new List<TmtFileInfo>
            {
                new TmtFileInfo(file2, "PlexA", 2, 1, anns),
                new TmtFileInfo(file3, "PlexA", 1, 2, anns),
                new TmtFileInfo(file1, "PlexA", 1, 1, anns)
            });

            var rows = File.ReadAllLines(written).Skip(1).Select(l => l.Split('\t')[0]).ToList();
            Assert.That(rows.SequenceEqual(new[] { file1, file3, file2 }));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Two TmtFileInfo naming the same file are the same file, whatever the casing -- and their
        /// hash codes agree, so they collapse in a set rather than both surviving it.
        /// </summary>
        [Test]
        public static void TmtFileInfoIdentityIsTheFilePath()
        {
            var a = new TmtFileInfo(@"C:\data\file.1.raw", "PlexA", 1, 1, null);
            var b = new TmtFileInfo(@"c:\DATA\FILE.1.RAW", "PlexB", 9, 9, null);
            var c = new TmtFileInfo(@"C:\data\file.2.raw", "PlexA", 1, 1, null);

            Assert.That(a.Equals(b), "Same path, different case -> the same file");
            Assert.That(a.GetHashCode() == b.GetHashCode(), "Equal instances must hash equally");
            Assert.That(!a.Equals(c));
            Assert.That(!a.Equals(null));
            Assert.That(!a.Equals("not a TmtFileInfo"));
            Assert.That(new HashSet<TmtFileInfo> { a, b, c }.Count == 2);
        }

        [Test]
        public static void TmtFileInfoNormalizesItsOptionalArguments()
        {
            var info = new TmtFileInfo(Path.Combine("some", "dir", "file.1.raw"), null, 1, 1, null);

            Assert.That(info.Plex == string.Empty, "A null plex is stored as empty, never null");
            Assert.That(info.Annotations.Count == 0, "Null annotations become an empty list");
            Assert.That(info.ToString() == "file.1.raw", "ToString is the bare file name");
        }
        /// <summary>
        /// A row is judged one cell at a time and reported once. Every validation guard in the row
        /// loop ends in `continue`, so the first bad cell is the only one named -- the row is not
        /// re-judged on the columns after it.
        ///
        /// This is a user-facing contract, not an internal one: Read's error list is what
        /// ResolveExperimentalDesign prints to the console. Without the skip, one typo in a Sample
        /// Type cell also reports an invalid Biological Replicate, because TryParse left the out
        /// parameter at its default -- an error naming a cell the user filled in correctly.
        ///
        /// Each case is asserted by what the errors must NOT say, since a count cannot distinguish
        /// them: the "did not contain the file(s)" error appears either way, because a row skipped
        /// for being invalid and a row never reached both leave the file undefined.
        /// </summary>
        [Test]
        public static void ARowIsReportedOnceAtItsFirstBadCell()
        {
            List<string> ErrorsForRow(string name, string row)
            {
                string folder = NewFolder(name);
                string rawPath = Path.Combine(folder, "sample.raw");
                File.WriteAllText(rawPath, string.Empty);

                string designPath = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
                File.WriteAllLines(designPath, new[] { TmtExperimentalDesign.Header, row });

                TmtExperimentalDesign.Read(designPath, new List<string> { rawPath }, out var errors);
                return errors;
            }

            // Sample Type is unrecognised, and the Biological Replicate after it is also unparseable.
            var badType = ErrorsForRow("TmtDesignFirstBadCellType",
                "sample.raw\tPlex1\tS1\t126\tControl\tnot-a-number\t1\t1\tnot-a-sample-type");

            // Biological Replicate is unparseable, and the Fraction after it is below one.
            var badBio = ErrorsForRow("TmtDesignFirstBadCellBio",
                "sample.raw\tPlex1\tS1\t126\tControl\tnot-a-number\t0\t1\tstudy sample");

            // Fraction is below one, and the Technical Replicate after it is too.
            var badFraction = ErrorsForRow("TmtDesignFirstBadCellFraction",
                "sample.raw\tPlex1\tS1\t126\tControl\t1\t0\t0\tstudy sample");

            Assert.That(badType.Any(e => e.Contains("Sample Type")), "the unrecognised Sample Type is reported");
            Assert.That(!badType.Any(e => e.Contains("Biological Replicate")),
                "the row was judged past its Sample Type, reporting a Biological Replicate the guard never parsed");

            Assert.That(badBio.Any(e => e.Contains("Biological Replicate")), "the invalid Biological Replicate is reported");
            Assert.That(!badBio.Any(e => e.Contains("Fraction")),
                "the row was judged past its Biological Replicate, reporting a Fraction as well");

            Assert.That(badFraction.Any(e => e.Contains("Fraction")), "the out-of-range Fraction is reported");
            Assert.That(!badFraction.Any(e => e.Contains("Technical Replicate")),
                "the row was judged past its Fraction, reporting a Technical Replicate as well");
        }

        #region Review findings (#2774)

        private static string WriteDesign(string name, params string[] rows)
        {
            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtDesignReview", name);
            Directory.CreateDirectory(dir);
            string path = Path.Combine(dir, "TmtDesign.txt");
            var lines = new List<string>
            {
                "File\tPlex\tSample Name\tTMT Channel\tCondition\tBiological Replicate\tFraction\tTechnical Replicate\tSample Type"
            };
            lines.AddRange(rows);
            File.WriteAllLines(path, lines);
            return path;
        }

        /// <summary>
        /// Two rows describing one channel differently must be reported, not silently collapsed.
        ///
        /// Fractions of a plex legitimately repeat the same channel-to-sample map, which is why a
        /// repeat is collapsed rather than refused. But collapsing rows that DISAGREE keeps the first
        /// and discards the second with nothing said -- and the uniqueness rule cannot catch it, since
        /// two different sample names are two different tuples, so it only fires when the duplicates
        /// agree, which is the harmless case.
        /// </summary>
        [Test]
        public static void ChannelDescribedTwiceWithDifferentSamplesIsReported()
        {
            string path = WriteDesign("disagree",
                "sample.raw\tPlex1\tTumour\t126\tDisease\t1\t1\t1\tstudy sample",
                "sample.raw\tPlex1\tNormal\t126\tControl\t2\t1\t1\tstudy sample");

            TmtExperimentalDesign.Read(path, new List<string>(), out var errors);

            Assert.That(errors, Is.Not.Empty, "the second description was discarded silently");
            Assert.That(string.Join(" ", errors), Does.Contain("126").And.Contain("disagree"));
        }

        /// <summary>
        /// A bridge channel is by definition the same material carried in more than one plex under one
        /// name -- which is what the Sample Type column added by this PR exists to describe. Collecting
        /// the uniqueness tuples across the whole file rejected exactly that design; the only way
        /// through was renaming the bridges apart, which throws away the fact that they are the same
        /// material and so defeats the point of having them.
        /// </summary>
        [Test]
        public static void ABridgeSharedAcrossPlexesIsNotADuplicate()
        {
            string path = WriteDesign("bridge",
                "p1.raw\tPlex1\tTumour\t126\tDisease\t1\t1\t1\tstudy sample",
                "p1.raw\tPlex1\tBridge\t127N\tBridge\t1\t1\t1\tbridge",
                "p2.raw\tPlex2\tNormal\t126\tControl\t2\t1\t1\tstudy sample",
                "p2.raw\tPlex2\tBridge\t127N\tBridge\t1\t1\t1\tbridge");

            var files = TmtExperimentalDesign.Read(path, new List<string>(), out var errors);

            Assert.That(errors, Is.Empty, string.Join(" | ", errors));
            Assert.That(files, Has.Count.EqualTo(2));
        }

        /// <summary>
        /// The same rule still fires within one plex, which is the case worth catching.
        /// </summary>
        [Test]
        public static void TheSameSampleTwiceInOnePlexIsStillADuplicate()
        {
            string path = WriteDesign("dupInPlex",
                "p1.raw\tPlex1\tTumour\t126\tDisease\t1\t1\t1\tstudy sample",
                "p1.raw\tPlex1\tTumour\t127N\tDisease\t1\t1\t1\tstudy sample");

            TmtExperimentalDesign.Read(path, new List<string>(), out var errors);

            Assert.That(string.Join(" ", errors), Does.Contain("Duplicate combination"));
        }

        /// <summary>
        /// Biological Replicate is validated the way Fraction and Technical Replicate are. A bare
        /// int.TryParse accepted 0 and -3, and ToMzLibDesign passes the value straight through.
        /// </summary>
        [TestCase("0")]
        [TestCase("-3")]
        public static void BiologicalReplicateBelowOneIsRejected(string bio)
        {
            string path = WriteDesign("bio" + bio.Replace("-", "neg"),
                $"sample.raw\tPlex1\tTumour\t126\tDisease\t{bio}\t1\t1\tstudy sample");

            TmtExperimentalDesign.Read(path, new List<string>(), out var errors);

            Assert.That(string.Join(" ", errors), Does.Contain("Biological Replicate must be >= 1"));
        }

        /// <summary>
        /// Null means "do not filter", the same meaning the guard in ToMzLibDesign gives it. It used to
        /// throw ArgumentNullException, so the two ends of the same parameter disagreed.
        /// </summary>
        [Test]
        public static void ANullFileListMeansDoNotFilterRatherThanThrowing()
        {
            string path = WriteDesign("nullFiles",
                "sample.raw\tPlex1\tTumour\t126\tDisease\t1\t1\t1\tstudy sample");

            List<TmtFileInfo> files = null;
            Assert.That(() => files = TmtExperimentalDesign.Read(path, null, out _), Throws.Nothing);
            Assert.That(files, Has.Count.EqualTo(1));
        }

        /// <summary>
        /// The placeholder row Write emits for a file with no annotations has to survive being read
        /// back, since that is the whole reason it is written -- so a design someone is part-way
        /// through authoring does not lose the file. It did not: the blank Biological Replicate failed
        /// int.TryParse, the row was skipped, and the file was then reported missing from the design.
        ///
        /// This matters most for the GUI PR, where a part-authored design is the normal state rather
        /// than an edge case.
        /// </summary>
        [Test]
        public static void APlaceholderRowForAnUnannotatedFileSurvivesARoundTrip()
        {
            string path = WriteDesign("placeholder",
                "sample.raw\tPlex1\t\t\t\t\t1\t1\t");

            var files = TmtExperimentalDesign.Read(path, new List<string>(), out var errors);

            Assert.That(errors, Is.Empty, string.Join(" | ", errors));
            Assert.That(files, Has.Count.EqualTo(1), "the file the placeholder exists to preserve was lost");
            Assert.That(files[0].Annotations, Is.Empty, "and it carries no annotations, which is the point");
        }

        #endregion

    }
}
