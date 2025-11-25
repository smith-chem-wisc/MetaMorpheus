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
    }
}