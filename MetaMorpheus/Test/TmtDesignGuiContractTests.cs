using EngineLayer;
using MetaMorpheusCommandLine;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The contract the TMT experimental design GUI has to keep with the command line: a design the
    /// window reports as saved must be one <see cref="TmtExperimentalDesign.Read"/> loads without
    /// error, <see cref="TmtExperimentalDesign.ToMzLibDesign"/> projects, and
    /// <see cref="Program.ResolveExperimentalDesign"/> accepts.
    /// </summary>
    /// <remarks>
    /// The window's own code is in GUI.csproj and Test.csproj does not reference it, so the WPF
    /// handlers cannot be driven from here. What these tests pin is the model-level shape the
    /// window's Save produces -- one <see cref="TmtFileInfo"/> per file, grouped by plex, sharing one
    /// annotation list per plex, ordered by fraction then technical replicate -- and the tag-table
    /// facts the window's channel grid depends on.
    /// </remarks>
    [TestFixture]
    internal static class TmtDesignGuiContractTests
    {
        private const int Continue = 0;

        /// <summary>
        /// The premise of the stale-label bug in the annotate window: rebuilding the channel rows only
        /// when the channel COUNT changes is not enough, because two tag types share a count.
        /// </summary>
        [Test]
        public static void Itraq4AndDiLeu4ShareAChannelCountButNotTheirLabels()
        {
            var itraq4 = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.iTRAQ4);
            var dileu4 = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.diLeu4);

            Assert.Multiple(() =>
            {
                Assert.That(itraq4.Count, Is.EqualTo(dileu4.Count),
                    "the two tags have the same channel count, which is what made a count check look sufficient");
                Assert.That(itraq4, Is.Not.EqualTo(dileu4),
                    "and different labels, so carrying rows across on a count match keeps the wrong ones");
                Assert.That(itraq4, Is.EqualTo(new[] { "114", "115", "116", "117" }));
                Assert.That(dileu4, Is.EqualTo(new[] { "115", "116", "117", "118" }));
            });
        }

        /// <summary>
        /// iTRAQ4/diLeu4 is the ONLY pair that collides on count, so a test that switches between
        /// those two covers the whole hazard. If a tag type is ever added that shares a count with
        /// another, this fails and says so.
        /// </summary>
        [Test]
        public static void FourIsTheOnlyChannelCountTwoTagTypesShare()
        {
            var countsByType = Enum.GetValues(typeof(IsobaricMassTagType))
                .Cast<IsobaricMassTagType>()
                .ToDictionary(t => t, t => IsobaricMassTag.GetReporterIonLabels(t)?.Count ?? 0);

            var shared = countsByType.GroupBy(kv => kv.Value)
                                     .Where(g => g.Count() > 1)
                                     .ToList();

            Assert.That(shared.Count, Is.EqualTo(1),
                "expected exactly one colliding channel count, found: "
                + string.Join("; ", shared.Select(g => g.Key + " -> " + string.Join(", ", g.Select(kv => kv.Key)))));
            Assert.That(shared.Single().Key, Is.EqualTo(4));
        }

        /// <summary>
        /// A tag carried over from a different tag type is not silently accepted downstream -- this is
        /// the failure the window would have produced, and it is loud rather than silent.
        /// </summary>
        [Test]
        public static void AnItraq4LabelInADiLeu4PlexIsRejectedByTheProjection()
        {
            var diLeu4 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.diLeu4);

            // "114" is iTRAQ4's first channel and is not a diLeu4 channel at all.
            var file = FileWith(SpectraPath("stale.mzML"), "Plex1", 1, 1, ("114", 1, TmtSampleType.StudySample));

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { file }, diLeu4, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(design, Is.Null);
                Assert.That(errors.Any(), Is.True, "a channel that is not part of the tag has to be reported");
            });
        }

        /// <summary>
        /// The whole two-plex design the design window builds on Save, written and then read back the
        /// way the window now reads it back and the way a command-line run reads it.
        /// </summary>
        [Test]
        public static void AWholeTwoPlexDesignRoundTripsThroughTheParser()
        {
            string folder = NewFolder("TmtGuiRoundTrip");
            string plex1File = Path.Combine(folder, "plex1.mzML");
            string plex2File = Path.Combine(folder, "plex2.mzML");

            var design = BuildDesignTheWayTheGuiDoes(
                (plex1File, "Plex1", 1, 1),
                (plex2File, "Plex2", 1, 1));

            string written = TmtExperimentalDesign.Write(design);

            var readBack = TmtExperimentalDesign.Read(written,
                new List<string> { plex1File, plex2File }, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors, Is.Empty, "a design the window would report as saved must load: "
                    + string.Join(" | ", errors));
                Assert.That(readBack.Count, Is.EqualTo(2), "both plexes' files survive the write");
                Assert.That(readBack.Select(f => f.Plex).OrderBy(p => p),
                    Is.EqualTo(new[] { "Plex1", "Plex2" }),
                    "annotating the second plex used to truncate the first out of the file");
            });

            var projected = TmtExperimentalDesign.ToMzLibDesign(readBack,
                IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT11), out var projectionErrors);

            Assert.Multiple(() =>
            {
                Assert.That(projectionErrors, Is.Empty);
                Assert.That(projected.FileNameSampleInfoDictionary.Count, Is.EqualTo(2));
                Assert.That(projected.FileNameSampleInfoDictionary.Values.All(v => v.Length == 11), Is.True,
                    "every file carries all eleven channels");
            });

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The same design put in front of the command line, which is the consumer that decides
        /// whether the file was worth writing.
        /// </summary>
        [Test]
        public static void AWholeTwoPlexDesignIsAcceptedByTheCommandLine()
        {
            string folder = NewFolder("TmtGuiCmdAccepts");
            string plex1File = Path.Combine(folder, "plex1.mzML");
            string plex2File = Path.Combine(folder, "plex2.mzML");

            TmtExperimentalDesign.Write(BuildDesignTheWayTheGuiDoes(
                (plex1File, "Plex1", 1, 1),
                (plex2File, "Plex2", 1, 1)));

            var written = new List<string>();
            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { plex1File, plex2File },
                normalizationRequested: true, reportToConsole: true, write: written.Add);

            Assert.Multiple(() =>
            {
                Assert.That(code, Is.EqualTo(Continue), "the command line refused a design the GUI wrote");
                Assert.That(written.Single(),
                    Is.EqualTo("Read " + GlobalVariables.TmtExperimentalDesignFileName + " successfully"));
            });

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A plex the user has not annotated yet is still written, as the file-only placeholder row,
        /// so saving a part-authored design does not lose its files.
        /// </summary>
        [Test]
        public static void APlexWithNoAnnotationsKeepsItsFiles()
        {
            string folder = NewFolder("TmtGuiPartAuthored");
            string annotated = Path.Combine(folder, "annotated.mzML");
            string bare = Path.Combine(folder, "bare.mzML");

            var design = new List<TmtFileInfo>
            {
                FileWith(annotated, "Plex1", 1, 1, ("126", 1, TmtSampleType.StudySample)),
                new TmtFileInfo(bare, "Plex2", 1, 1, Array.Empty<TmtPlexAnnotation>())
            };

            string written = TmtExperimentalDesign.Write(design);
            var readBack = TmtExperimentalDesign.Read(written,
                new List<string> { annotated, bare }, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors, Is.Empty);
                Assert.That(readBack.Count, Is.EqualTo(2), "the unannotated plex's file must survive");
                Assert.That(readBack.Single(f => f.Plex == "Plex2").Annotations, Is.Empty);
            });

            Directory.Delete(folder, true);
        }

        #region Helpers

        /// <summary>
        /// One TmtFileInfo per file, grouped by plex, every file in a plex sharing that plex's single
        /// annotation list, ordered by fraction then technical replicate. This mirrors the design
        /// window's BuildDesign.
        /// </summary>
        private static List<TmtFileInfo> BuildDesignTheWayTheGuiDoes(
            params (string path, string plex, int fraction, int techrep)[] rows)
        {
            var labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.TMT11);
            var files = new List<TmtFileInfo>();

            foreach (var group in rows.GroupBy(r => r.plex, StringComparer.OrdinalIgnoreCase)
                                      .OrderBy(g => g.Key, StringComparer.OrdinalIgnoreCase))
            {
                var annotations = labels
                    .Select((label, i) => new TmtPlexAnnotation
                    {
                        Tag = label,
                        SampleName = group.Key + "_S" + (i + 1),
                        Condition = "Control",
                        BiologicalReplicate = 1,
                        SampleType = TmtSampleType.StudySample
                    })
                    .ToList();

                foreach (var r in group.OrderBy(r => r.fraction).ThenBy(r => r.techrep))
                    files.Add(new TmtFileInfo(r.path, group.Key, r.fraction, r.techrep, annotations));
            }

            return files;
        }

        private static TmtFileInfo FileWith(string path, string plex, int fraction, int techrep,
            params (string tag, int bio, TmtSampleType type)[] channels)
        {
            var annotations = channels
                .Select(c => new TmtPlexAnnotation
                {
                    Tag = c.tag,
                    SampleName = "S" + c.bio,
                    Condition = "Control",
                    BiologicalReplicate = c.bio,
                    SampleType = c.type
                })
                .ToList();

            return new TmtFileInfo(path, plex, fraction, techrep, annotations);
        }

        private static string SpectraPath(string name) =>
            Path.DirectorySeparatorChar + Path.Combine("data", name);

        private static string NewFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder))
                Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        #endregion
    }
}
