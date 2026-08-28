using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    /// <summary>
    /// Covers the Sample Type column of TmtDesign.txt and the projection of a TMT design onto
    /// mzLib's <see cref="IExperimentalDesign"/>, which is what the quantification engine consumes.
    /// </summary>
    [TestFixture]
    internal static class TmtSampleTypeAndMzLibDesignTests
    {
        #region Sample Type parsing

        [TestCase("", TmtSampleType.StudySample)]
        [TestCase("   ", TmtSampleType.StudySample)]
        [TestCase("study sample", TmtSampleType.StudySample)]
        [TestCase("Study Sample", TmtSampleType.StudySample)]
        [TestCase("reference", TmtSampleType.Reference)]
        [TestCase("REFERENCE", TmtSampleType.Reference)]
        [TestCase("pooled", TmtSampleType.Reference)]
        [TestCase("bridge", TmtSampleType.Bridge)]
        [TestCase("carrier", TmtSampleType.Carrier)]
        [TestCase("empty", TmtSampleType.Empty)]
        public static void TryParseSampleType_AcceptsTheSdrfSpellings(string value, TmtSampleType expected)
        {
            Assert.IsTrue(TmtExperimentalDesign.TryParseSampleType(value, out var parsed));
            Assert.AreEqual(expected, parsed);
        }

        [TestCase("nonsense")]
        [TestCase("ref")]
        public static void TryParseSampleType_RejectsAnythingElse(string value)
        {
            Assert.IsFalse(TmtExperimentalDesign.TryParseSampleType(value, out _));
        }

        [Test]
        public static void SampleType_RoundTripsThroughItsDesignFileSpelling()
        {
            foreach (TmtSampleType type in Enum.GetValues(typeof(TmtSampleType)))
            {
                string written = TmtExperimentalDesign.ToDesignFileValue(type);
                Assert.IsTrue(TmtExperimentalDesign.TryParseSampleType(written, out var readBack), written);
                Assert.AreEqual(type, readBack);
                Assert.Contains(written, TmtExperimentalDesign.SampleTypeNames.ToList());
            }
        }

        /// <summary>
        /// Only a reference or a bridge normalizes the other channels, and that is the single bit
        /// mzLib's IsobaricQuantSampleInfo carries.
        /// </summary>
        [TestCase(TmtSampleType.Reference, true)]
        [TestCase(TmtSampleType.Bridge, true)]
        [TestCase(TmtSampleType.StudySample, false)]
        [TestCase(TmtSampleType.Carrier, false)]
        [TestCase(TmtSampleType.Empty, false)]
        public static void IsReferenceChannel_IsTrueOnlyForReferenceAndBridge(TmtSampleType type, bool expected)
        {
            Assert.AreEqual(expected, TmtExperimentalDesign.IsReferenceChannel(type));
        }

        #endregion

        #region Sample Type in the design file

        private static string WriteDesign(string folder, string header, params string[] rows)
        {
            Directory.CreateDirectory(folder);
            string path = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(path, new[] { header }.Concat(rows));
            return path;
        }

        /// <summary>
        /// A design file written before the Sample Type column existed must still load, with every
        /// channel reading as a study sample.
        /// </summary>
        [Test]
        public static void Read_DesignWithoutSampleTypeColumn_StillLoads()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtSampleType_Legacy");
            string file = Path.Combine(folder, "run1.raw");

            string legacyHeader = "File\tPlex\tSample Name\tTMT Channel\tCondition\tBiological Replicate\tFraction\tTechnical Replicate";
            string path = WriteDesign(folder, legacyHeader,
                $"{file}\tPlexA\tS1\t126\tControl\t1\t1\t1",
                $"{file}\tPlexA\tS2\t127N\tTreated\t1\t1\t1");

            var files = TmtExperimentalDesign.Read(path, new List<string> { file }, out var errors);

            Assert.IsEmpty(errors);
            Assert.AreEqual(1, files.Count);
            Assert.IsTrue(files[0].Annotations.All(a => a.SampleType == TmtSampleType.StudySample));

            Directory.Delete(folder, true);
        }

        [Test]
        public static void Read_ParsesSampleTypeWhenPresent()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtSampleType_Present");
            string file = Path.Combine(folder, "run1.raw");

            string path = WriteDesign(folder, TmtExperimentalDesign.Header,
                $"{file}\tPlexA\tPool\t126\tReference\t1\t1\t1\treference",
                $"{file}\tPlexA\tS1\t127N\tTreated\t1\t1\t1\tstudy sample",
                $"{file}\tPlexA\t\t127C\t\t1\t1\t1\tempty");

            var files = TmtExperimentalDesign.Read(path, new List<string> { file }, out var errors);

            Assert.IsEmpty(errors);
            var byTag = files[0].Annotations.ToDictionary(a => a.Tag);
            Assert.AreEqual(TmtSampleType.Reference, byTag["126"].SampleType);
            Assert.AreEqual(TmtSampleType.StudySample, byTag["127N"].SampleType);
            Assert.AreEqual(TmtSampleType.Empty, byTag["127C"].SampleType);

            Directory.Delete(folder, true);
        }

        [Test]
        public static void Read_UnrecognisedSampleType_IsReportedAndTheRowIsSkipped()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtSampleType_Bad");
            string file = Path.Combine(folder, "run1.raw");

            string path = WriteDesign(folder, TmtExperimentalDesign.Header,
                $"{file}\tPlexA\tS1\t126\tControl\t1\t1\t1\tbanana");

            TmtExperimentalDesign.Read(path, new List<string> { file }, out var errors);

            Assert.IsTrue(errors.Any(e => e.Contains("not a recognised Sample Type")), string.Join(" | ", errors));

            Directory.Delete(folder, true);
        }

        /// <summary>The header this class writes must be one the same class accepts.</summary>
        [Test]
        public static void WrittenHeader_IsAcceptedOnRead()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtSampleType_RoundTrip");
            string file = Path.Combine(folder, "run1.raw");
            Directory.CreateDirectory(folder);

            var annotations = new List<TmtPlexAnnotation>
            {
                new() { Tag = "126", SampleName = "Pool", Condition = "Reference", BiologicalReplicate = 1, SampleType = TmtSampleType.Reference },
                new() { Tag = "127N", SampleName = "S1", Condition = "Treated", BiologicalReplicate = 1, SampleType = TmtSampleType.StudySample }
            };

            TmtExperimentalDesign.Write(new List<TmtFileInfo> { new(file, "PlexA", 1, 1, annotations) });

            var readBack = TmtExperimentalDesign.Read(
                Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName),
                new List<string> { file }, out var errors);

            Assert.IsEmpty(errors);
            var byTag = readBack[0].Annotations.ToDictionary(a => a.Tag);
            Assert.AreEqual(TmtSampleType.Reference, byTag["126"].SampleType);
            Assert.AreEqual(TmtSampleType.StudySample, byTag["127N"].SampleType);

            Directory.Delete(folder, true);
        }

        #endregion

        #region Projection onto mzLib

        private static IsobaricMassTag Tmt10() => IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

        private static TmtFileInfo FileWith(string path, string plex, int fraction, int techrep,
            params (string tag, string condition, int bio, TmtSampleType type)[] channels)
        {
            var annotations = channels
                .Select(c => new TmtPlexAnnotation
                {
                    Tag = c.tag,
                    SampleName = c.condition + "_" + c.bio,
                    Condition = c.condition,
                    BiologicalReplicate = c.bio,
                    SampleType = c.type
                })
                .ToList();
            return new TmtFileInfo(path, plex, fraction, techrep, annotations);
        }

        /// <summary>
        /// The load-bearing contract: the ISampleInfo array must be in ascending reporter m/z order,
        /// because MetaMorpheus fills ISpectralMatch.Intensities in that order and mzLib maps the two
        /// positionally. A design file listing channels in any other order must not change it.
        /// </summary>
        [Test]
        public static void ToMzLibDesign_OrdersChannelsByReporterMz_NotByDesignFileOrder()
        {
            var tag = Tmt10();
            Assert.IsNotNull(tag, "TMT10 must resolve from the loaded modifications");

            // Deliberately scrambled relative to the tag's own channel order.
            var file = FileWith(@"C:\data\run1.raw", "PlexA", 1, 1,
                ("129N", "D", 1, TmtSampleType.StudySample),
                ("126", "A", 1, TmtSampleType.Reference),
                ("127C", "C", 1, TmtSampleType.StudySample),
                ("127N", "B", 1, TmtSampleType.StudySample));

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { file }, tag, out var errors);

            Assert.IsEmpty(errors);
            var samples = design.FileNameSampleInfoDictionary["run1.raw"];

            var expectedLabels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.TMT10);
            Assert.AreEqual(expectedLabels.Count, samples.Length,
                "one entry per tag channel, so the array lines up with Intensities");

            for (int i = 0; i < expectedLabels.Count; i++)
            {
                var isobaric = (IsobaricQuantSampleInfo)samples[i];
                Assert.AreEqual(expectedLabels[i], isobaric.ChannelLabel);
                Assert.AreEqual(tag.ReporterIonMzs[i], isobaric.ReporterIonMz, 1e-6);
            }

            // Reporter m/z must be strictly increasing, which is what "tag order" means.
            for (int i = 1; i < samples.Length; i++)
            {
                Assert.Less(((IsobaricQuantSampleInfo)samples[i - 1]).ReporterIonMz,
                            ((IsobaricQuantSampleInfo)samples[i]).ReporterIonMz);
            }
        }

        [Test]
        public static void ToMzLibDesign_CarriesTheAnnotationOntoTheMatchingChannel()
        {
            var tag = Tmt10();
            var file = FileWith(@"C:\data\run1.raw", "PlexA", 2, 3,
                ("126", "Reference", 1, TmtSampleType.Reference),
                ("127N", "Treated", 2, TmtSampleType.StudySample));

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { file }, tag, out var errors);
            Assert.IsEmpty(errors);

            var samples = design.FileNameSampleInfoDictionary["run1.raw"]
                .Cast<IsobaricQuantSampleInfo>()
                .ToDictionary(s => s.ChannelLabel);

            Assert.AreEqual("Reference", samples["126"].Condition);
            Assert.AreEqual(1, samples["126"].BiologicalReplicate);
            Assert.IsTrue(samples["126"].IsReferenceChannel);

            Assert.AreEqual("Treated", samples["127N"].Condition);
            Assert.AreEqual(2, samples["127N"].BiologicalReplicate);
            Assert.IsFalse(samples["127N"].IsReferenceChannel);

            // Fraction and technical replicate are per-file, not per-channel.
            Assert.IsTrue(samples.Values.All(s => s.Fraction == 2 && s.TechnicalReplicate == 3));
        }

        /// <summary>
        /// A channel the design does not mention still gets an entry, because MetaMorpheus reports an
        /// intensity for every channel of the tag. Skipping it would shift every later channel by one.
        /// </summary>
        [Test]
        public static void ToMzLibDesign_UnannotatedChannelsBecomeEmptyRatherThanBeingSkipped()
        {
            var tag = Tmt10();
            var file = FileWith(@"C:\data\run1.raw", "PlexA", 1, 1,
                ("126", "Reference", 1, TmtSampleType.Reference));

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { file }, tag, out var errors);
            Assert.IsEmpty(errors);

            var samples = design.FileNameSampleInfoDictionary["run1.raw"];
            Assert.AreEqual(10, samples.Length);

            var unannotated = samples.Cast<IsobaricQuantSampleInfo>().Where(s => s.ChannelLabel != "126").ToList();
            Assert.AreEqual(9, unannotated.Count);
            Assert.IsTrue(unannotated.All(s => s.Condition == string.Empty));
            Assert.IsTrue(unannotated.All(s => !s.IsReferenceChannel));
        }

        [Test]
        public static void ToMzLibDesign_ChannelThatIsNotPartOfTheTag_IsAnError()
        {
            var tag = Tmt10();
            var file = FileWith(@"C:\data\run1.raw", "PlexA", 1, 1,
                ("135N", "Treated", 1, TmtSampleType.StudySample));   // TMT18 channel, not TMT10

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { file }, tag, out var errors);

            Assert.IsNull(design, "a design with an impossible channel must not be half-usable");
            Assert.IsTrue(errors.Any(e => e.Contains("135N")), string.Join(" | ", errors));
        }

        /// <summary>Plex is free text in the file and an int on the sample info; the mapping must be stable.</summary>
        [Test]
        public static void ToMzLibDesign_AssignsStablePlexIds()
        {
            var tag = Tmt10();
            var files = new[]
            {
                FileWith(@"C:\data\runB.raw", "PlexB", 1, 1, ("126", "A", 1, TmtSampleType.StudySample)),
                FileWith(@"C:\data\runA.raw", "PlexA", 1, 1, ("126", "A", 1, TmtSampleType.StudySample)),
                FileWith(@"C:\data\runA2.raw", "PlexA", 1, 2, ("126", "A", 1, TmtSampleType.StudySample))
            };

            var design = TmtExperimentalDesign.ToMzLibDesign(files, tag, out var errors);
            Assert.IsEmpty(errors);

            int PlexOf(string fileName) =>
                ((IsobaricQuantSampleInfo)design.FileNameSampleInfoDictionary[fileName][0]).PlexId;

            Assert.AreEqual(PlexOf("runA.raw"), PlexOf("runA2.raw"), "same plex, same id");
            Assert.AreNotEqual(PlexOf("runA.raw"), PlexOf("runB.raw"));
            Assert.AreEqual(1, PlexOf("runA.raw"), "ids follow the sorted plex names");
            Assert.AreEqual(2, PlexOf("runB.raw"));
        }

        /// <summary>The dictionary key must be the file name WITH extension, which is what PivotByFile looks up.</summary>
        [Test]
        public static void ToMzLibDesign_KeysByFileNameWithExtension()
        {
            var tag = Tmt10();
            var file = FileWith(@"C:\some\deep\path\run1.raw", "PlexA", 1, 1,
                ("126", "A", 1, TmtSampleType.StudySample));

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { file }, tag, out _);

            Assert.IsTrue(design.FileNameSampleInfoDictionary.ContainsKey("run1.raw"));
            Assert.AreEqual(1, design.FileNameSampleInfoDictionary.Count);
        }

        [Test]
        public static void ToMzLibDesign_WithoutATag_IsAnError()
        {
            var file = FileWith(@"C:\data\run1.raw", "PlexA", 1, 1,
                ("126", "A", 1, TmtSampleType.StudySample));

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { file }, null, out var errors);

            Assert.IsNull(design);
            Assert.IsTrue(errors.Any(e => e.Contains("isobaric mass tag")), string.Join(" | ", errors));
        }

        [Test]
        public static void ToMzLibDesign_WithNoFiles_IsAnError()
        {
            var design = TmtExperimentalDesign.ToMzLibDesign(new List<TmtFileInfo>(), Tmt10(), out var errors);

            Assert.IsNull(design);
            Assert.IsNotEmpty(errors);
        }

        /// <summary>
        /// The window offering the channel drop-down and the tag that validates it must agree, or the
        /// GUI can write a design that cannot be projected. They were two separate tables until now.
        /// </summary>
        [Test]
        public static void EveryTagsLabelsMatchItsReporterIonCount()
        {
            foreach (IsobaricMassTagType type in Enum.GetValues(typeof(IsobaricMassTagType)))
            {
                var tag = IsobaricMassTag.GetIsobaricMassTag(type);
                if (tag == null) continue;   // modification not loaded in this environment

                var labels = IsobaricMassTag.GetReporterIonLabels(type);
                Assert.IsNotNull(labels, type.ToString());
                Assert.AreEqual(labels.Count, tag.ReporterIonMzs.Length,
                    $"{type} has {labels.Count} labels but {tag.ReporterIonMzs.Length} reporter ions");
            }
        }

        /// <summary>
        /// The stronger half of the same invariant, and the one that actually protects the projection:
        /// ToMzLibDesign pairs channelLabels[i] with ReporterIonMzs[i], so label i has to name the
        /// channel whose m/z sits at i. A count cannot catch a mislabelled channel -- iTRAQ 8-plex
        /// carried the name "120" for the 121 reagent for exactly that reason, extracting the right ion
        /// under a name no kit has.
        ///
        /// Every label begins with its nominal mass, so the assertion is available cheaply.
        /// </summary>
        [Test]
        public static void EveryTagsLabelsNameTheChannelAtTheirOwnIndex()
        {
            foreach (IsobaricMassTagType type in Enum.GetValues(typeof(IsobaricMassTagType)))
            {
                var tag = IsobaricMassTag.GetIsobaricMassTag(type);
                if (tag == null) continue;   // modification not loaded in this environment

                var labels = IsobaricMassTag.GetReporterIonLabels(type);
                Assert.IsNotNull(labels, type.ToString());

                for (int i = 0; i < labels.Count; i++)
                {
                    string digits = new string(labels[i].TakeWhile(char.IsDigit).ToArray());
                    Assert.IsNotEmpty(digits, $"{type} label '{labels[i]}' does not begin with a nominal mass");

                    int nominal = int.Parse(digits);
                    int observed = (int)Math.Round(tag.ReporterIonMzs[i]);

                    Assert.AreEqual(nominal, observed,
                        $"{type} label '{labels[i]}' at index {i} names channel {nominal}, " +
                        $"but the reporter ion at that index is {tag.ReporterIonMzs[i]:F4} (channel {observed})");
                }
            }
        }

        /// <summary>
        /// ToMzLibDesign pairs channelLabels[i] with ReporterIonMzs[i] positionally. If the two ever
        /// differ in length that pairing is meaningless, and the failure is silent: every channel past
        /// the shorter array would be assigned the wrong reporter ion, producing a design that looks
        /// complete and reports the wrong sample's abundance. The guard has to refuse rather than
        /// truncate, and this is what proves it does.
        ///
        /// The two lengths come from independent sources -- the DI block in Mods/tmt.txt and the label
        /// table in IsobaricMassTag -- which is exactly how the iTRAQ 8-plex name drifted. A real tag
        /// cannot currently reach this state (EveryTagsLabelsMatchItsReporterIonCount asserts as much),
        /// so the mismatch is constructed by reflection.
        /// </summary>
        [Test]
        public static void MismatchedChannelAndReporterIonCountsAreRefused()
        {
            var tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT11);
            Assert.That(tag, Is.Not.Null, "TMT11 must load, or this test proves nothing");

            int realCount = tag.ReporterIonMzs.Length;
            typeof(IsobaricMassTag)
                .GetProperty(nameof(IsobaricMassTag.ReporterIonMzs))!
                .SetValue(tag, tag.ReporterIonMzs.Take(realCount - 1).ToArray());

            var design = TmtExperimentalDesign.ToMzLibDesign(new[] { OneFile("a.raw", "Plex1", "126") }, tag, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(design, Is.Null, "a design must not be produced from a misaligned tag");
                Assert.That(errors.Any(e => e.Contains("disagree")), Is.True);
                Assert.That(errors.Any(e => e.Contains($"{realCount} channels")), Is.True,
                    "the error should name both counts, since the point is which source to fix");
            });
        }

        /// <summary>
        /// The design is keyed by file name, so two entries resolving to the same name would have one
        /// silently overwrite the other's channels.
        ///
        /// The projection is all-or-nothing: any error yields a null design rather than a partial one.
        /// That is the right call here -- a design missing a file quantifies that file's channels
        /// against nothing and says so nowhere -- and the per-entry `continue` exists to collect every
        /// problem in one pass rather than to salvage the rest.
        /// </summary>
        [Test]
        public static void TwoEntriesSharingAFileNameAreRefusedWithoutLosingTheRest()
        {
            var tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT11);
            Assert.That(tag, Is.Not.Null);

            var files = new[]
            {
                OneFile(Path.Combine("runA", "sample.raw"), "Plex1", "126"),
                OneFile(Path.Combine("runB", "sample.raw"), "Plex1", "127N"),   // same name, different folder
                OneFile("other.raw", "Plex1", "127C"),
            };

            var design = TmtExperimentalDesign.ToMzLibDesign(files, tag, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors.Any(e => e.Contains("share the file name")), Is.True);
                Assert.That(design, Is.Null,
                    "a colliding key must fail the whole projection, not yield a design missing a file");
            });
        }

        /// <summary>
        /// An entry whose path has no file-name component cannot be keyed at all. Reported, and the
        /// whole projection refused, rather than producing an entry under an empty key that would
        /// never match a spectra file.
        /// </summary>
        [Test]
        public static void EntryWithNoFileNameIsReported()
        {
            var tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT11);
            Assert.That(tag, Is.Not.Null);

            var files = new[]
            {
                OneFile(Path.Combine("some", "folder") + Path.DirectorySeparatorChar, "Plex1", "126"),
                OneFile("good.raw", "Plex1", "127N"),
            };

            var design = TmtExperimentalDesign.ToMzLibDesign(files, tag, out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors.Any(e => e.Contains("no file name")), Is.True);
                Assert.That(design, Is.Null);
            });
        }

        /// <summary>
        /// Write round-trips a sample type through ToDesignFileValue, so an unrecognised value has to
        /// produce something Read will accept rather than an empty cell that fails to parse on reload.
        /// </summary>
        [Test]
        public static void UnknownSampleTypeWritesAsStudySample()
        {
            string written = TmtExperimentalDesign.ToDesignFileValue((TmtSampleType)999);

            Assert.That(written, Is.EqualTo("study sample"));
            Assert.That(TmtExperimentalDesign.TryParseSampleType(written, out var parsed), Is.True,
                "whatever is written must survive a read back");
            Assert.That(parsed, Is.EqualTo(TmtSampleType.StudySample));
        }

        /// <summary>A single-file, single-channel design entry, for the checks above.</summary>
        private static TmtFileInfo OneFile(string path, string plex, string tag) =>
            new TmtFileInfo(path, plex, 1, 1, new[]
            {
                new TmtPlexAnnotation
                {
                    Tag = tag,
                    SampleName = "S1",
                    Condition = "Control",
                    BiologicalReplicate = 1,
                    SampleType = TmtSampleType.StudySample
                }
            });

        #endregion
    }
}
