using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Proteomics;

namespace Test
{
    /// <summary>
    /// ProteinGroup exposes two SpectraFileInfo-shaped views over the base class's sample-keyed
    /// quantification store: FilesForQuantification and IntensitiesByFile. Isobaric quantification
    /// puts IsobaricQuantSampleInfo into that store -- one entry per channel rather than per file --
    /// so both views have to answer for a sample they cannot represent, and they have to answer the
    /// same way.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class ProteinGroupIsobaricViewTests
    {
        private const string File = "tmt.raw";

        private static ProteinGroup GroupWith(params ISampleInfo[] samples)
        {
            var protein = new Protein("PEPTIDEK", "P1");
            var group = new ProteinGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            group.SamplesForQuantification = samples.ToList();
            group.IntensitiesBySample = samples
                .Select((s, i) => (s, v: 100.0 * (i + 1)))
                .ToDictionary(t => t.s, t => t.v);

            return group;
        }

        private static IsobaricQuantSampleInfo Channel(string label, double mz) =>
            new IsobaricQuantSampleInfo(File, "Control", 0, 0, 0, 0, label, mz, false);

        [Test]
        public void IntensitiesByFile_DropsIsobaricChannelsInsteadOfThrowing()
        {
            var group = GroupWith(Channel("126", 126.0), Channel("127N", 127.1));

            // Previously a hard cast: reading this after an isobaric search threw
            // InvalidCastException from inside the ordinary output path.
            IReadOnlyDictionary<SpectraFileInfo, double> byFile = null;
            Assert.DoesNotThrow(() => byFile = group.IntensitiesByFile);
            Assert.That(byFile, Is.Empty);
        }

        [Test]
        public void TheTwoFileViewsAgreeOnWhatTheyDrop()
        {
            var labelFree = new SpectraFileInfo(File, "Control", 0, 0, 0);
            var group = GroupWith(labelFree, Channel("126", 126.0), Channel("127N", 127.1));

            // FilesForQuantification has always filtered; IntensitiesByFile now does too, so a caller
            // reading both sees the same set of files rather than one view succeeding and the other
            // throwing on identical data.
            Assert.Multiple(() =>
            {
                Assert.That(group.FilesForQuantification, Is.EqualTo(new[] { labelFree }));
                Assert.That(group.IntensitiesByFile.Keys, Is.EqualTo(new[] { labelFree }));
                Assert.That(group.IntensitiesByFile[labelFree], Is.EqualTo(100.0));
            });
        }

        [Test]
        public void LabelFreeQuantificationIsUnaffected()
        {
            var a = new SpectraFileInfo("a.raw", "Control", 0, 0, 0);
            var b = new SpectraFileInfo("b.raw", "Treated", 1, 0, 0);
            var group = GroupWith(a, b);

            Assert.Multiple(() =>
            {
                Assert.That(group.IntensitiesByFile, Has.Count.EqualTo(2));
                Assert.That(group.IntensitiesByFile[a], Is.EqualTo(100.0));
                Assert.That(group.IntensitiesByFile[b], Is.EqualTo(200.0));
            });
        }
    }
}
