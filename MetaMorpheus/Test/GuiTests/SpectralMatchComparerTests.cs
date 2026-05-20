using System.IO;
using GuiFunctions;
using GuiFunctions.Util;
using NUnit.Framework;
using Readers;

namespace Test.GuiTests;

[TestFixture]
public class SpectralMatchComparerTests
{
    private PsmFromTsv _psm;
    private OsmFromTsv _osm;

    [OneTimeSetUp]
    public void OneTimeSetUp()
    {
        var psmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SequenceCoverageTestPSM.psmtsv");
        var osmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "AllOligos.osmtsv");

        _psm = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(psmPath, out _)[0];
        _osm = SpectrumMatchTsvReader.ReadTsv<OsmFromTsv>(osmPath, out _)[0];
    }

    [Test]
    public void Equals_ReturnsTrue_ForSameReference()
    {
        Assert.That(SpectralMatchComparer.Instance.Equals(_psm, _psm), Is.True);
    }

    [Test]
    public void Equals_ReturnsFalse_WhenLeftIsNull()
    {
        Assert.That(SpectralMatchComparer.Instance.Equals(null, _psm), Is.False);
    }

    [Test]
    public void Equals_ReturnsFalse_WhenRightIsNull()
    {
        Assert.That(SpectralMatchComparer.Instance.Equals(_psm, null), Is.False);
    }

    [Test]
    public void Equals_ReturnsFalse_WhenTypesDiffer()
    {
        Assert.That(SpectralMatchComparer.Instance.Equals(_psm, _osm), Is.False);
    }

    [Test]
    public void Equals_ReturnsTrue_WhenAllComparedFieldsMatch()
    {
        var psmCopy = _psm.ReplaceFullSequence(_psm.FullSequence, _psm.BaseSeq);

        Assert.That(SpectralMatchComparer.Instance.Equals(_psm, psmCopy), Is.True);
        Assert.That(SpectralMatchComparer.Instance.GetHashCode(_psm), Is.EqualTo(SpectralMatchComparer.Instance.GetHashCode(psmCopy)));
    }

    [Test]
    public void Equals_ReturnsFalse_WhenAnyComparedFieldDiffers()
    {
        var psmDifferent = _psm.ReplaceFullSequence(_psm.FullSequence + "_ALT", _psm.BaseSeq);

        Assert.That(SpectralMatchComparer.Instance.Equals(_psm, psmDifferent), Is.False);
    }
}
