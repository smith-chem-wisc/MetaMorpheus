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

    [TestCase("FullSequence")]
    [TestCase("Ms2ScanNumber")]
    [TestCase("FileNameWithoutExtension")]
    [TestCase("PrecursorScanNum")]
    public void Equals_ReturnsFalse_WhenFieldDiffers(string fieldName)
    {
        var basePsm = new MutablePsmFromTsv(_psm);
        var variantPsm = new MutablePsmFromTsv(_psm);

        switch (fieldName)
        {
            case "FullSequence":
                variantPsm = new MutablePsmFromTsv(_psm, _psm.FullSequence + "_ALT", _psm.BaseSeq);
                break;
            case "Ms2ScanNumber":
                variantPsm.SetMs2ScanNumber(_psm.Ms2ScanNumber + 1);
                break;
            case "FileNameWithoutExtension":
                variantPsm.SetFileNameWithoutExtension(_psm.FileNameWithoutExtension + "_ALT");
                break;
            case "PrecursorScanNum":
                variantPsm.SetPrecursorScanNum(_psm.PrecursorScanNum + 1);
                break;
        }

        Assert.That(SpectralMatchComparer.Instance.Equals(basePsm, variantPsm), Is.False);
    }

    /// <summary>
    /// Test helper that exposes protected setters on PsmFromTsv for field-level isolation.
    /// </summary>
    internal class MutablePsmFromTsv : PsmFromTsv
    {
        public MutablePsmFromTsv(PsmFromTsv source) : base(source, source.FullSequence) { }

        public MutablePsmFromTsv(PsmFromTsv source, string fullSequence, string baseSequence = "")
            : base(source, fullSequence, baseSequence: baseSequence) { }

        public void SetMs2ScanNumber(int value) { Ms2ScanNumber = value; }
        public void SetFileNameWithoutExtension(string value) { FileNameWithoutExtension = value; }
        public void SetPrecursorScanNum(int value) { PrecursorScanNum = value; }
    }
}
