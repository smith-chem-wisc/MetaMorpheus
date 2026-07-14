using GuiFunctions;
using GuiFunctions.MetaDraw;
using NUnit.Framework;
using System;
using System.Windows.Media;

namespace Test.MetaDraw;

[TestFixture]
public class BioPolymerCoverageColorMapperFactoryTests
{
    private static BioPolymerCoverageResultModel MakeResult(double? precursorIntensity, double score, string fileName = "file1.psmtsv")
    {
        var psm = new DummySpectralmatch(filePath: fileName);
        psm.SetScore(score);
        if (precursorIntensity.HasValue)
        {
            var backing = typeof(Readers.SpectrumMatchFromTsv)
                .GetField("<PrecursorIntensity>k__BackingField",
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            backing.SetValue(psm, precursorIntensity.Value);
        }
        return new BioPolymerCoverageResultModel(psm, "ABC", 1, 2, BioPolymerCoverageType.Unique);
    }

    [Test]
    public void Create_None_ReturnsCategoricalMapper()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.None);
        Assert.That(mapper, Is.InstanceOf<CategoricalBioPolymerCoverageColorMapper>());
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.None));
        Assert.That(mapper.IsNumeric, Is.False);
    }

    [Test]
    public void Create_None_GetBrush_ReturnsGray()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.None);
        var result = MakeResult(null, 10);
        Assert.That(mapper.GetBrush(result, null).Color, Is.EqualTo(Colors.Gray));
    }

    [Test]
    public void Create_CoverageType_ReturnsCategoricalMapper()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.CoverageType);
        Assert.That(mapper, Is.InstanceOf<CategoricalBioPolymerCoverageColorMapper>());
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.CoverageType));
        Assert.That(mapper.IsNumeric, Is.False);
    }

    [Test]
    public void Create_CoverageType_GetBrush_UsesCoverageColor()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.CoverageType);
        var result = MakeResult(null, 10);
        var expected = MetaDrawSettings.BioPolymerCoverageColors[BioPolymerCoverageType.Unique];
        Assert.That(mapper.GetBrush(result, null), Is.EqualTo(expected));
    }

    [Test]
    public void Create_FileOrigin_ReturnsCategoricalMapper()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.FileOrigin);
        Assert.That(mapper, Is.InstanceOf<FileOriginColorMapper>());
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.FileOrigin));
        Assert.That(mapper.IsNumeric, Is.False);
    }

    [Test]
    public void Create_FileOrigin_GetBrush_UsesStablePerFileColor()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.FileOrigin);
        var result = MakeResult(null, 10, fileName: "file1.psmtsv");
        var first = mapper.GetBrush(result, null);
        var second = mapper.GetBrush(result, null);
        Assert.That(second, Is.EqualTo(first));
    }

    [Test]
    public void Create_FileOrigin_EmptyFileName_ReturnsGray()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.FileOrigin);
        var result = MakeResult(null, 10, fileName: "");
        Assert.That(mapper.GetBrush(result, null).Color, Is.EqualTo(Colors.Gray));
    }

    [Test]
    public void Create_FileOrigin_DifferentFiles_GetDifferentColors()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.FileOrigin);
        var result1 = MakeResult(null, 10, fileName: "file1.psmtsv");
        var result2 = MakeResult(null, 10, fileName: "file2.psmtsv");
        Assert.That(mapper.GetBrush(result1, null), Is.Not.EqualTo(mapper.GetBrush(result2, null)));
    }

    [Test]
    public void Create_PrecursorIntensity_ReturnsNumericMapper()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.PrecursorIntensity);
        Assert.That(mapper, Is.InstanceOf<PrecursorIntensityColorMapper>());
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.PrecursorIntensity));
        Assert.That(mapper.IsNumeric, Is.True);
    }

    [Test]
    public void Create_Score_ReturnsNumericMapper()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.Score);
        Assert.That(mapper, Is.InstanceOf<ScoreColorMapper>());
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.Score));
        Assert.That(mapper.IsNumeric, Is.True);
    }

    [Test]
    public void Create_UnknownColorBy_Throws()
    {
        var unknown = (ColorResultsBy)9999;
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BioPolymerCoverageColorMapperFactory.Create(unknown));
    }
}
