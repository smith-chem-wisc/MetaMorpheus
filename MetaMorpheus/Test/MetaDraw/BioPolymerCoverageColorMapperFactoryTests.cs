using GuiFunctions;
using GuiFunctions.MetaDraw;
using GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;
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
        mapper.Prepare([result], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result).Color, Is.EqualTo(Colors.Gray));
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
        // Settings store OxyColor now, so the mapper converts on the way out and cannot hand back the
        // stored instance. Compare the colour rather than brush identity.
        var expected = MetaDrawSettings.BioPolymerCoverageColors[BioPolymerCoverageType.Unique];
        mapper.Prepare([result], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result).Color, Is.EqualTo(DrawnSequence.ParseColorFromOxyColor(expected)));
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
        mapper.Prepare([result], ColorGradientType.Viridis, false);
        var first = mapper.GetBrush(result);
        var second = mapper.GetBrush(result);
        Assert.That(second, Is.EqualTo(first));
    }

    [Test]
    public void Create_FileOrigin_EmptyFileName_ReturnsGray()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.FileOrigin);
        var result = MakeResult(null, 10, fileName: "");
        mapper.Prepare([result], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result).Color, Is.EqualTo(Colors.Gray));
    }

    [Test]
    public void Create_FileOrigin_DifferentFiles_GetDifferentColors()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.FileOrigin);
        var result1 = MakeResult(null, 10, fileName: "file1.psmtsv");
        var result2 = MakeResult(null, 10, fileName: "file2.psmtsv");
        mapper.Prepare([result1, result2], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result1), Is.Not.EqualTo(mapper.GetBrush(result2)));
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
