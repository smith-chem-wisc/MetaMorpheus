using NUnit.Framework;
using TaskLayer.ParallelSearch.Util.Converters;

namespace Test.ParallelSearchTask.Util;

[TestFixture]
public class ConverterTypeConverterTests
{
    [Test]
    public void SemiColonDelimitedToDoubleArray_ParsesAndRounds()
    {
        var converter = new SemiColonDelimitedToDoubleArrayTypeConverter();
        SemiColonDelimitedToDoubleArrayTypeConverter.RoundingPlaces = 3;

        var parsed = (double[])converter.ConvertFromString("1.25;;2.5;3", null!, null!);
        string written = converter.ConvertToString(new[] { 1.23456, 9.87654 }, null!, null!);

        Assert.Multiple(() =>
        {
            Assert.That(parsed, Is.EquivalentTo(new[] { 1.25, 2.5, 3.0 }));
            Assert.That(written, Is.EqualTo("1.235;9.877"));
        });
    }

    [Test]
    public void SemiColonDelimitedToDoubleArray_ConvertToStringWithWrongType_WritesEmptyField()
    {
        var converter = new SemiColonDelimitedToDoubleArrayTypeConverter();
        Assert.That(converter.ConvertToString("not-a-list", null!, null!), Is.Empty);
    }
}
