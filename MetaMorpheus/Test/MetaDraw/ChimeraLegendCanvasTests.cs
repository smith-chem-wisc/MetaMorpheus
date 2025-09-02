using System;
using System.Linq;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using NUnit.Framework;

namespace Test.MetaDraw;

[TestFixture, Apartment(System.Threading.ApartmentState.STA)]
public class ChimeraLegendCanvasTests
{
    [Test]
    public void Constructor_BuildsLegendWithExpectedChildren()
    {
        // Arrange
        var group = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;

        // Act
        var legend = new ChimeraLegendCanvas(group);

        // Assert
        // Should have at least one ellipse and one text block per proteoform
        int proteoformCount = group.ChimericPsms.Count;
        int ellipseCount = legend.Children.OfType<Ellipse>().Count();
        int textCount = legend.Children.OfType<TextBlock>().Count();
        Assert.That(ellipseCount, Is.GreaterThanOrEqualTo(1));
        Assert.That(textCount, Is.GreaterThanOrEqualTo(1));
        Assert.That(legend.Width, Is.GreaterThan(0));
        Assert.That(legend.Height, Is.GreaterThan(0));
    }

    [Test]
    [NonParallelizable]
    public void BuildLegend_WrapsLongText()
    {
        // Arrange
        var originalParam = MetaDrawSettings.ChimeraLegendMaxWidth;
        MetaDrawSettings.ChimeraLegendMaxWidth = 20;
        var group = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        // Simulate a long protein name
        foreach (var psm in group.ChimericPsms)
        {
            var nameProp = psm.Psm.GetType().GetProperty("Name", System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance);
            nameProp?.SetValue(psm.Psm, new string('A', 5000));
        }

        // Act
        var legend = new ChimeraLegendCanvas(group);

        // Assert
        var textBlock = legend.Children.OfType<TextBlock>().FirstOrDefault();
        Assert.That(textBlock, Is.Not.Null);
        Assert.That(textBlock.TextWrapping, Is.EqualTo(System.Windows.TextWrapping.Wrap));
        Assert.That(textBlock.MaxWidth, Is.GreaterThan(0));
        // Should wrap: height should be more than one line
        textBlock.Measure(new System.Windows.Size(textBlock.MaxWidth, double.PositiveInfinity));
        Assert.That(textBlock.DesiredSize.Height, Is.GreaterThan(16));
        MetaDrawSettings.ChimeraLegendMaxWidth = originalParam;
    }

    [Test]
    public void BuildLegend_HandlesAmbiguousNames()
    {
        // Arrange
        var group = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        foreach (var psm in group.ChimericPsms)
        {
            var nameProp = psm.Psm.GetType().GetProperty("Name", System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance);
            nameProp?.SetValue(psm.Psm, "A|B|C");
        }

        MetaDrawSettings.ChimeraLegendTakeFirstIfAmbiguous = true;

        // Act
        var legend = new ChimeraLegendCanvas(group);

        // Assert
        var textBlock = legend.Children.OfType<TextBlock>().FirstOrDefault();
        Assert.That(textBlock, Is.Not.Null);
        Assert.That(textBlock.Text, Is.EqualTo("A"));
    }

    [Test]
    public void BuildLegend_HandlesNullOrEmptyModString()
    {
        // Arrange
        var group = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        foreach (var psm in group.ChimericPsms)
            psm.ModString = null;

        MetaDrawSettings.ChimeraLegendMainTextType = LegendDisplayProperty.Modifications;

        // Act
        var legend = new ChimeraLegendCanvas(group);

        // Assert
        var textBlock = legend.Children.OfType<TextBlock>().FirstOrDefault();
        Assert.That(textBlock, Is.Not.Null);
        Assert.That(textBlock.Text, Is.EqualTo("Unmodified"));
    }

    [Test]
    public void MeasureTextWidth_ReturnsPositiveWidth()
    {
        // Arrange
        var group = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var legend = new ChimeraLegendCanvas(group);

        // Act
        double width = legend.GetType()
            .GetMethod("MeasureTextWidth", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            .Invoke(legend, new object[] { "TestText", 12.0, System.Windows.FontWeights.Regular }) as double? ?? 0;

        // Assert
        Assert.That(width, Is.GreaterThan(0));
    }

    [Test]
    public void SanitizeIfAmbiguous_ReturnsFirstPartIfAmbiguous()
    {
        // Arrange
        var group = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var legend = new ChimeraLegendCanvas(group);

        // Act
        var method = legend.GetType().GetMethod("SanitizeIfAmbiguous", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        string result = method.Invoke(legend, new object[] { "A|B|C" }) as string;

        // Assert
        Assert.That(result, Is.EqualTo("A"));
    }

    [Test]
    public void BuildLegend_HandlesAllLegendDisplayPropertyCases_ForMainAndSubText_UsingCachedPsmValues()
    {
        var group = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var psms = group.ChimericPsms.Select(p => p.Psm).ToList();
        var modStrings = group.ChimericPsms.Select(p => p.ModString ?? "Unmodified").ToList();

        // Cache all property values for each PSM
        var cached = group.ChimericPsms.Select(p => new
        {
            Name = p.Psm.Name,
            Accession = p.Psm.Accession,
            BaseSeq = p.Psm.BaseSeq,
            FullSequence = p.Psm.FullSequence,
            ModString = p.ModString ?? "Unmodified"
        }).ToList();

        // Test all main text types
        foreach (LegendDisplayProperty mainType in Enum.GetValues(typeof(LegendDisplayProperty)))
        {
            MetaDrawSettings.ChimeraLegendMainTextType = mainType;
            var legend = new ChimeraLegendCanvas(group);

            // The main text is the header for the group (first proteoform)
            var mainTextBlock = legend.Children.OfType<TextBlock>().FirstOrDefault();
            Assert.That(mainTextBlock, Is.Not.Null);

            string expectedMain = mainType switch
            {
                LegendDisplayProperty.ProteinName => cached[0].Name,
                LegendDisplayProperty.ProteinAccession => cached[0].Accession,
                LegendDisplayProperty.BaseSequence => cached[0].BaseSeq,
                LegendDisplayProperty.FullSequence => cached[0].FullSequence,
                LegendDisplayProperty.Modifications => cached[0].ModString,
                _ => throw new NotImplementedException()
            };
            Assert.That(mainTextBlock.Text, Is.EqualTo(expectedMain));
        }

        // Test all sub text types
        foreach (LegendDisplayProperty subType in Enum.GetValues(typeof(LegendDisplayProperty)))
        {
            MetaDrawSettings.ChimeraLegendSubTextType = subType;
            // Use a main type that triggers subtexts (e.g., ProteinName)
            MetaDrawSettings.ChimeraLegendMainTextType = LegendDisplayProperty.ProteinName;
            var legend = new ChimeraLegendCanvas(group);

            // The subtexts are after the header, so skip the first TextBlock (header)
            var subTextBlocks = legend.Children.OfType<TextBlock>().Skip(1).ToList();
            Assert.That(subTextBlocks.Count, Is.EqualTo(group.ChimericPsms.Count));

            string[] expectedSubs = cached.Select(c =>
                subType switch
                {
                    LegendDisplayProperty.ProteinName => c.Name,
                    LegendDisplayProperty.ProteinAccession => c.Accession,
                    LegendDisplayProperty.BaseSequence => c.BaseSeq,
                    LegendDisplayProperty.FullSequence => c.FullSequence,
                    LegendDisplayProperty.Modifications => c.ModString,
                    _ => throw new NotImplementedException()
                }
            ).ToArray();

            for (int i = 0; i < subTextBlocks.Count; i++)
            {
                Assert.That(subTextBlocks[i].Text, Is.EqualTo(expectedSubs[i]));
            }
        }
    }
}
