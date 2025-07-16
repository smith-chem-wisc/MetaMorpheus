using System.Linq;
using System.Threading;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;
using GuiFunctions;
using GuiFunctions.MetaDraw.Chimeras;
using NUnit.Framework;

namespace Test.MetaDraw;

[TestFixture, Apartment(ApartmentState.STA)]
public class ChimeraDrawnSequenceTests
{
    [Test]
    public void ChimeraDrawnSequence_CreatesCanvasWithCorrectDimensions()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var canvas = new Canvas();

        // Act
        var drawnSequence = new ChimeraDrawnSequence(canvas, chimeraGroup);

        // Assert
        int longestSequenceLength = chimeraGroup.ChimericPsms.Max(psm => psm.Psm.BaseSeq.Split('|')[0].Length);
        double expectedWidth = (longestSequenceLength + 4) * MetaDrawSettings.AnnotatedSequenceTextSpacing + 20;
        double expectedHeight = 40 * chimeraGroup.ChimericPsms.Count + 20;
        Assert.That(canvas.Width, Is.EqualTo(expectedWidth));
        Assert.That(canvas.Height, Is.EqualTo(expectedHeight));
    }

    [Test]
    public void ChimeraDrawnSequence_DrawsAllSequences()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        var canvas = new Canvas();

        // Act
        var drawnSequence = new ChimeraDrawnSequence(canvas, chimeraGroup);

        // Assert
        // At least one text element per residue per sequence
        int minTextBlocks = chimeraGroup.ChimericPsms.Sum(psm => psm.Psm.BaseSeq.Split('|')[0].Length);
        int actualTextBlocks = canvas.Children.OfType<System.Windows.Controls.TextBlock>().Count();
        Assert.That(actualTextBlocks, Is.GreaterThanOrEqualTo(minTextBlocks));
    }

    [Test]
    public void ChimeraDrawnSequence_AddsMatchedIonsWithCorrectColors()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var expectedIons = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ExpectedIonsByColor;
        var canvas = new Canvas();

        // Act
        var drawnSequence = new ChimeraDrawnSequence(canvas, chimeraGroup);

        // Assert
        // For each expected color, check that at least one shape with that color exists
        foreach (var color in expectedIons.Keys)
        {
            var wpfColor = DrawnSequence.ParseColorFromOxyColor(color);
            bool found = canvas.Children.OfType<Polyline>().Any(line => (line.Stroke as SolidColorBrush)?.Color == wpfColor);
            Assert.That(found, Is.True, $"Expected to find a shape with color {wpfColor} for ion color {color}.");
        }
    }
}