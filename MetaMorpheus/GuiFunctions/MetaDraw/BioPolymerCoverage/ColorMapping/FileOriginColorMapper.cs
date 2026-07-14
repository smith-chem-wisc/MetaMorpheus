using GuiFunctions.Util;
using System.Collections.Generic;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class FileOriginColorMapper : BioPolymerCoverageColorMapper
{
    private readonly CyclicalQueue<SolidColorBrush> _colorQueue = new(new[]
    {
        new SolidColorBrush(Colors.Red),
        new SolidColorBrush(Colors.Blue),
        new SolidColorBrush(Colors.Green),
        new SolidColorBrush(Colors.Orange),
        new SolidColorBrush(Colors.Purple),
        new SolidColorBrush(Colors.Teal),
        new SolidColorBrush(Colors.Brown),
        new SolidColorBrush(Colors.Pink),
        new SolidColorBrush(Colors.Yellow),
        new SolidColorBrush(Colors.Gray),
        new SolidColorBrush(Colors.Cyan),
        new SolidColorBrush(Colors.Magenta),
        new SolidColorBrush(Colors.LimeGreen),
        new SolidColorBrush(Colors.DarkBlue),
        new SolidColorBrush(Colors.DarkRed),
        new SolidColorBrush(Colors.DarkGreen),
        new SolidColorBrush(Colors.Gold),
        new SolidColorBrush(Colors.Indigo),
        new SolidColorBrush(Colors.Olive),
        new SolidColorBrush(Colors.Maroon),
        new SolidColorBrush(Colors.Navy),
        new SolidColorBrush(Colors.Turquoise),
        new SolidColorBrush(Colors.Violet),
        new SolidColorBrush(Colors.Sienna),
        new SolidColorBrush(Colors.Salmon),
        new SolidColorBrush(Colors.Coral),
        new SolidColorBrush(Colors.Khaki),
        new SolidColorBrush(Colors.Plum),
        new SolidColorBrush(Colors.Peru),
        new SolidColorBrush(Colors.SteelBlue),
        new SolidColorBrush(Colors.MediumPurple),
        new SolidColorBrush(Colors.MediumSeaGreen),
        new SolidColorBrush(Colors.MediumSlateBlue),
        new SolidColorBrush(Colors.MediumVioletRed),
        new SolidColorBrush(Colors.MediumOrchid),
        new SolidColorBrush(Colors.MediumTurquoise),
        new SolidColorBrush(Colors.MediumSpringGreen),
        new SolidColorBrush(Colors.MediumAquamarine)
    });

    private readonly Dictionary<string, SolidColorBrush> _identifierToColor = [];

    public override ColorResultsBy ColorBy => ColorResultsBy.FileOrigin;
    public override bool IsNumeric => false;
    public override bool SupportsGradientSelection => false;
    public override bool SupportsLogScale => false;

    public override SolidColorBrush GetBrush(BioPolymerCoverageResultModel result, BioPolymerCoverageColorScale? scale)
    {
        var identifier = result?.Match?.FileName;
        if (string.IsNullOrEmpty(identifier))
            return new SolidColorBrush(Colors.Gray);

        if (!_identifierToColor.TryGetValue(identifier, out var brush))
        {
            brush = _colorQueue.Dequeue();
            _identifierToColor[identifier] = brush;
        }

        return brush;
    }
}
