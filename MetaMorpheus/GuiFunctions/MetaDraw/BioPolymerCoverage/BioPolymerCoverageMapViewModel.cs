using GuiFunctions.Util;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public enum ColorResultsBy
{
    None,
    CoverageType,
    FileOrigin,
}

public class BioPolymerCoverageMapViewModel : BaseViewModel
{

    #region Color Handling

    private static readonly CyclicalQueue<SolidColorBrush> ColorQueue = new(new[]
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

    private static readonly Dictionary<string, SolidColorBrush> IdentifierToColor = [];

    public ColorResultsBy[] AllColorByTypes { get; } = Enum.GetValues<ColorResultsBy>();

    private ColorResultsBy _colorBy = ColorResultsBy.CoverageType;
    public ColorResultsBy ColorBy
    {
        get => _colorBy;
        set
        {
            if (_colorBy != value)
            {
                _colorBy = value;
                OnPropertyChanged(nameof(ColorBy));
                IdentifierToColor.Clear();
                ColorQueue.Reset();
                Redraw();
            }
        }
    }

    private SolidColorBrush GetColorBrush(BioPolymerCoverageResultModel result)
    {
        switch (_colorBy)
        {
            case ColorResultsBy.CoverageType:
                return MetaDrawSettings.BioPolymerCoverageColors[result.CoverageType];
            case ColorResultsBy.FileOrigin:
                // Use FileName as identifier since FullFilePath is not available
                return GetColorByIdentifier(result.Match.FileName);
            default:
                return new SolidColorBrush(Colors.Gray);
        }
    }

    // Helper method for FileOrigin coloring
    private SolidColorBrush GetColorByIdentifier(string identifier)
    {
        if (string.IsNullOrEmpty(identifier))
            return new SolidColorBrush(Colors.Gray);

        if (!IdentifierToColor.TryGetValue(identifier, out var brush))
        {
            brush = ColorQueue.Dequeue();
            IdentifierToColor[identifier] = brush;
        }
        return brush;
    }

    #endregion

    private DrawingImage _coverageDrawing;
    public DrawingImage CoverageDrawing
    {
        get => _coverageDrawing;
        internal set { _coverageDrawing = value; OnPropertyChanged(nameof(CoverageDrawing)); }
    }

    private BioPolymerGroupViewModel _group;
    public BioPolymerGroupViewModel Group
    {
        get => _group;
        set { _group = value; OnPropertyChanged(nameof(Group)); Redraw(); }
    }

    private int _lettersPerRow = 50;
    public int LettersPerRow
    {
        get => _lettersPerRow;
        set { _lettersPerRow = value; OnPropertyChanged(nameof(LettersPerRow)); }
    }

    private double _availableWidth = 800; // Default, should be set by view on resize
    public double AvailableWidth
    {
        get => _availableWidth;
        set { _availableWidth = value; OnPropertyChanged(nameof(AvailableWidth)); UpdateLettersPerRow(value); }
    }

    // Call this from the view when the canvas size changes
    public void UpdateLettersPerRow(double availableWidth)
    {
        _availableWidth = availableWidth;
        int newLetters = Math.Max(10, (int)(availableWidth / (MetaDrawSettings.BioPolymerCoverageFontSize * 0.70)));
        if (newLetters != LettersPerRow)
            LettersPerRow = newLetters;
        else
            Redraw();
    }

    private void Redraw()
    {
        if (Group == null)
        {
            CoverageDrawing = null;
            return;
        }

        var seq = Group.Sequence;
        var results = Group.CoverageResults;
        var filteredResults = results.Where(p => MetaDrawSettings.FilterAcceptsPsm(p.Match)).ToList();
        int lettersPerRow = LettersPerRow;
        double fontSize = MetaDrawSettings.BioPolymerCoverageFontSize;

        double plotMargin = fontSize * 2;
        double usableWidth = Math.Max(1, AvailableWidth - 2 * plotMargin);
        double letterWidth = usableWidth / lettersPerRow;
        double letterHeight = fontSize * 1.2;
        double rectHeight = fontSize * 0.8;
        double yPad = fontSize * 0.2;
        double rectJitter = rectHeight * 0.5;

        int nRows = (seq.Length + lettersPerRow - 1) / lettersPerRow;
        var rowRectangles = new Dictionary<int, List<(int startCol, int endCol)>>();
        var rowOccupancy = new Dictionary<int, int[]>(); // row -> per-column track usage

        var dv = new DrawingVisual();
        using (var dc = dv.RenderOpen())
        {
            double legendSpacing = fontSize * 0.5;
            double dpi = MetaDrawSettings.CanvasPdfExportDpi;
            double headerTop = legendSpacing;

            // --- Top line: global coverage metrics ---
            int seqLen = seq.Length;
            var anyCovered = new bool[seqLen];
            var uniqueCovered = new bool[seqLen];

            foreach (var r in filteredResults)
            {
                bool isUnique = r.CoverageType == BioPolymerCoverageType.Unique
                             || r.CoverageType == BioPolymerCoverageType.UniqueMissedCleavage;

                int s = Math.Max(1, r.Start) - 1;
                int e = Math.Min(seqLen, r.End) - 1;
                for (int i = s; i <= e; i++)
                {
                    anyCovered[i] = true;
                    if (isUnique) uniqueCovered[i] = true;
                }
            }
            int uniqueCoveredCount = uniqueCovered.Count(b => b);
            int maxCoveredCount = anyCovered.Count(b => b);

            var metricsText = $"Unique Coverage: {uniqueCoveredCount}/{seqLen} ({Math.Round(100 * uniqueCoveredCount / (double)seqLen, 2)}%)   |   Maximum Coverage: {maxCoveredCount}/{seqLen} ({Math.Round(100 * maxCoveredCount / (double)seqLen, 2)}%)";
            var metricsFt = new FormattedText(
                metricsText,
                System.Globalization.CultureInfo.CurrentCulture,
                FlowDirection.LeftToRight,
                new Typeface("Segoe UI"),
                fontSize * 0.9,
                Brushes.Black,
                dpi);

            double metricsX = plotMargin + (usableWidth - metricsFt.Width) / 2;
            double metricsY = headerTop;
            dc.DrawText(metricsFt, new Point(metricsX, metricsY));

            // --- Second line: legend, wrapped ---
            double legendY = metricsY + metricsFt.Height + legendSpacing;
            var legendItems = CreateLegendItems(filteredResults, fontSize, dpi, ColorBy);

            var legendLines = WrapLegendItems(legendItems, fontSize, dpi, usableWidth);

            double legendLineHeight = legendItems.Count == 0 ? 0 : Math.Max(fontSize * 0.8, legendItems.Max(li => li.Text.Height));
            double sepWidth = new FormattedText("\t", System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Segoe UI"), fontSize * 0.8, Brushes.DimGray, dpi).Width;

            foreach (var line in legendLines)
            {
                double lineWidth = line.Sum(item => fontSize * 0.8 + 6 + item.Text.Width) + sepWidth * (line.Count - 1);
                double legendX = plotMargin + (usableWidth - lineWidth) / 2;
                for (int i = 0; i < line.Count; i++)
                {
                    var (brush, text) = line[i];
                    var rect = new Rect(legendX, legendY + (legendLineHeight - fontSize * 0.8) / 2, fontSize * 0.8, fontSize * 0.8);
                    dc.DrawRoundedRectangle(brush, null, rect, 4, 4);
                    legendX += fontSize * 0.8 + 6;

                    double textY = legendY + (legendLineHeight - text.Height) / 2;
                    dc.DrawText(text, new Point(legendX, textY));
                    legendX += text.Width;

                    if (i < line.Count - 1)
                    {
                        double sepY = legendY + (legendLineHeight - text.Height) / 2;
                        dc.DrawText(new FormattedText("\t", System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Segoe UI"), fontSize * 0.8, Brushes.DimGray, dpi), new Point(legendX, sepY));
                        legendX += sepWidth;
                    }
                }
                legendY += legendLineHeight + legendSpacing * 0.2;
            }

            // --- Offset plot below header block (metrics + legend) ---
            double headerBlockHeight = metricsFt.Height + legendSpacing + legendLines.Count * (legendLineHeight + legendSpacing * 0.2);
            double plotYOffset = headerTop + headerBlockHeight + legendSpacing;

            // --- Draw rectangles for each peptide, with jitter for overlaps ---
            foreach (var res in filteredResults)
            {
                int startIdx = res.Start - 1;
                int endIdx = res.End - 1;
                int len = endIdx - startIdx + 1;

                int idx = startIdx;
                int remaining = len;
                int row = idx / lettersPerRow;
                int col = idx % lettersPerRow;

                bool isFirstSegment = true;
                while (remaining > 0)
                {
                    int drawLen = Math.Min(remaining, lettersPerRow - col);
                    int thisStartCol = col;
                    int thisEndCol = col + drawLen - 1;

                    // Keep rects if you need them elsewhere (optional)
                    if (!rowRectangles.TryGetValue(row, out var rects))
                    {
                        rects = new List<(int, int)>();
                        rowRectangles[row] = rects;
                    }
                    rects.Add((thisStartCol, thisEndCol));

                    // Per-row, per-column occupancy: first plotted at a column gets track 0.
                    // Only overlapping residues force indent > 0.
                    if (!rowOccupancy.TryGetValue(row, out var occ))
                    {
                        occ = new int[lettersPerRow]; // defaults to 0
                        rowOccupancy[row] = occ;
                    }

                    // Find the lowest track that can host this span = current max occupancy across span
                    int trackIndex = 0;
                    for (int c = thisStartCol; c <= thisEndCol && c < lettersPerRow; c++)
                        trackIndex = Math.Max(trackIndex, occ[c]);

                    // Reserve that track across the span
                    for (int c = thisStartCol; c <= thisEndCol && c < lettersPerRow; c++)
                        occ[c] = trackIndex + 1;

                    // Now draw using this track’s vertical offset
                    double x = plotMargin + thisStartCol * letterWidth;
                    double y = plotYOffset
                               + row * (letterHeight + yPad + rectHeight * 2)
                               + letterHeight + yPad
                               + trackIndex * rectJitter;
                    double width = drawLen * letterWidth;
                    double height = rectHeight;

                    // Determine which corners to round
                    double r = 6;
                    bool roundLeft = false, roundRight = false;

                    if (isFirstSegment)
                        roundLeft = true;
                    if (remaining == drawLen)
                        roundRight = true;

                    // Use a slightly opaque brush for fill
                    var baseBrush = GetColorBrush(res);
                    var color = (baseBrush as SolidColorBrush)?.Color ?? Colors.Gray;
                    var brush = new SolidColorBrush(color) { Opacity = 0.75 };

                    // Create a fully opaque, darkened pen for the outline
                    Color outlineColor = Color.Multiply(color, 0.9f); // darken for contrast
                    var outlinePen = new Pen(new SolidColorBrush(outlineColor) { Opacity = 1.0 }, 2.0);
                    outlinePen.Freeze();

                    if (roundLeft || roundRight)
                    {
                        // Custom geometry for per-corner rounding
                        var geometry = new StreamGeometry();
                        using (var ctx = geometry.Open())
                        {
                            // Top left
                            if (roundLeft)
                                ctx.BeginFigure(new Point(x + r, y), true, true);
                            else
                                ctx.BeginFigure(new Point(x, y), true, true);

                            // Top edge
                            ctx.LineTo(new Point(x + width - (roundRight ? r : 0), y), true, false);

                            // Top right
                            if (roundRight)
                                ctx.ArcTo(new Point(x + width, y + r), new Size(r, r), 0, false, SweepDirection.Clockwise, true, false);
                            else
                                ctx.LineTo(new Point(x + width, y), true, false);

                            // Right edge
                            ctx.LineTo(new Point(x + width, y + height - (roundRight ? r : 0)), true, false);

                            // Bottom right
                            if (roundRight)
                                ctx.ArcTo(new Point(x + width - r, y + height), new Size(r, r), 0, false, SweepDirection.Clockwise, true, false);
                            else
                                ctx.LineTo(new Point(x + width, y + height), true, false);

                            // Bottom edge
                            ctx.LineTo(new Point(x + (roundLeft ? r : 0), y + height), true, false);

                            // Bottom left
                            if (roundLeft)
                                ctx.ArcTo(new Point(x, y + height - r), new Size(r, r), 0, false, SweepDirection.Clockwise, true, false);
                            else
                                ctx.LineTo(new Point(x, y + height), true, false);

                            // Left edge
                            ctx.LineTo(new Point(x, y + (roundLeft ? r : 0)), true, false);

                            // Top left
                            if (roundLeft)
                                ctx.ArcTo(new Point(x + r, y), new Size(r, r), 0, false, SweepDirection.Clockwise, true, false);
                            else
                                ctx.LineTo(new Point(x, y), true, false);
                        }
                        geometry.Freeze();
                        dc.DrawGeometry(brush, outlinePen, geometry);
                    }
                    else
                    {
                        // No rounded corners, no outline
                        dc.DrawRectangle(brush, null, new Rect(x, y, width, height));
                    }

                    remaining -= drawLen;
                    idx += drawLen;
                    row++;
                    col = 0;
                    isFirstSegment = false;
                }
            }

            // --- Draw letters, bold and evenly spaced, with margin ---
            // Add row index labels on both sides
            // Draw the letters for this row and determine the vertical alignment
            for (int row = 0; row < nRows; row++)
            {
                int rowStartIdx = row * lettersPerRow;
                int rowEndIdx = Math.Min(seq.Length, (row + 1) * lettersPerRow) - 1;
                double y = plotYOffset + row * (letterHeight + yPad + rectHeight * 2);

                // Prepare left/right labels
                var leftLabelText = (rowStartIdx + 1).ToString();
                var rightLabelText = (rowEndIdx + 1).ToString();
                var labelFontSize = Math.Max(6, fontSize - 4);

                var leftLabel = new FormattedText(
                    leftLabelText,
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface("Segoe UI"),
                    labelFontSize,
                    Brushes.DimGray,
                    dpi);

                var rightLabel = new FormattedText(
                    rightLabelText,
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface("Segoe UI"),
                    labelFontSize,
                    Brushes.DimGray,
                    dpi);

                // Use the first letter in the row for vertical alignment
                int firstIdx = rowStartIdx;
                var letterText = new FormattedText(
                    seq[firstIdx].ToString(),
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface(new FontFamily("Segoe UI"), FontStyles.Normal, FontWeights.Bold, FontStretches.Normal),
                    fontSize,
                    Brushes.Black,
                    dpi);

                // Align baselines: offset so that the baseline of the label matches the baseline of the letter
                double baselineOffset = y + letterText.Baseline - leftLabel.Baseline;

                double leftLabelX = plotMargin - leftLabel.Width - 6;
                double rightLabelX = plotMargin + usableWidth + 6;

                dc.DrawText(leftLabel, new Point(leftLabelX, baselineOffset));
                dc.DrawText(rightLabel, new Point(rightLabelX, baselineOffset));

                // Draw the letters for this row
                for (int col = 0; col < lettersPerRow; col++)
                {
                    int i = rowStartIdx + col;
                    if (i >= seq.Length)
                        break;

                    double x = plotMargin + col * letterWidth;
                    var formattedText = new FormattedText(
                        seq[i].ToString(),
                        System.Globalization.CultureInfo.CurrentCulture,
                        FlowDirection.LeftToRight,
                        new Typeface(new FontFamily("Segoe UI"), FontStyles.Normal, FontWeights.Bold, FontStretches.Normal),
                        fontSize,
                        Brushes.Black,
                        dpi);

                    double centeredX = x + (letterWidth - formattedText.Width) / 2;
                    dc.DrawText(formattedText, new Point(centeredX, y));
                }
            }
        }

        CoverageDrawing = new DrawingImage(dv.Drawing);
    }

    // Preferred legend order
    private static readonly BioPolymerCoverageType[] LegendOrder =
    {
        BioPolymerCoverageType.Unique,
        BioPolymerCoverageType.UniqueMissedCleavage,
        BioPolymerCoverageType.TandemRepeat,
        BioPolymerCoverageType.TandemRepeatMissedCleavage,
        BioPolymerCoverageType.Shared,
        BioPolymerCoverageType.SharedMissedCleavage
    };

    private List<(SolidColorBrush Brush, FormattedText Text)> CreateLegendItems( List<BioPolymerCoverageResultModel> filteredResults, double fontSize, double dpi, ColorResultsBy colorBy)
    {
        var items = new List<(SolidColorBrush, FormattedText)>();
        if (colorBy == ColorResultsBy.None)
            return items;

        if (colorBy == ColorResultsBy.FileOrigin)
        {
            var countsByFile = filteredResults
                .GroupBy(r => r.Match.FileName)
                .ToDictionary(g => g.Key, g => g.Count());

            string unitLabel = filteredResults.FirstOrDefault()?.Match.GetDigestionProductLabel();

            foreach (var kvp in countsByFile)
            {
                string fileName = kvp.Key;
                int count = kvp.Value;
                if (string.IsNullOrEmpty(fileName) || count == 0) continue;

                string label = $"{fileName} {unitLabel}s: {count}";
                var ft = new FormattedText(
                    label,
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface("Segoe UI"),
                    fontSize * 0.8,
                    Brushes.Black,
                    dpi);
                items.Add((GetColorByIdentifier(fileName), ft));
            }
            return items;
        }

        var countsByType = filteredResults.GroupBy(r => r.CoverageType)
            .ToDictionary(g => g.Key, g => g.Count());

        string unitLabelType = filteredResults.FirstOrDefault()?.Match.GetDigestionProductLabel();

        foreach (var t in LegendOrder)
        {
            if (!countsByType.TryGetValue(t, out int count) || count == 0) continue;

            string label = $"{BaseViewModel.AddSpaces(t.ToString())} {unitLabelType}s: {count}";
            var ft = new FormattedText(
                label,
                System.Globalization.CultureInfo.CurrentCulture,
                FlowDirection.LeftToRight,
                new Typeface("Segoe UI"),
                fontSize * 0.8,
                Brushes.Black,
                dpi);
            items.Add((MetaDrawSettings.BioPolymerCoverageColors[t], ft));
        }
        return items;
    }

    // Helper: wrap legend items into lines that fit within usableWidth
    private List<List<(SolidColorBrush Brush, FormattedText Text)>> WrapLegendItems(
        List<(SolidColorBrush Brush, FormattedText Text)> items,
        double fontSize,
        double dpi,
        double usableWidth)
    {
        var lines = new List<List<(SolidColorBrush Brush, FormattedText Text)>>();
        var currentLine = new List<(SolidColorBrush Brush, FormattedText Text)>();
        double sepWidth = new FormattedText("\t", System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Segoe UI"), fontSize * 0.8, Brushes.DimGray, dpi).Width;
        double currentWidth = 0;

        foreach (var item in items)
        {
            double itemWidth = fontSize * 0.8 + 6 + item.Text.Width + (currentLine.Count > 0 ? sepWidth : 0);
            if (currentLine.Count > 0 && currentWidth + itemWidth > usableWidth)
            {
                lines.Add(currentLine);
                currentLine = new List<(SolidColorBrush Brush, FormattedText Text)>();
                currentWidth = 0;
                itemWidth = fontSize * 0.8 + 6 + item.Text.Width;
            }
            currentLine.Add(item);
            currentWidth += itemWidth;
        }
        if (currentLine.Count > 0)
            lines.Add(currentLine);

        return lines;
    }
}