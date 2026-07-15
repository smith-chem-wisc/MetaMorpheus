using CsvHelper.Configuration.Attributes;
using GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;
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
    PrecursorIntensity,
    Score,
}
public static class ColorResultsByExtensions
{
    public static bool IsNumeric(this ColorResultsBy colorBy)
    {
        switch (colorBy)
        {
            case ColorResultsBy.PrecursorIntensity:
            case ColorResultsBy.Score:
                return true;
            default:
                return false;
        }
    }
}

public class BioPolymerCoverageMapViewModel : BaseViewModel
{
    #region Color Handling

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
                OnPropertyChanged(nameof(IsNumericColorMode));
                ApplyNumericModeDefaultsIfNeeded(value);
                Redraw();
            }
        }
    }

    public ColorGradientType[] AllGradients { get; } = Enum.GetValues<ColorGradientType>();
    public ColorGradientType SelectedGradientType
    {
        get => MetaDrawSettings.BioPolymerCoverageGradientType;
        set
        {
            if (MetaDrawSettings.BioPolymerCoverageGradientType != value)
            {
                MetaDrawSettings.BioPolymerCoverageGradientType = value;
                OnPropertyChanged(nameof(SelectedGradientType));
                Redraw();
            }
        }
    }

    private bool _useLogColorScale;
    public bool UseLogColorScale
    {
        get => _useLogColorScale;
        set
        {
            if (_useLogColorScale != value)
            {
                _useLogColorScale = value;
                OnPropertyChanged(nameof(UseLogColorScale));
                Redraw();
            }
        }
    }

    public bool IsNumericColorMode => ColorBy.IsNumeric();


    private Dictionary<ColorResultsBy, bool> _hasAppliedLogDefault = Enum.GetValues<ColorResultsBy>().ToDictionary(p => p, p => false);
    private void ApplyNumericModeDefaultsIfNeeded(ColorResultsBy color)
    {
        if (_hasAppliedLogDefault[color]) 
            return;

        if (!color.IsNumeric()) 
            return;

        var numeric = BioPolymerCoverageColorMapperFactory.Create(color);
        _useLogColorScale = numeric.DefaultUseLogScale;
        _hasAppliedLogDefault[color] = true;
        OnPropertyChanged(nameof(UseLogColorScale));
    }

    #endregion

    #region Drawing Properties

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

    #endregion

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
        var bioPolymerName = Group.ProteinName;
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

        var colorMap = BioPolymerCoverageColorMapperFactory.Create(_colorBy);
        colorMap.Prepare(filteredResults, SelectedGradientType, UseLogColorScale);

        var dv = new DrawingVisual();
        using (var dc = dv.RenderOpen())
        {
            double legendSpacing = fontSize * 0.5;
            double dpi = MetaDrawSettings.CanvasPdfExportDpi;
            double headerTop = legendSpacing;

            // --- Draw BioPolymer Name at the very top ---
            var bioPolymerNameText = new FormattedText(
                bioPolymerName,
                System.Globalization.CultureInfo.CurrentCulture,
                FlowDirection.LeftToRight,
                new Typeface("Segoe UI Bold"),
                fontSize * 1.1,
                Brushes.Black,
                dpi);

            double nameX = plotMargin + (usableWidth - bioPolymerNameText.Width) / 2;
            double nameY = legendSpacing; // Start at the top margin
            dc.DrawText(bioPolymerNameText, new Point(nameX, nameY));

            // --- Global coverage metrics ---
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
            double metricsY = nameY + bioPolymerNameText.Height + legendSpacing * 0.5;
            dc.DrawText(metricsFt, new Point(metricsX, metricsY));

            // --- Legend, wrapped ---
            var sep = new FormattedText("\t", System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Segoe UI"), fontSize * 0.8, Brushes.DimGray, dpi);
            var sepWidth = sep.Width;

            double legendY = metricsY + metricsFt.Height + legendSpacing;
            var legendItems = CreateLegendItems(fontSize, dpi, colorMap);

            var legendLines = WrapLegendItems(legendItems, fontSize, usableWidth, sepWidth);

            double legendLineHeight = legendItems.Count == 0 ? 0 : Math.Max(fontSize * 0.8, legendItems.Max(li => li.Text.Height));


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
                        dc.DrawText(sep, new Point(legendX, sepY));
                        legendX += sepWidth;
                    }
                }
                legendY += legendLineHeight + legendSpacing * 0.2;
            }

            double gradientBlockH = 0;
            if (colorMap.IsNumeric)
                gradientBlockH = DrawNumericLegend(dc, legendY, plotMargin, usableWidth, legendSpacing, fontSize, dpi, colorMap, out legendY);

            // --- Offset plot below header block (metrics + legend) ---
            double headerBlockHeight = bioPolymerNameText.Height + metricsFt.Height + legendSpacing
                + legendLines.Count * (legendLineHeight + legendSpacing * 0.2)
                + gradientBlockH;
            double plotYOffset = headerTop + headerBlockHeight + legendSpacing;
            var drawnRects = new HashSet<(int row, int startCol, int endCol)>();

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

                    // Cache drawn rectangles to avoid duplicates
                    var rectKey = (row, thisStartCol, thisEndCol);
                    if (!drawnRects.Add(rectKey))
                    {
                        // Skip duplicate rectangle
                        remaining -= drawLen;
                        idx += drawLen;
                        row++;
                        col = 0;
                        isFirstSegment = false;
                        continue;
                    }

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
                    var baseBrush = colorMap.GetBrush(res);
                    var color = (baseBrush as SolidColorBrush)?.Color ?? Colors.Gray;
                    var brush = new SolidColorBrush(color) { Opacity = 0.75 };
                    brush.Freeze();

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

    private double DrawNumericLegend(
        DrawingContext dc,
        double legendY,
        double plotMargin,
        double usableWidth,
        double legendSpacing,
        double fontSize,
        double dpi,
        BioPolymerCoverageColorMapper colorMapper,
        out double newLegendY)
    {
        newLegendY = legendY;
        var brushes = colorMapper.GradientBrushes;
        if (brushes is null || brushes.Count == 0)
            return 0;

        double barWidth = Math.Min(250, usableWidth * 0.6);
        double barHeight = fontSize * 0.8;
        double barX = plotMargin + (usableWidth - barWidth) / 2;
        double barY = legendY + legendSpacing;
        double bandW = barWidth / brushes.Count;
        for (int i = 0; i < brushes.Count; i++)
            dc.DrawRectangle(brushes[i], null, new Rect(barX + i * bandW, barY, bandW + 0.5, barHeight));

        var metricTxt = new FormattedText(colorMapper.GradientLegendTitle!, System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Segoe UI"), fontSize * 0.7, Brushes.DimGray, dpi);
        dc.DrawText(metricTxt, new Point(barX + (barWidth - metricTxt.Width) / 2, barY - metricTxt.Height - 2));

        var minTxt = new FormattedText($"{colorMapper.GradientMinValue:G3}", System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Segoe UI"), fontSize * 0.7, Brushes.Black, dpi);
        var maxTxt = new FormattedText($"{colorMapper.GradientMaxValue:G3}", System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Segoe UI"), fontSize * 0.7, Brushes.Black, dpi);
        dc.DrawText(minTxt, new Point(barX, barY + barHeight + 2));
        dc.DrawText(maxTxt, new Point(barX + barWidth - maxTxt.Width, barY + barHeight + 2));

        double blockH = barHeight + Math.Max(minTxt.Height, maxTxt.Height) + legendSpacing * 2;
        newLegendY = legendY + blockH;
        return blockH;
    }

    private List<(SolidColorBrush Brush, FormattedText Text)> CreateLegendItems(double fontSize, double dpi, BioPolymerCoverageColorMapper colorMapper)
    {
        var items = new List<(SolidColorBrush, FormattedText)>();
        foreach (var item in colorMapper.LegendItems)
        {
            var ft = new FormattedText(
                item.Label,
                System.Globalization.CultureInfo.CurrentCulture,
                FlowDirection.LeftToRight,
                new Typeface("Segoe UI"),
                fontSize * 0.8,
                Brushes.Black,
                dpi);
            items.Add((item.Brush, ft));
        }
        return items;
    }

    // Helper: wrap legend items into lines that fit within usableWidth
    private List<List<(SolidColorBrush Brush, FormattedText Text)>> WrapLegendItems(
        List<(SolidColorBrush Brush, FormattedText Text)> items,
        double fontSize,
        double usableWidth,
        double sepWidth)
    {
        var lines = new List<List<(SolidColorBrush Brush, FormattedText Text)>>();
        var currentLine = new List<(SolidColorBrush Brush, FormattedText Text)>();
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
