using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public class BioPolymerCoverageMapViewModel : BaseViewModel
{

    private DrawingImage _coverageDrawing;
    public DrawingImage CoverageDrawing
    {
        get => _coverageDrawing;
        private set { _coverageDrawing = value; OnPropertyChanged(nameof(CoverageDrawing)); }
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
        int size;
        try
        {
            size = MetaDrawSettings.BioPolymerCoverageFontSize;
        }
        catch
        {
            size = 16; // default
        }

        _availableWidth = availableWidth;
        int newLetters = Math.Max(10, (int)(availableWidth / (size * 0.70)));
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

        double plotMargin = fontSize * 2; // margin on left and right
        double usableWidth = Math.Max(1, AvailableWidth - 2 * plotMargin);
        double letterWidth = usableWidth / lettersPerRow;
        double letterHeight = fontSize * 1.2;
        double rectHeight = fontSize * 0.8;
        double yPad = fontSize * 0.2;
        double rectJitter = rectHeight * 0.5;

        int nRows = (seq.Length + lettersPerRow - 1) / lettersPerRow;
        var rowRectangles = new Dictionary<int, List<(int startCol, int endCol)>>();

        var typesInUse = results.Select(r => r.CoverageType).Distinct().ToList();

        var dv = new DrawingVisual();
        using (var dc = dv.RenderOpen())
        {
            // --- Header metrics + legend (centered) ---
            double legendBoxSize = fontSize * 0.8;
            double legendSpacing = fontSize * 0.5;
            double dpi = MetaDrawSettings.CanvasPdfExportDpi;
            double headerTop = legendSpacing;

            // === Metrics: Unique Coverage / Maximum Coverage (counts only) ===
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

            var metricsText = $"Unique Coverage: {uniqueCoveredCount}/{seqLen} ({Math.Round(100 * uniqueCoveredCount / (double)seqLen, 2)}%)   |   Maximum Coverage: {maxCoveredCount}/{seqLen} ({Math.Round(100 * maxCoveredCount / (double)seqLen ,2)}%)";
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

            // === Legend ===
            // 1) counts by type from the filtered (visible) set
            var countsByType = filteredResults.GroupBy(r => r.CoverageType)
                                              .ToDictionary(g => g.Key, g => g.Count());

            // 2) Oligos/Peptides/Matches label
            string unitLabel = GuessUnitLabel(filteredResults.FirstOrDefault()?.Match);

            // 3) Build ordered legend items that actually exist
            var legendItems = new List<(BioPolymerCoverageType Type, FormattedText Text)>();
            foreach (var t in LegendOrder)
            {
                if (!countsByType.TryGetValue(t, out int count) || count == 0) continue;
                string label = $"{HumanizeCoverageType(t)} {unitLabel}: {count}";
                var ft = new FormattedText(
                    label,
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface("Segoe UI"),
                    fontSize * 0.8,
                    Brushes.Black,
                    dpi);
                legendItems.Add((t, ft));
            }

            double legendLineHeight = legendItems.Count == 0
                ? 0
                : Math.Max(legendBoxSize, legendItems.Max(li => li.Text.Height));

            // Separator " | "
            var sepFt = new FormattedText(
                "\t",
                System.Globalization.CultureInfo.CurrentCulture,
                FlowDirection.LeftToRight,
                new Typeface("Segoe UI"),
                fontSize * 0.8,
                Brushes.DimGray,
                dpi);

            // Total legend width
            double legendTotalWidth = 0;
            for (int i = 0; i < legendItems.Count; i++)
            {
                legendTotalWidth += legendBoxSize + 6 + legendItems[i].Text.Width;
                if (i < legendItems.Count - 1) legendTotalWidth += sepFt.Width;
            }

            double legendY = metricsY + metricsFt.Height + legendSpacing;
            double legendX = plotMargin + (usableWidth - legendTotalWidth) / 2;

            // Draw legend
            for (int i = 0; i < legendItems.Count; i++)
            {
                var (type, text) = legendItems[i];
                var brush = MetaDrawSettings.BioPolymerCoverageColors[type];

                // color square
                var rect = new Rect(legendX, legendY + (legendLineHeight - legendBoxSize) / 2, legendBoxSize, legendBoxSize);
                dc.DrawRoundedRectangle(brush, null, rect, 4, 4);
                legendX += legendBoxSize + 6;

                // label
                double textY = legendY + (legendLineHeight - text.Height) / 2;
                dc.DrawText(text, new Point(legendX, textY));
                legendX += text.Width;

                if (i < legendItems.Count - 1)
                {
                    double sepY = legendY + (legendLineHeight - sepFt.Height) / 2;
                    dc.DrawText(sepFt, new Point(legendX, sepY));
                    legendX += sepFt.Width;
                }
            }

            // --- Offset plot below header block (metrics + optional legend) ---
            double headerBlockHeight = metricsFt.Height + (legendItems.Count > 0 ? (legendSpacing + legendLineHeight) : 0);
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

                    if (!rowRectangles.TryGetValue(row, out var rects))
                    {
                        rects = new List<(int, int)>();
                        rowRectangles[row] = rects;
                    }
                    int overlap = rects.Count(r => !(thisEndCol < r.startCol || thisStartCol > r.endCol));
                    rects.Add((thisStartCol, thisEndCol));

                    double x = plotMargin + thisStartCol * letterWidth;
                    double y = plotYOffset + row * (letterHeight + yPad + rectHeight * 2) + letterHeight + yPad + overlap * rectJitter;
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
                    var baseBrush = MetaDrawSettings.BioPolymerCoverageColors[res.CoverageType];
                    var color = (baseBrush as SolidColorBrush)?.Color ?? Colors.Gray;
                    var brush = new SolidColorBrush(color) { Opacity = 0.75 };

                    // Create a fully opaque, darkened pen for the outline
                    Color outlineColor = Color.Multiply(color, 0.6f); // darken for contrast
                    var outlinePen = new Pen(new SolidColorBrush(outlineColor) { Opacity = 1.0 }, 1.0);
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

    private static string HumanizeCoverageType(BioPolymerCoverageType t) => t switch
    {
        BioPolymerCoverageType.Unique => "Unique",
        BioPolymerCoverageType.UniqueMissedCleavage => "Unique Missed Cleavages",
        BioPolymerCoverageType.TandemRepeat => "Tandem Repeat",
        BioPolymerCoverageType.TandemRepeatMissedCleavage => "Tandem Repeat Missed Cleavages",
        BioPolymerCoverageType.Shared => "Shared",
        BioPolymerCoverageType.SharedMissedCleavage => "Shared Missed Cleavages",
        _ => t.ToString()
    };

    // Safer than hard-referencing OsmFromTsv / PsmFromTsv types:
    private static string GuessUnitLabel(object match)
    {
        if (match is null) return "Matches";
        var n = match.GetType().Name;
        if (n.IndexOf("Osm", StringComparison.OrdinalIgnoreCase) >= 0) return "Oligos";
        if (n.IndexOf("Psm", StringComparison.OrdinalIgnoreCase) >= 0) return "Peptides";
        return "Matches";
    }
}
