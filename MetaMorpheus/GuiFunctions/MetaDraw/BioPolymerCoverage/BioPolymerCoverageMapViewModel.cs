using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Media;
using System.Windows;

namespace GuiFunctions;

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
        int newLetters = Math.Max(10, (int)(availableWidth / (size * 0.8)));
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
            // --- Legend at the top, centered horizontally ---
            double legendBoxSize = fontSize * 0.8;
            double legendSpacing = fontSize * 0.5;
            double legendY = legendSpacing;

            // Calculate legend total width
            double legendTotalWidth = 0;
            var legendEntries = new List<(BioPolymerCoverageType type, double textWidth, FormattedText text)>();
            foreach (var type in typesInUse)
            {
                var text = AddSpaces(type.ToString());
                var formattedText = new FormattedText(
                    text,
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface("Segoe UI"),
                    fontSize * 0.8,
                    Brushes.Black,
                    VisualTreeHelper.GetDpi(Application.Current.MainWindow).PixelsPerDip);

                legendEntries.Add((type, formattedText.Width, formattedText));
                legendTotalWidth += legendBoxSize + 6 + formattedText.Width + 24;
            }
            if (legendTotalWidth > 0)
                legendTotalWidth -= 24; // Remove last extra spacing

            double legendX = plotMargin + (usableWidth - legendTotalWidth) / 2;

            // Draw the legend
            foreach (var (type, textWidth, formattedText) in legendEntries)
            {
                var brush = MetaDrawSettings.BioPolymerCoverageColors[type];
                var rect = new Rect(legendX, legendY, legendBoxSize, legendBoxSize);
                dc.DrawRoundedRectangle(brush, null, rect, 4, 4);

                // Vertically center the text with respect to the color box
                double textY = legendY + (legendBoxSize - formattedText.Height) / 2;
                dc.DrawText(formattedText, new Point(legendX + legendBoxSize + 6, textY));
                legendX += legendBoxSize + 6 + textWidth + 24;
            }

            // --- Offset plot below the legend ---
            double plotYOffset = legendBoxSize + legendSpacing * 2;

            // --- Draw rectangles for each peptide, with jitter for overlaps ---
            foreach (var res in results.Where(p => MetaDrawSettings.FilterAcceptsPsm(p.Match)))
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
            for (int i = 0; i < seq.Length; i++)
            {
                int row = i / lettersPerRow;
                int col = i % lettersPerRow;
                double x = plotMargin + col * letterWidth;
                double y = plotYOffset + row * (letterHeight + yPad + rectHeight * 2);

                var formattedText = new FormattedText(
                    seq[i].ToString(),
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface(new FontFamily("Segoe UI"), FontStyles.Normal, FontWeights.Bold, FontStretches.Normal),
                    fontSize,
                    Brushes.Black,
                    VisualTreeHelper.GetDpi(Application.Current.MainWindow).PixelsPerDip);

                double centeredX = x + (letterWidth - formattedText.Width) / 2;
                dc.DrawText(formattedText, new Point(centeredX, y));
            }
        }

        CoverageDrawing = new DrawingImage(dv.Drawing);
    }
}
