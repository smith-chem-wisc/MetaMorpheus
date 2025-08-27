using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Windows.Media;
using System.Windows;

namespace GuiFunctions;

public class BioPolymerCoverageMapViewModel : BaseViewModel
{
    private BioPolymerGroupViewModel _group;
    private int _lettersPerRow = 50;
    private double _fontSize = 16;
    private Dictionary<BioPolymerCoverageType, Brush> _coverageTypeBrushes;

    public BioPolymerCoverageMapViewModel()
    {
        // Default colors, can be exposed as settings
        _coverageTypeBrushes = new()
        {
            { BioPolymerCoverageType.Unique, Brushes.LightGreen },
            { BioPolymerCoverageType.UniqueMissedCleavage, Brushes.YellowGreen },
            { BioPolymerCoverageType.TandemRepeat, Brushes.LightBlue },
            { BioPolymerCoverageType.TandemRepeatMissedCleavage, Brushes.SkyBlue },
            { BioPolymerCoverageType.Shared, Brushes.Orange },
            { BioPolymerCoverageType.SharedMissedCleavage, Brushes.OrangeRed }
        };
    }

    private DrawingImage _coverageDrawing;
    public DrawingImage CoverageDrawing
    {
        get => _coverageDrawing;
        private set { _coverageDrawing = value; OnPropertyChanged(nameof(CoverageDrawing)); }
    }

    public BioPolymerGroupViewModel Group
    {
        get => _group;
        set { _group = value; OnPropertyChanged(nameof(Group)); Redraw(); }
    }

    public int LettersPerRow
    {
        get => _lettersPerRow;
        set { _lettersPerRow = value; OnPropertyChanged(nameof(LettersPerRow)); }
    }

    public double FontSize
    {
        get => _fontSize;
        set { _fontSize = value; OnPropertyChanged(nameof(FontSize)); }
    }

    public Dictionary<BioPolymerCoverageType, Brush> CoverageTypeBrushes => _coverageTypeBrushes;

    // Call this from the view when the canvas size changes
    public void UpdateLettersPerRow(double availableWidth)
    {
        // Estimate: each letter ~fontSize*0.7 wide, add some padding
        int newLetters = Math.Max(10, (int)(availableWidth / (FontSize * 0.7)));
        if (newLetters != LettersPerRow)
            LettersPerRow = newLetters;
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
        double fontSize = FontSize;
        double letterWidth = fontSize * 0.7;
        double letterHeight = fontSize * 1.2;
        double rectHeight = fontSize * 0.8;
        double yPad = fontSize * 0.2;

        var dv = new DrawingVisual();
        using (var dc = dv.RenderOpen())
        {
            // Draw rectangles for each peptide
            foreach (var res in results)
            {
                int startIdx = res.Start - 1;
                int endIdx = res.End - 1;
                int row = startIdx / lettersPerRow;
                int col = startIdx % lettersPerRow;
                int len = endIdx - startIdx + 1;

                double x = col * letterWidth;
                double y = row * (letterHeight + yPad) + letterHeight * 0.2;
                double width = Math.Min(len, lettersPerRow - col) * letterWidth;
                double height = rectHeight;

                var brush = CoverageTypeBrushes[res.CoverageType];
                var rect = new Rect(x, y, width, height);
                dc.DrawRoundedRectangle(brush, null, rect, 6, 6);

                // If peptide spans multiple rows, draw additional rectangles
                int remaining = len - (lettersPerRow - col);
                int nextStart = startIdx + (lettersPerRow - col);
                while (remaining > 0)
                {
                    int nextRow = nextStart / lettersPerRow;
                    int nextCol = nextStart % lettersPerRow;
                    int drawLen = Math.Min(remaining, lettersPerRow);

                    double nx = nextCol * letterWidth;
                    double ny = nextRow * (letterHeight + yPad) + letterHeight * 0.2;
                    double nwidth = drawLen * letterWidth;

                    var nrect = new Rect(nx, ny, nwidth, height);
                    dc.DrawRoundedRectangle(brush, null, nrect, 6, 6);

                    remaining -= drawLen;
                    nextStart += drawLen;
                }
            }

            // Draw letters
            for (int i = 0; i < seq.Length; i++)
            {
                int row = i / lettersPerRow;
                int col = i % lettersPerRow;
                double x = col * letterWidth;
                double y = row * (letterHeight + yPad);

                var formattedText = new FormattedText(
                    seq[i].ToString(),
                    System.Globalization.CultureInfo.CurrentCulture,
                    FlowDirection.LeftToRight,
                    new Typeface("Segoe UI"),
                    fontSize,
                    Brushes.Black,
                    VisualTreeHelper.GetDpi(Application.Current.MainWindow).PixelsPerDip);

                dc.DrawText(formattedText, new Point(x, y));
            }
        }

        CoverageDrawing = new DrawingImage(dv.Drawing);
    }
}
