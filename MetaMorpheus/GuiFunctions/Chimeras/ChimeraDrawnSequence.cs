using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;
using Omics.Fragmentation;

namespace GuiFunctions;

public class ChimeraDrawnSequence
{
    public Canvas SequenceDrawingCanvas;
    public ChimeraGroupViewModel ChimeraGroupViewModel;
    private int _numSequences;
    private static readonly int _yStep = 40;
    private static readonly int _canvasBuffer = 20;
    private ChimeraAnalysisTabViewModel? _parent;

    // Pools for UI elements
    private readonly List<TextBlock> _textBlockPool = new();
    private int _textBlockPoolIndex = 0;
    private readonly List<Ellipse> _ellipsePool = new();
    private int _ellipsePoolIndex = 0;
    private readonly List<Polyline> _polyLinePool = new();
    private int _polyLinePoolIndex = 0;

    public ChimeraDrawnSequence(Canvas sequenceDrawingCanvas, ChimeraGroupViewModel chimeraGroupViewModel,
        ChimeraAnalysisTabViewModel parent = null)
    {
        SequenceDrawingCanvas = sequenceDrawingCanvas;
        ChimeraGroupViewModel = chimeraGroupViewModel;
        _numSequences = ChimeraGroupViewModel.ChimericPsms.Count;
        _parent = parent;

        Draw();
    }

    public ChimeraDrawnSequence UpdateData(ChimeraGroupViewModel newGroup)
    {
        ChimeraGroupViewModel = newGroup;
        _numSequences = ChimeraGroupViewModel.ChimericPsms.Count;
        Draw();

        return this;
    }

    private void Draw()
    {
        DrawnSequence.ClearCanvas(SequenceDrawingCanvas);
        _textBlockPoolIndex = 0;
        _ellipsePoolIndex = 0;
        _polyLinePoolIndex = 0;

        // Set dimensions
        var longestSequenceLength = ChimeraGroupViewModel.ChimericPsms.Max(psm => psm.Psm.BaseSeq.Split('|')[0].Length);
        SequenceDrawingCanvas.Width = (longestSequenceLength + 4) * MetaDrawSettings.AnnotatedSequenceTextSpacing + _canvasBuffer;
        SequenceDrawingCanvas.Height = _yStep * _numSequences + _canvasBuffer;

        DrawSequences();

        // Remove unused pooled elements from the canvas
        for (int i = _textBlockPoolIndex; i < _textBlockPool.Count; i++)
        {
            if (_textBlockPool[i].Parent == SequenceDrawingCanvas)
                SequenceDrawingCanvas.Children.Remove(_textBlockPool[i]);
        }
        for (int i = _ellipsePoolIndex; i < _ellipsePool.Count; i++)
        {
            if (_ellipsePool[i].Parent == SequenceDrawingCanvas)
                SequenceDrawingCanvas.Children.Remove(_ellipsePool[i]);
        }
        for (int i = _polyLinePoolIndex; i < _polyLinePool.Count; i++)
        {
            if (_polyLinePool[i].Parent == SequenceDrawingCanvas)
                SequenceDrawingCanvas.Children.Remove(_polyLinePool[i]);
        }
    }

    private void DrawSequences()
    {
        int maxBaseSeqLength = ChimeraGroupViewModel.ChimericPsms.Max(p => p.Psm.BaseSeq.Length);
        for (var index = 0; index < ChimeraGroupViewModel.ChimericPsms.Count; index++)
        {
            var psm = ChimeraGroupViewModel.ChimericPsms[index];
            DrawBaseSequence(psm, index);
            AddModifications(psm, index);
            AddMatchedIons(psm, index);
            AddCircles(psm, index, maxBaseSeqLength);
        }
    }

    private void AddMatchedIons(ChimericSpectralMatchModel psm, int row)
    {
        Color color;

        // Shared Ions
        if (ChimeraGroupViewModel.MatchedFragmentIonsByColor.ContainsKey(ChimeraSpectrumMatchPlot.MultipleProteinSharedColor))
        {
            color = DrawnSequence.ParseColorFromOxyColor(ChimeraSpectrumMatchPlot.MultipleProteinSharedColor);
            foreach (var ion in ChimeraGroupViewModel
                         .MatchedFragmentIonsByColor[ChimeraSpectrumMatchPlot.MultipleProteinSharedColor]
                         .Select(p => p.Item1))
            {
                var matchedIon = psm.Psm.MatchedIons.FirstOrDefault(mi =>
                    mi.Charge == ion.Charge &&
                    mi.NeutralTheoreticalProduct.ProductType == ion.NeutralTheoreticalProduct.ProductType &&
                    mi.NeutralTheoreticalProduct.FragmentNumber == ion.NeutralTheoreticalProduct.FragmentNumber &&
                    Math.Abs(mi.Intensity - ion.Intensity) < 1e-6 &&
                    Math.Abs(mi.Mz - ion.Mz) < 1e-6);

                if (matchedIon == null)
                    continue;
                AddMatchedIon(ion, color, row, psm.Psm.BaseSeq.Length);
            }
        }

        // Protein Shared Ions
        if (ChimeraGroupViewModel.MatchedFragmentIonsByColor.ContainsKey(psm.ProteinColor))
        {
            color = DrawnSequence.ParseColorFromOxyColor(psm.ProteinColor);
            foreach (var ion in ChimeraGroupViewModel
                         .MatchedFragmentIonsByColor[psm.ProteinColor]
                         .Select(p => p.Item1))
            {
                if (!psm.Psm.MatchedIons.Contains(ion)) continue;
                AddMatchedIon(ion, color, row, psm.Psm.BaseSeq.Length);
            }
        }

        // Unique Ions
        if (ChimeraGroupViewModel.MatchedFragmentIonsByColor.ContainsKey(psm.Color))
        {
            color = DrawnSequence.ParseColorFromOxyColor(psm.Color);
            foreach (var ion in ChimeraGroupViewModel
                         .MatchedFragmentIonsByColor[psm.Color]
                         .Select(p => p.Item1))
            {
                AddMatchedIon(ion, color, row, psm.Psm.BaseSeq.Length);
            }
        }
    }

    private void AddMatchedIon(MatchedFragmentIon ion, Color color, int row, int sequenceLength)
    {
        if (ion.IsInternalFragment) 
            return;

        double x, y;
        var residueNum = ion.NeutralTheoreticalProduct.ProductType == ProductType.y
            ? sequenceLength - ion.NeutralTheoreticalProduct.FragmentNumber
            : ion.NeutralTheoreticalProduct.AminoAcidPosition;
        x = GetX(residueNum);
        y = GetY(row) + MetaDrawSettings.ProductTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType];

        
        if (ion.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime)
        {
            DrawCTermIonPooled(new Point(x, y), color,2);
        }
        else if (ion.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime)
        {
            DrawNTermIonPooled(new Point(x, y), color, 2);
        }
    }

    private void AddCircles(ChimericSpectralMatchModel psm, int row, int maxBaseSeqLength)
    {
        var color = DrawnSequence.ParseColorBrushFromOxyColor(psm.Color);
        DrawCirclePooled(new Point(GetX(maxBaseSeqLength + 1), GetY(row)), color);
        DrawCirclePooled(new Point(GetX(-2), GetY(row)), color);
    }

    private void DrawBaseSequence(ChimericSpectralMatchModel psm, int row)
    {
        var baseSeq = psm.Psm.BaseSeq.Split('|')[0];
        int index = 0;
        for (; index < baseSeq.Length; index++)
        {
            var x = GetX(index);
            var y = GetY(row);
            DrawTextPooled(new Point(x, y), baseSeq[index].ToString(), Brushes.Black);
        }
    }

    private void AddModifications(ChimericSpectralMatchModel psm, int row)
    {
        foreach (var mod in psm.AllModsOneIsNterminus)
        {
            var x = GetX(mod.Key - 2);
            var y = GetY(row);
            var color = DrawnSequence.ParseColorBrushFromOxyColor(psm.Color);
            DrawCirclePooled(new Point(x, y), color);
        }
    }

    private static double GetY(int row)
    {
        return (row * _yStep) + 10;
    }

    private static double GetX(int residueIndex)
    {
        return (residueIndex + 1) * MetaDrawSettings.AnnotatedSequenceTextSpacing + 22;
    }

    private void DrawTextPooled(Point loc, string txt, Brush clr)
    {
        TextBlock tb;
        if (_textBlockPoolIndex < _textBlockPool.Count)
        {
            tb = _textBlockPool[_textBlockPoolIndex];
        }
        else
        {
            tb = new TextBlock();
            _textBlockPool.Add(tb);
        }
        _textBlockPoolIndex++;

        tb.Foreground = clr;
        tb.Text = txt;
        tb.Height = 30;
        tb.FontSize = 25;
        tb.FontWeight = FontWeights.Bold;
        tb.FontFamily = new FontFamily("Arial");
        tb.TextAlignment = TextAlignment.Center;
        tb.HorizontalAlignment = HorizontalAlignment.Center;
        tb.Width = 24;

        Canvas.SetTop(tb, loc.Y);
        Canvas.SetLeft(tb, loc.X);
        Panel.SetZIndex(tb, 2);

        if (tb.Parent != SequenceDrawingCanvas)
            SequenceDrawingCanvas.Children.Add(tb);
    }

    private void DrawCirclePooled(Point loc, SolidColorBrush clr)
    {
        Ellipse circle;
        if (_ellipsePoolIndex < _ellipsePool.Count)
        {
            circle = _ellipsePool[_ellipsePoolIndex];
        }
        else
        {
            circle = new Ellipse();
            _ellipsePool.Add(circle);
        }
        _ellipsePoolIndex++;

        circle.Width = 24;
        circle.Height = 24;
        circle.Stroke = clr;
        circle.StrokeThickness = 1;
        circle.Fill = clr;
        circle.Opacity = 0.7;

        Canvas.SetLeft(circle, loc.X);
        Canvas.SetTop(circle, loc.Y);
        Panel.SetZIndex(circle, 1);

        if (circle.Parent != SequenceDrawingCanvas)
            SequenceDrawingCanvas.Children.Add(circle);
    }
    private void DrawCTermIonPooled(Point topLoc, Color clr, int thickness = 1)
    {
        Polyline bot;
        if (_polyLinePoolIndex < _polyLinePool.Count)
        {
            bot = _polyLinePool[_polyLinePoolIndex];
        }
        else
        {
            bot = new Polyline();
            _polyLinePool.Add(bot);
        }
        _polyLinePoolIndex++;


        double x = topLoc.X, y = topLoc.Y;
        bot.Points = new PointCollection() { new Point(x + 10, y + 10), new Point(x, y + 10), new Point(x, y + 24) };
        bot.Stroke = new SolidColorBrush(clr);
        bot.StrokeThickness = thickness;
        Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
        SequenceDrawingCanvas.Children.Add(bot);
    }

    private void DrawNTermIonPooled(Point botLoc, Color clr, int thickness = 1)
    {
        Polyline bot;
        if (_polyLinePoolIndex < _polyLinePool.Count)
        {
            bot = _polyLinePool[_polyLinePoolIndex];
        }
        else
        {
            bot = new Polyline();
            _polyLinePool.Add(bot);
        }
        _polyLinePoolIndex++;

        double x = botLoc.X, y = botLoc.Y;
        bot.Points = new PointCollection() { new Point(x - 10, y - 10), new Point(x, y - 10), new Point(x, y - 24) };
        bot.Stroke = new SolidColorBrush(clr);
        bot.StrokeThickness = thickness;
        Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
        SequenceDrawingCanvas.Children.Add(bot);
    }
}
