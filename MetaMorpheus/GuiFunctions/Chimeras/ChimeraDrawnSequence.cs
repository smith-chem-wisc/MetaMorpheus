using System;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using Omics.Fragmentation;

namespace GuiFunctions;

public class ChimeraDrawnSequence
{
    public Canvas SequenceDrawingCanvas;
    public ChimeraGroupViewModel ChimeraGroupViewModel;
    private readonly int _numSequences;
    private static readonly int _yStep = 40;
    private static readonly int _canvasBuffer = 20;
    private ChimeraAnalysisTabViewModel? _parent;

    public ChimeraDrawnSequence(Canvas sequenceDrawingCanvas, ChimeraGroupViewModel chimeraGroupViewModel,
        ChimeraAnalysisTabViewModel parent = null)
    {
        SequenceDrawingCanvas = sequenceDrawingCanvas;
        ChimeraGroupViewModel = chimeraGroupViewModel;
        _numSequences = ChimeraGroupViewModel.ChimericPsms.Count;
        _parent = parent;

        DrawnSequence.ClearCanvas(SequenceDrawingCanvas);
        SetDrawingDimensions();
        DrawSequences();
    }

    private void SetDrawingDimensions()
    {
        var longestSequenceLength = ChimeraGroupViewModel.ChimericPsms.Max(psm => psm.Psm.BaseSeq.Split('|')[0].Length);
        SequenceDrawingCanvas.Width = (longestSequenceLength + 4) * MetaDrawSettings.AnnotatedSequenceTextSpacing + _canvasBuffer;
        SequenceDrawingCanvas.Height = _yStep * _numSequences + _canvasBuffer;
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
        double x, y;
        var residueNum = ion.NeutralTheoreticalProduct.ProductType == ProductType.y
            ? sequenceLength - ion.NeutralTheoreticalProduct.FragmentNumber
            : ion.NeutralTheoreticalProduct.AminoAcidPosition;
        x = GetX(residueNum);
        y = GetY(row) + MetaDrawSettings.ProductTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType];

        // is internal
        if (!ion.IsInternalFragment)
        {
            if (ion.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime)
            {
                DrawnSequence.DrawCTermIon(SequenceDrawingCanvas, new Point(x, y), color, "", 2);
            }
            else if (ion.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime)
            {
                DrawnSequence.DrawNTermIon(SequenceDrawingCanvas, new Point(x, y), color, "", 2);
            }
        }
    }

    private void AddCircles(ChimericSpectralMatchModel psm, int row, int maxBaseSeqLength)
    {
        var color = DrawnSequence.ParseColorBrushFromOxyColor(psm.Color);
        DrawnSequence.DrawCircle(SequenceDrawingCanvas, new Point(GetX(maxBaseSeqLength + 1), GetY(row)), color);
        DrawnSequence.DrawCircle(SequenceDrawingCanvas, new Point(GetX(-2), GetY(row)), color);
    }

    private void DrawBaseSequence(ChimericSpectralMatchModel psm, int row)
    {
        var baseSeq = psm.Psm.BaseSeq.Split('|')[0];
        int index = 0;
        for (; index < baseSeq.Length; index++)
        {
            var x = GetX(index);
            var y = GetY(row);
            DrawnSequence.DrawText(SequenceDrawingCanvas, new Point(x, y), baseSeq[index].ToString(), Brushes.Black);
        }
    }

    private void AddModifications(ChimericSpectralMatchModel psm, int row)
    {
        foreach (var mod in psm.AllModsOneIsNterminus)
        {
            var x = GetX(mod.Key - 2);
            var y = GetY(row);
            var color = DrawnSequence.ParseColorBrushFromOxyColor(psm.Color);
            DrawnSequence.DrawCircle(SequenceDrawingCanvas, new Point(x, y), color);
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
}
