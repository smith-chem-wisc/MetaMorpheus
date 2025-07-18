using System.Linq;
using System.Windows;
using System;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace GuiFunctions.MetaDraw;
public enum LegendDisplayProperty
{
    ProteinName,
    ProteinAccession,
    BaseSequence,
    FullSequence,
    Modifications
}

public class ChimeraLegendCanvas : Canvas
{
    public ChimeraGroupViewModel GroupViewModel { get; }

    public ChimeraLegendCanvas(ChimeraGroupViewModel groupViewModel)
    {
        GroupViewModel = groupViewModel;
        BuildLegend();
    }

    private void BuildLegend()
    {
        Children.Clear();
        double y = 10;
        double x = 10;
        double rowHeight = 24;
        double maxTextWidth = 0;

        // Group by legend key (refactor to allow dynamic grouping if needed)
        var proteinGroups = GroupViewModel.ChimericPsms
            .GroupBy(psm => psm.Psm.Name)
            .OrderByDescending(g => g.Count());

        int proteinIndex = 0;
        foreach (var group in proteinGroups)
        {
            var proteoforms = group.ToList();
            var proteinColor = proteoforms[0].ProteinColor;

            if (proteoforms.Count == 1)
            {
                maxTextWidth = Math.Max(maxTextWidth, MeasureTextWidth(GetMainText(proteoforms[0]), 12, FontWeights.Regular));
                var ellipse = new Ellipse
                {
                    Width = 12,
                    Height = 12,
                    Fill = new SolidColorBrush(DrawnSequence.ParseColorFromOxyColor(proteoforms[0].Color)),
                    Stroke = Brushes.Black,
                    StrokeThickness = 1
                };
                SetLeft(ellipse, x);
                SetTop(ellipse, y);

                var text = new TextBlock
                {
                    Text = GetMainText(proteoforms[0]),
                    FontWeight = FontWeights.DemiBold,
                    FontSize = 12,
                    Margin = new Thickness(6, 0, 0, 0)
                };
                SetLeft(text, x + 18);
                SetTop(text, y - 2);

                Children.Add(ellipse);
                Children.Add(text);

                y += rowHeight;
            }
            else
            {
                maxTextWidth = Math.Max(maxTextWidth, MeasureTextWidth("Shared Ions", 12, FontWeights.Regular));
                var sharedEllipse = new Ellipse
                {
                    Width = 12,
                    Height = 12,
                    Fill = new SolidColorBrush(DrawnSequence.ParseColorFromOxyColor(proteinColor)),
                    Stroke = Brushes.Black,
                    StrokeThickness = 1
                };
                SetLeft(sharedEllipse, x);
                SetTop(sharedEllipse, y);

                var header = new TextBlock
                {
                    Text = GetMainText(proteoforms[0]),
                    FontWeight = FontWeights.DemiBold,
                    FontSize = 12
                };
                SetLeft(header, x + 18);
                SetTop(header, y - 2);

                Children.Add(sharedEllipse);
                Children.Add(header);
                y += rowHeight;

                for (int i = 0; i < proteoforms.Count; i++)
                {
                    var color = proteoforms[i].Color;
                    var ellipse = new Ellipse
                    {
                        Width = 12,
                        Height = 12,
                        Fill = new SolidColorBrush(DrawnSequence.ParseColorFromOxyColor(color)),
                        Stroke = Brushes.Black,
                        StrokeThickness = 1
                    };
                    SetLeft(ellipse, x + 10);
                    SetTop(ellipse, y);

                    var text = new TextBlock
                    {
                        Text = GetSubText(proteoforms[i]),
                        FontWeight = FontWeights.Regular,
                        FontSize = 12,
                        Margin = new Thickness(6, 0, 0, 0)
                    };
                    SetLeft(text, x + 28);
                    SetTop(text, y - 2);

                    Children.Add(ellipse);
                    Children.Add(text);

                    y += rowHeight;
                    maxTextWidth = Math.Max(maxTextWidth, MeasureTextWidth(GetSubText(proteoforms[i]), 12, FontWeights.Regular));
                }
            }
            proteinIndex++;
        }

        double leftMargin = 10;
        double ellipseAndSpacing = 18;
        double rightMargin = 10;
        Width = leftMargin + ellipseAndSpacing + maxTextWidth + rightMargin;
        Height = y + 10;
    }

    private string GetMainText(ChimericSpectralMatchModel psm)
    {
        switch (MetaDrawSettings.ChimeraLegendMainTextType)
        {
            case LegendDisplayProperty.ProteinAccession: return SanitizeIfAmbiguous(psm.Psm.Accession);
            case LegendDisplayProperty.BaseSequence: return SanitizeIfAmbiguous(psm.Psm.BaseSeq);
            case LegendDisplayProperty.FullSequence: return SanitizeIfAmbiguous(psm.Psm.FullSequence ?? psm.Psm.BaseSeq);
            case LegendDisplayProperty.Modifications: return string.IsNullOrEmpty(psm.ModString) ? "Unmodified" : psm.ModString;
            case LegendDisplayProperty.ProteinName: 
            default: return SanitizeIfAmbiguous(psm.Psm.Name);
        }
    }

    private string GetSubText(ChimericSpectralMatchModel psm)
    {
        switch (MetaDrawSettings.ChimeraLegendSubTextType)
        {
            case LegendDisplayProperty.ProteinAccession: return SanitizeIfAmbiguous(psm.Psm.Accession);
            case LegendDisplayProperty.BaseSequence: return SanitizeIfAmbiguous(psm.Psm.BaseSeq);
            case LegendDisplayProperty.FullSequence: return SanitizeIfAmbiguous(psm.Psm.FullSequence ?? psm.Psm.BaseSeq);
            case LegendDisplayProperty.ProteinName: return SanitizeIfAmbiguous(psm.Psm.Name);
            case LegendDisplayProperty.Modifications: 
            default: return string.IsNullOrEmpty(psm.ModString) ? "Unmodified" : psm.ModString;
        }
    }

    private double MeasureTextWidth(string text, double fontSize = 12, FontWeight? fontWeight = null)
    {
        var typeface = new Typeface("Segoe UI");
        var formattedText = new FormattedText(
            text ?? "",
            System.Globalization.CultureInfo.CurrentCulture,
            FlowDirection.LeftToRight,
            typeface,
            fontSize,
            Brushes.Black,
            new NumberSubstitution(),
            1.0);

        if (fontWeight.HasValue)
            formattedText.SetFontWeight(fontWeight.Value);

        return formattedText.WidthIncludingTrailingWhitespace;
    }

    private string SanitizeIfAmbiguous(string input)
    {
        if (MetaDrawSettings.ChimeraLegendTakeFirstIfAmbiguous)
        {
            var parts = input.Split('|');
            if (parts.Length > 1)
            {
                return parts[0]; // Return the first part if ambiguous
            }
        }
        return input; // Return the original input if not ambiguous or if setting is false
    }
}