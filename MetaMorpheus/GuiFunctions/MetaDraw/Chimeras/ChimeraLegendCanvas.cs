using System.Linq;
using System.Windows;
using System;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;
using System.Collections.Generic;

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
    private const int TextFontSize = 12;
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

        var proteinGroups = GroupViewModel.ChimericPsms
            .GroupBy(psm => psm.Psm.Name)
            .OrderByDescending(g => g.Count());

        int proteinIndex = 0;
        foreach (var group in proteinGroups)
        {
            // Avoid repeated ToList() calls and repeated property access
            var proteoforms = group.ToList();
            int proteoformCount = proteoforms.Count;
            var firstProteoform = proteoforms[0];
            var proteinColor = firstProteoform.ProteinColor;

            if (proteoformCount == 1)
            {
                // Cache main text and its width
                string mainText = GetMainText(firstProteoform);
                double mainTextWidth = MeasureTextWidth(mainText, TextFontSize, FontWeights.Regular);
                maxTextWidth = Math.Max(maxTextWidth, mainTextWidth);

                var ellipse = new Ellipse
                {
                    Width = TextFontSize,
                    Height = TextFontSize,
                    Fill = new SolidColorBrush(DrawnSequence.ParseColorFromOxyColor(firstProteoform.Color)),
                    Stroke = Brushes.Black,
                    StrokeThickness = 1
                };
                SetLeft(ellipse, x);
                SetTop(ellipse, y);

                var text = new TextBlock
                {
                    Text = mainText,
                    FontWeight = FontWeights.DemiBold,
                    FontSize = TextFontSize,
                    Margin = new Thickness(6, 0, 0, 0),
                    TextWrapping = TextWrapping.Wrap,
                    MaxWidth = MetaDrawSettings.ChimeraLegendMaxWidth
                };

                text.Measure(new Size(MetaDrawSettings.ChimeraLegendMaxWidth, double.PositiveInfinity));
                double textHeight = Math.Max(rowHeight, text.DesiredSize.Height);

                SetLeft(text, x + 18);
                SetTop(text, y - 2);

                Children.Add(ellipse);
                Children.Add(text);

                y += textHeight;
            }
            else
            {
                // Only measure "Shared Ions" once
                const string sharedIonsText = "Shared Ions";
                double sharedIonsWidth = MeasureTextWidth(sharedIonsText, TextFontSize, FontWeights.Regular);
                maxTextWidth = Math.Max(maxTextWidth, sharedIonsWidth);

                var sharedEllipse = new Ellipse
                {
                    Width = TextFontSize,
                    Height = TextFontSize,
                    Fill = new SolidColorBrush(DrawnSequence.ParseColorFromOxyColor(proteinColor)),
                    Stroke = Brushes.Black,
                    StrokeThickness = 1
                };
                SetLeft(sharedEllipse, x);
                SetTop(sharedEllipse, y);

                string headerText = GetMainText(firstProteoform);
                var header = new TextBlock
                {
                    Text = headerText,
                    FontWeight = FontWeights.DemiBold,
                    FontSize = TextFontSize,
                    TextWrapping = TextWrapping.Wrap,
                    MaxWidth = MetaDrawSettings.ChimeraLegendMaxWidth
                };

                header.Measure(new Size(MetaDrawSettings.ChimeraLegendMaxWidth, double.PositiveInfinity));
                double textHeight = Math.Max(rowHeight, header.DesiredSize.Height);

                SetLeft(header, x + 18);
                SetTop(header, y - 2);
                Children.Add(sharedEllipse);
                Children.Add(header);

                y += textHeight;

                // Cache subtexts and their widths to avoid repeated calls
                for (int i = 0; i < proteoformCount; i++)
                {
                    var pf = proteoforms[i];
                    var color = pf.Color;
                    string subText = GetSubText(pf);
                    double subTextWidth = MeasureTextWidth(subText, TextFontSize, FontWeights.Regular);

                    var ellipse = new Ellipse
                    {
                        Width = TextFontSize,
                        Height = TextFontSize,
                        Fill = new SolidColorBrush(DrawnSequence.ParseColorFromOxyColor(color)),
                        Stroke = Brushes.Black,
                        StrokeThickness = 1
                    };
                    SetLeft(ellipse, x + 10);
                    SetTop(ellipse, y);

                    var text = new TextBlock
                    {
                        Text = subText,
                        FontWeight = FontWeights.Regular,
                        FontSize = TextFontSize,
                        Margin = new Thickness(6, 0, 0, 0),
                        TextWrapping = TextWrapping.Wrap,
                        MaxWidth = MetaDrawSettings.ChimeraLegendMaxWidth
                    };

                    text.Measure(new Size(MetaDrawSettings.ChimeraLegendMaxWidth, double.PositiveInfinity));
                    textHeight = Math.Max(rowHeight, text.DesiredSize.Height);

                    SetLeft(text, x + 28);
                    SetTop(text, y - 2);

                    Children.Add(ellipse);
                    Children.Add(text);

                    y += textHeight;
                    maxTextWidth = Math.Max(maxTextWidth, subTextWidth);
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
            case LegendDisplayProperty.FullSequence: return SanitizeIfAmbiguous(psm.Psm.FullSequence);
            case LegendDisplayProperty.Modifications: return string.IsNullOrEmpty(psm.ModString) ? "Unmodified" : SanitizeIfAmbiguous(psm.ModString);
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
            case LegendDisplayProperty.FullSequence: return SanitizeIfAmbiguous(psm.Psm.FullSequence);
            case LegendDisplayProperty.ProteinName: return SanitizeIfAmbiguous(psm.Psm.Name);
            case LegendDisplayProperty.Modifications: 
            default: return string.IsNullOrEmpty(psm.ModString) ? "Unmodified" : SanitizeIfAmbiguous(psm.ModString);
        }
    }

    private double MeasureTextWidth(string text, double fontSize = TextFontSize, FontWeight? fontWeight = null)
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