using System.Linq;
using System.Windows;
using System;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace GuiFunctions;
public enum LegendDisplayProperty
{
    ProteinName,
    ProteinAccession,
    BaseSequence,
    FullSequence
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
        this.Children.Clear();
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
            var proteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][0];

            if (proteoforms.Count == 1)
            {
                maxTextWidth = Math.Max(maxTextWidth, MeasureTextWidth(group.Key, 12, FontWeights.Regular));
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
                    Text = group.Key,
                    FontWeight = FontWeights.DemiBold,
                    FontSize = 12,
                    Margin = new Thickness(6, 0, 0, 0)
                };
                SetLeft(text, x + 18);
                SetTop(text, y - 2);

                this.Children.Add(ellipse);
                this.Children.Add(text);

                y += rowHeight;
            }
            else
            {
                maxTextWidth = Math.Max(maxTextWidth, MeasureTextWidth("Shared Ions", 12, FontWeights.Regular));

                var header = new TextBlock
                {
                    Text = group.Key,
                    FontWeight = FontWeights.DemiBold,
                    FontSize = 12
                };
                SetLeft(header, x);
                SetTop(header, y);
                this.Children.Add(header);
                y += rowHeight;

                var sharedEllipse = new Ellipse
                {
                    Width = 12,
                    Height = 12,
                    Fill = new SolidColorBrush(DrawnSequence.ParseColorFromOxyColor(proteinColor)),
                    Stroke = Brushes.Black,
                    StrokeThickness = 1
                };
                SetLeft(sharedEllipse, x + 10);
                SetTop(sharedEllipse, y);

                var sharedText = new TextBlock
                {
                    Text = "Shared Ions",
                    FontWeight = FontWeights.DemiBold,
                    FontSize = 12,
                    Margin = new Thickness(6, 0, 0, 0)
                };
                SetLeft(sharedText, x + 28);
                SetTop(sharedText, y - 2);

                this.Children.Add(sharedEllipse);
                this.Children.Add(sharedText);

                y += rowHeight;

                for (int i = 0; i < proteoforms.Count; i++)
                {
                    var color = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][i + 1];
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
                        Text = string.IsNullOrEmpty(proteoforms[i].ModString) ? "" : proteoforms[i].ModString,
                        FontWeight = FontWeights.Regular,
                        FontSize = 12,
                        Margin = new Thickness(6, 0, 0, 0)
                    };
                    SetLeft(text, x + 28);
                    SetTop(text, y - 2);

                    this.Children.Add(ellipse);
                    this.Children.Add(text);

                    y += rowHeight;
                    if (!string.IsNullOrEmpty(proteoforms[i].ModString))
                        maxTextWidth = Math.Max(maxTextWidth, MeasureTextWidth(proteoforms[i].ModString, 12, FontWeights.Regular));
                }
            }
            proteinIndex++;
        }

        double leftMargin = 10;
        double ellipseAndSpacing = 28;
        double rightMargin = 10;
        this.Width = leftMargin + ellipseAndSpacing + maxTextWidth + rightMargin;
        this.Height = y + 10;
    }

    private string GetLegendKey(ChimericSpectralMatchModel psm)
    {
        switch (MetaDrawSettings.LegendDisplayMode)
        {
            case LegendDisplayProperty.ProteinAccession:
                return psm.Psm.Accession;
            case LegendDisplayProperty.BaseSequence:
                return psm.Psm.BaseSeq;
            case LegendDisplayProperty.FullSequence:
                return psm.Psm.FullSequence ?? psm.Psm.BaseSeq;
            case LegendDisplayProperty.ProteinName:
            default:
                return psm.Psm.Name;
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
}