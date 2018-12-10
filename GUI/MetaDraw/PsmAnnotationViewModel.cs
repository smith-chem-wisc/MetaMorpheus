using Chemistry;
using EngineLayer;
using MassSpectrometry;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;

namespace ViewModels
{
    public class PsmAnnotationViewModel : INotifyPropertyChanged
    {
        private const double STROKE_THICKNESS_UNANNOTATED = 0.5;
        private const double STROKE_THICKNESS_ANNOTATED = 2.0;
        private PlotModel privateModel;

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        {
          { ProductType.b, OxyColors.Blue },
          { ProductType.y, OxyColors.Purple },
          { ProductType.c, OxyColors.Gold },
          { ProductType.zDot, OxyColors.Orange },
          { ProductType.D, OxyColors.DodgerBlue },
          { ProductType.M, OxyColors.Firebrick }
        };

        private static Dictionary<ProductType, OxyColor> betaPeptideProductTypeDrawColors = new Dictionary<ProductType, OxyColor>
        {
          { ProductType.b, OxyColors.LightBlue },
          { ProductType.y, OxyColors.MediumPurple },
          { ProductType.c, OxyColors.LightGoldenrodYellow },
          { ProductType.zDot, OxyColors.OrangeRed },
          { ProductType.D, OxyColors.AliceBlue },
          { ProductType.M, OxyColors.LightCoral }
        };

        public PlotModel Model
        {
            get
            {
                return this.privateModel;
            }
            set
            {
                this.privateModel = value;
                NotifyPropertyChanged("Model");
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected void NotifyPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(propertyName));
            }
        }

        public PsmAnnotationViewModel()
        {
            // Create and Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            privateModel = new PlotModel { Title = "Spectrum Annotation", Subtitle = "using OxyPlot" };
        }

        // single peptides (not crosslink)
        public void DrawPeptideSpectralMatch(MsDataScan msDataScan, PsmFromTsv psmToDraw)
        {
            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.Model = Draw(msDataScan, psmToDraw);
        }

        private PlotModel Draw(MsDataScan msDataScan, PsmFromTsv psmToDraw)
        {
            // x is m/z, y is intensity
            var spectrumMzs = msDataScan.MassSpectrum.XArray;
            var spectrumIntensities = msDataScan.MassSpectrum.YArray;

            string subtitle = psmToDraw.FullSequence;
            if (psmToDraw != null && psmToDraw.BetaPeptideBaseSequence != null)
            {
                subtitle = psmToDraw.FullSequence + "\n" + psmToDraw.BetaPeptideFullSequence;
            }
            PlotModel model = new PlotModel { Title = "Spectrum Annotation of Scan #" + msDataScan.OneBasedScanNumber, DefaultFontSize = 15, Subtitle = subtitle };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = spectrumMzs.Max() * 1.02, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumMzs.Max() * 5 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity", Minimum = 0, Maximum = spectrumIntensities.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumIntensities.Max() * 1.3 });
            model.Axes[1].Zoom(0, spectrumIntensities.Max() * 1.1);

            LineSeries[] allIons = new LineSeries[spectrumMzs.Length];

            // draw the matched peaks; if the PSM is null, we're just drawing the peaks in the scan without annotation, so skip this part
            if (psmToDraw != null)
            {
                foreach (var peak in psmToDraw.MatchedIons)
                {
                    OxyColor ionColor;

                    if (productTypeDrawColors.ContainsKey(peak.NeutralTheoreticalProduct.ProductType))
                    {
                        ionColor = productTypeDrawColors[peak.NeutralTheoreticalProduct.ProductType];
                    }
                    else
                    {
                        ionColor = OxyColors.Turquoise;
                    }

                    int i = msDataScan.MassSpectrum.GetClosestPeakIndex(peak.NeutralTheoreticalProduct.NeutralMass.ToMz(peak.Charge)).Value;

                    // peak line
                    allIons[i] = new LineSeries();
                    allIons[i].Color = ionColor;
                    allIons[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
                    allIons[i].Points.Add(new DataPoint(peak.Mz, 0));
                    allIons[i].Points.Add(new DataPoint(peak.Mz, spectrumIntensities[i]));

                    // peak annotation
                    string peakAnnotationText = peak.NeutralTheoreticalProduct.ProductType.ToString().ToLower() + peak.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + " (" + peak.Mz.ToString("F3") + ")";
                    if (peak.NeutralTheoreticalProduct.NeutralLoss != 0)
                    {
                        peakAnnotationText = peak.NeutralTheoreticalProduct.ProductType.ToString().ToLower() + peak.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + "-" + peak.NeutralTheoreticalProduct.NeutralLoss.ToString("F2") + " (" + peak.Mz.ToString("F3") + ")";
                    }

                    var peakAnnotation = new TextAnnotation();
                    peakAnnotation.TextRotation = -60;
                    peakAnnotation.Font = "Arial";
                    peakAnnotation.FontSize = 12;
                    peakAnnotation.FontWeight = 2.0;
                    peakAnnotation.TextColor = ionColor;
                    peakAnnotation.StrokeThickness = 0;
                    peakAnnotation.Text = peakAnnotationText;
                    peakAnnotation.TextPosition = new DataPoint(allIons[i].Points[1].X, allIons[i].Points[1].Y + peakAnnotation.Text.Length * 1.5 / 4);
                    peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Left;
                    model.Annotations.Add(peakAnnotation);

                    model.Series.Add(allIons[i]);
                }

                if (psmToDraw.BetaPeptideBaseSequence != null)
                {
                    foreach (var peak in psmToDraw.BetaPeptideMatchedIons)
                    {
                        OxyColor ionColor;

                        if (productTypeDrawColors.ContainsKey(peak.NeutralTheoreticalProduct.ProductType))
                        {
                            ionColor = betaPeptideProductTypeDrawColors[peak.NeutralTheoreticalProduct.ProductType];
                        }
                        else
                        {
                            ionColor = OxyColors.Turquoise;
                        }

                        int i = msDataScan.MassSpectrum.GetClosestPeakIndex(peak.NeutralTheoreticalProduct.NeutralMass.ToMz(peak.Charge)).Value;

                        // peak line
                        allIons[i] = new LineSeries();
                        allIons[i].Color = ionColor;
                        allIons[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
                        allIons[i].Points.Add(new DataPoint(peak.Mz, 0));
                        allIons[i].Points.Add(new DataPoint(peak.Mz, spectrumIntensities[i]));

                        // peak annotation
                        string peakAnnotationText = "beta-" + peak.NeutralTheoreticalProduct.ProductType.ToString().ToLower() + peak.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + " (" + peak.Mz.ToString("F3") + ")";
                        if (peak.NeutralTheoreticalProduct.NeutralLoss != 0)
                        {
                            peakAnnotationText = "beta-" + peak.NeutralTheoreticalProduct.ProductType.ToString().ToLower() + peak.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + "-" + peak.NeutralTheoreticalProduct.NeutralLoss.ToString("F2") + " (" + peak.Mz.ToString("F3") + ")";
                        }

                        var peakAnnotation = new TextAnnotation();
                        peakAnnotation.TextRotation = -60;
                        peakAnnotation.Font = "Arial";
                        peakAnnotation.FontSize = 12;
                        peakAnnotation.FontWeight = 2.0;
                        peakAnnotation.TextColor = ionColor;
                        peakAnnotation.StrokeThickness = 0;
                        peakAnnotation.Text = peakAnnotationText;
                        peakAnnotation.TextPosition = new DataPoint(allIons[i].Points[1].X, allIons[i].Points[1].Y + peakAnnotation.Text.Length * 1.5 / 4);
                        peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Left;
                        model.Annotations.Add(peakAnnotation);

                        model.Series.Add(allIons[i]);
                    }
                }
            }

            // draw the remaining unmatched peaks
            for (int i = 0; i < spectrumMzs.Length; i++)
            {
                // peak has already been drawn (it is a matched peak)
                if (allIons[i] != null)
                {
                    continue;
                }

                allIons[i] = new LineSeries();
                allIons[i].Color = OxyColors.DimGray;
                allIons[i].StrokeThickness = STROKE_THICKNESS_UNANNOTATED;
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], 0));
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], spectrumIntensities[i]));
                model.Series.Add(allIons[i]);
            }

            // Axes are created automatically if they are not defined      
            return model;
        }

        private PlotModel DrawPdf(MsDataScan msDataScan, PropertyInfo[] properties, PsmFromTsv psm, bool redraw)
        {
            if (redraw)
            {
                this.Model = Draw(msDataScan, psm);
            }

            var y = Model.DefaultYAxis.ActualMaximum - Model.DefaultYAxis.ActualMaximum * 0.03;
            var x = Model.DefaultXAxis.ActualMaximum - Model.DefaultXAxis.ActualMaximum * 0.01;
            var diff = (y - (Model.DefaultYAxis.ActualMaximum * 0.1)) / properties.Length;

            // properties to exclude
            string[] exclude = { "PrecursorScanNum", "PrecursorMz", "MatchedIons",
                "TotalIonCurrent", "Notch", " EssentialSeq", "MissedCleavage",
                "OrganismName", "StartAndEndResiduesInProtein", "PreviousAminoAcid",
                "NextAminoAcid", "Rank", "BetaPeptideMatchedIons", "BetaPeptideRank" };

            var propertiesList = properties.ToList();
            foreach (PropertyInfo property in properties.ToList())
            {
                if (exclude.Any(e => e.Contains(property.Name)))
                {
                    propertiesList.Remove(property);
                }
            }

            var displayedProperties = propertiesList.Where(p => p.GetValue(psm) != null); // only display non-null properties

            foreach (PropertyInfo property in displayedProperties)
            {
                var propertyAnnotation = new TextAnnotation
                {
                    Text = property.Name + ": " + property.GetValue(psm),
                    TextPosition = new DataPoint(x, y),
                    FontSize = 9,
                    StrokeThickness = 0,
                    TextHorizontalAlignment = HorizontalAlignment.Right
                };

                y -= diff;
                Model.Annotations.Add(propertyAnnotation);
            }

            return Model;
        }

        public void DrawPeptideSpectralMatchPdf(MsDataScan msDataScan, PsmFromTsv psm, string fileName, bool redraw)
        {
            var properties = psm.GetType().GetProperties();
            var pdfModel = DrawPdf(msDataScan, properties, psm, redraw);

            string dir = Path.GetDirectoryName(fileName);
            if (!Directory.Exists(dir))
            {
                Directory.CreateDirectory(dir);
            }

            using (var stream = File.Create(fileName))
            {
                PdfExporter pdf = new PdfExporter { Width = 800, Height = 500 };
                pdf.Export(pdfModel, stream);
            }
        }
    }
}
