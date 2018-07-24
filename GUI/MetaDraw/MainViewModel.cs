using System.Linq;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;
using EngineLayer.CrosslinkSearch;
using EngineLayer;
using System.ComponentModel;
using MetaMorpheusGUI;
using System.IO;
using MassSpectrometry;
using System.Collections.Generic;

namespace ViewModels
{
    public class MainViewModel : INotifyPropertyChanged
    {
        private const double STROKE_THICKNESS_UNANNOTATED = 0.5;
        private const double STROKE_THICKNESS_ANNOTATED = 2.0;
        private PlotModel privateModel;

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        { { ProductType.B, OxyColors.Blue },
          { ProductType.BnoB1ions, OxyColors.Blue },
          { ProductType.Y, OxyColors.Purple },
          { ProductType.C, OxyColors.Gold },
          { ProductType.Zdot, OxyColors.Orange } };

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

        public MainViewModel()
        {
            // Create and Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            privateModel = new PlotModel { Title = "Spectrum Annotation", Subtitle = "using OxyPlot" };
        }

        // single peptides (not crosslink)
        public void DrawPeptideSpectralMatch(MsDataScan msDataScan, MetaDrawPsm psmToDraw, DrawParams drawParams)
        {
            // x is m/z, y is intensity
            var spectrumMzs = msDataScan.MassSpectrum.XArray;
            var spectrumIntensities = msDataScan.MassSpectrum.YArray;

            PlotModel model = new PlotModel { Title = "Spectrum Annotation of Scan #" + msDataScan.OneBasedScanNumber, DefaultFontSize = 15, Subtitle = psmToDraw.FullSequence };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = spectrumMzs.Max() * 1.02, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumMzs.Max() * 5 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity", Minimum = 0, Maximum = spectrumIntensities.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumIntensities.Max() * 1.3 });
            model.Axes[1].Zoom(0, spectrumIntensities.Max() * 1.1);

            LineSeries[] allIons = new LineSeries[spectrumMzs.Length];

            // draw the matched peaks; if the PSM is null, we're just drawing the peaks in the scan without annotation, so skip this part
            if (psmToDraw != null)
            {
                foreach (var peak in psmToDraw.FragmentIons)
                {
                    OxyColor ionColor = productTypeDrawColors[peak.ProductType];

                    int i = msDataScan.MassSpectrum.GetClosestPeakIndex(peak.Mz).Value;

                    // peak line
                    allIons[i] = new LineSeries();
                    allIons[i].Color = ionColor;
                    allIons[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
                    allIons[i].Points.Add(new DataPoint(peak.Mz, 0));
                    allIons[i].Points.Add(new DataPoint(peak.Mz, spectrumIntensities[i]));

                    // peak annotation
                    var peakAnnotation = new TextAnnotation();
                    peakAnnotation.TextRotation = -60;
                    peakAnnotation.Font = "Arial";
                    peakAnnotation.FontSize = 12;
                    peakAnnotation.FontWeight = 2.0;
                    peakAnnotation.TextColor = ionColor;
                    peakAnnotation.StrokeThickness = 0;
                    peakAnnotation.Text = "(" + peak.Mz.ToString("F3") + ") " + peak.ProductType.ToString().ToLower() + "-" + peak.IonNumber;
                    peakAnnotation.TextPosition = new DataPoint(allIons[i].Points[1].X, allIons[i].Points[1].Y + peakAnnotation.Text.Length * 1.5 / 4);
                    peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Left;
                    model.Annotations.Add(peakAnnotation);

                    model.Series.Add(allIons[i]);
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

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.Model = model;
        }

        // TODO: Crosslink annotation
        public void DrawPsmCross(Ms2ScanWithSpecificMass MsScanForDraw, PsmCross psmParentsForDraw)
        {
            var x = MsScanForDraw.TheScan.MassSpectrum.XArray;
            var y = MsScanForDraw.TheScan.MassSpectrum.YArray;

            string scanNum = psmParentsForDraw.ScanNumber.ToString();
            string sequence1 = psmParentsForDraw.FullSequence + "-" + psmParentsForDraw.XlPos.ToString();
            string sequence2 = psmParentsForDraw.BetaPsmCross.FullSequence + "-" + psmParentsForDraw.BetaPsmCross.XlPos.ToString();

            var matchedIonDic1 = psmParentsForDraw.MatchedIonInfo;
            var matchedIonDic2 = psmParentsForDraw.BetaPsmCross.MatchedIonInfo;

            PlotModel model = new PlotModel { Title = "Spectrum anotation of Scan " + scanNum + " for Crosslinked Peptide", DefaultFontSize = 15 };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = x.Max() * 1.02, AbsoluteMinimum = 0 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity(counts)", Minimum = 0, Maximum = y.Max() * 1.2, AbsoluteMinimum = 0 });
            var textAnnoSeq1 = new TextAnnotation() { };
            //textAnnoSeq1.TextRotation=90;
            textAnnoSeq1.FontSize = 12; textAnnoSeq1.TextColor = OxyColors.Red; textAnnoSeq1.StrokeThickness = 0; textAnnoSeq1.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.15); textAnnoSeq1.Text = sequence1;
            var textAnnoSeq2 = new TextAnnotation() { };
            //textAnnoSeq2.TextRotation=90;
            textAnnoSeq2.FontSize = 12; textAnnoSeq2.TextColor = OxyColors.Blue; textAnnoSeq2.StrokeThickness = 0; textAnnoSeq2.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.1); textAnnoSeq2.Text = sequence2;
            model.Annotations.Add(textAnnoSeq1);
            model.Annotations.Add(textAnnoSeq2);

            LineSeries[] s0 = new LineSeries[x.Length];
            LineSeries[] s1 = new LineSeries[x.Length];
            LineSeries[] s2 = new LineSeries[x.Length];

            //Draw the ms/ms scan peaks
            for (int i = 0; i < x.Length; i++)
            {
                s0[i] = new LineSeries();
                s0[i].Color = OxyColors.DimGray;
                s0[i].StrokeThickness = STROKE_THICKNESS_UNANNOTATED;
                s0[i].Points.Add(new DataPoint(x[i], 0));
                s0[i].Points.Add(new DataPoint(x[i], y[i]));
                model.Series.Add(s0[i]);
            }
            //Draw the ms/ms scan matched peaks

            for (int i = 0; i < matchedIonDic1.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Red;
                if (matchedIonDic1.MatchedIonMz[i] > 0)
                {
                    s1[i] = new LineSeries();
                    s1[i].Color = ionColor;
                    s1[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, 0));
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, matchedIonDic1.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    //textAnno1.TextRotation=90;
                    textAnno1.FontSize = 12;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s1[i].Points[1];
                    textAnno1.Text = matchedIonDic1.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    //textAnno2.TextRotation=90;
                    textAnno2.FontSize = 12;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s1[i].Points[1].X, s1[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic1.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s1[i]);
                }

            }

            for (int i = 0; i < matchedIonDic2.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Blue;
                if (matchedIonDic2.MatchedIonMz[i] > 0)
                {
                    s2[i] = new LineSeries();
                    s2[i].Color = ionColor;
                    s2[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i] + 1.007277, 0));
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i] + 1.007277, matchedIonDic2.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    //textAnno1.TextRotation=90;
                    textAnno1.FontSize = 12;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s2[i].Points[1];
                    textAnno1.Text = matchedIonDic2.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    //textAnno2.TextRotation=90;
                    textAnno2.FontSize = 12;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s2[i].Points[1].X, s2[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic2.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s2[i]);
                }
            }

            // Axes are created automatically if they are not defined

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.Model = model;
        }

        // TODO: Draw to pdf file
        public void XLDrawMSMatchToPdf(Ms2ScanWithSpecificMass MsScanForDraw, PsmCross psmParentsForDraw, int order, string OutputFolder)
        {
            var x = MsScanForDraw.TheScan.MassSpectrum.XArray;
            var y = MsScanForDraw.TheScan.MassSpectrum.YArray;

            string scanNum = psmParentsForDraw.ScanNumber.ToString();
            string sequence1 = psmParentsForDraw.FullSequence + "-" + psmParentsForDraw.XlPos.ToString();
            string sequence2 = psmParentsForDraw.BetaPsmCross.FullSequence + "-" + psmParentsForDraw.BetaPsmCross.XlPos.ToString();

            var matchedIonDic1 = psmParentsForDraw.MatchedIonInfo;
            var matchedIonDic2 = psmParentsForDraw.BetaPsmCross.MatchedIonInfo;

            PlotModel model = new PlotModel { Title = "Spectrum anotation of Scan " + scanNum + " for Crosslinked Peptide", DefaultFontSize = 15 };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = x.Max() * 1.02, AbsoluteMinimum = 0 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity(counts)", Minimum = 0, Maximum = y.Max() * 1.2, AbsoluteMinimum = 0 });
            var textAnnoSeq1 = new TextAnnotation() { };
            //textAnnoSeq1.TextRotation=90;
            textAnnoSeq1.FontSize = 12; textAnnoSeq1.TextColor = OxyColors.Red; textAnnoSeq1.StrokeThickness = 0; textAnnoSeq1.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.15); textAnnoSeq1.Text = sequence1;
            var textAnnoSeq2 = new TextAnnotation() { };
            //textAnnoSeq2.TextRotation=90;
            textAnnoSeq2.FontSize = 12; textAnnoSeq2.TextColor = OxyColors.Blue; textAnnoSeq2.StrokeThickness = 0; textAnnoSeq2.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.1); textAnnoSeq2.Text = sequence2;
            model.Annotations.Add(textAnnoSeq1);
            model.Annotations.Add(textAnnoSeq2);

            LineSeries[] s0 = new LineSeries[x.Length];
            LineSeries[] s1 = new LineSeries[x.Length];
            LineSeries[] s2 = new LineSeries[x.Length];

            //Draw the ms/ms scan peaks
            for (int i = 0; i < x.Length; i++)
            {
                s0[i] = new LineSeries();
                s0[i].Color = OxyColors.DimGray;
                s0[i].StrokeThickness = STROKE_THICKNESS_UNANNOTATED;
                s0[i].Points.Add(new DataPoint(x[i], 0));
                s0[i].Points.Add(new DataPoint(x[i], y[i]));
                model.Series.Add(s0[i]);
            }
            //Draw the ms/ms scan matched peaks

            for (int i = 0; i < matchedIonDic1.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Red;
                if (matchedIonDic1.MatchedIonMz[i] > 0)
                {
                    s1[i] = new LineSeries();
                    s1[i].Color = ionColor;
                    s1[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, 0));
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, matchedIonDic1.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    //textAnno1.TextRotation=90;
                    textAnno1.FontSize = 12;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s1[i].Points[1];
                    textAnno1.Text = matchedIonDic1.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    //textAnno2.TextRotation=90;
                    textAnno2.FontSize = 12;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s1[i].Points[1].X, s1[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic1.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s1[i]);
                }

            }

            for (int i = 0; i < matchedIonDic2.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Blue;
                if (matchedIonDic2.MatchedIonMz[i] > 0)
                {
                    s2[i] = new LineSeries();
                    s2[i].Color = ionColor;
                    s2[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i] + 1.007277, 0));
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i] + 1.007277, matchedIonDic2.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    //textAnno1.TextRotation=90;
                    textAnno1.FontSize = 12;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s2[i].Points[1];
                    textAnno1.Text = matchedIonDic2.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    //textAnno2.TextRotation=90;
                    textAnno2.FontSize = 12;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s2[i].Points[1].X, s2[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic2.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s2[i]);
                }
            }

            using (var stream = File.Create(OutputFolder + "\\" + order.ToString() + "_Scan" + scanNum + ".pdf"))
            {
                PdfExporter pdf = new PdfExporter { Width = 500, Height = 210 };
                pdf.Export(model, stream);
            }
        }
    }
}
