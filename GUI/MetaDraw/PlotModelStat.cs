using EngineLayer;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Globalization;

namespace MetaMorpheusGUI
{
    public class PlotModelStat : INotifyPropertyChanged, IPlotModel
    {
        private PlotModel privateModel;
        private ObservableCollection<PsmFromTsv> allPsms;

        public static List<string> PlotNames = new List<string> {
            "Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)",
            //"Histogram of Fragment PPM Errors", //TODO: implement fragment PPM error reading in MetaDraw
            "Histogram of Precursor Charges",
            "Histogram of Fragment Charges",
            "Precursor PPM Error vs. RT",
            //"Fragment PPM Error vs. RT", //TODO: implement fragment PPM error reading in MetaDraw
            "Histogram of PTM Spectral Counts",
            "Predicted RT vs. Observed RT"
        };

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        {
            { ProductType.b, OxyColors.Blue },
            { ProductType.y, OxyColors.Purple },
            { ProductType.c, OxyColors.Gold },
            { ProductType.zPlusOne, OxyColors.Orange },
            { ProductType.D, OxyColors.DodgerBlue },
            { ProductType.M, OxyColors.Firebrick }
        };

        public PlotModel Model
        {
            get
            {
                return privateModel;
            }
            private set
            {
                privateModel = value;
                NotifyPropertyChanged("Model");
            }
        }

        public OxyColor Background => OxyColors.White;

        public event PropertyChangedEventHandler PropertyChanged;

        protected void NotifyPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(propertyName));
            }
        }

        public PlotModelStat(string plotName, ObservableCollection<PsmFromTsv> psms)
        {
            privateModel = new PlotModel { Title = plotName };
            allPsms = psms;
            createPlot(plotName);
        }

        private void createPlot(string plotType)
        {
            if (plotType.Equals("Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)"))
            {
                histogramPlot(1);
            }
            else if (plotType.Equals("Histogram of Fragment PPM Errors"))
            {
                histogramPlot(2);
            }
            else if (plotType.Equals("Histogram of Precursor Charges"))
            {
                histogramPlot(3);
            }
            else if (plotType.Equals("Histogram of Fragment Charges"))
            {
                histogramPlot(4);
            }
            else if (plotType.Equals("Precursor PPM Error vs. RT"))
            {
                linePlot(1);
            }
            else if (plotType.Equals("Fragment PPM Error vs. RT"))
            {
                linePlot(2);
            }
            else if (plotType.Equals("Histogram of PTM Spectral Counts"))
            {
                histogramPlot(5);
            }
            else if (plotType.Equals("Predicted RT vs. Observed RT"))
            {
                linePlot(3);
            }
        }

        private void histogramPlot(int plotType)
        {
            double binSize = -1;
            SortedList<double, double> numCategory = new SortedList<double, double>();
            IEnumerable<double> numbers = new List<double>();
            Dictionary<string, int> dict = new Dictionary<string, int>();

            List<string> axes = new List<string>();
            var s1 = new ColumnSeries { ColumnWidth = 200, IsStacked = false, FillColor = OxyColors.Blue };

            switch (plotType)
            {
                case 1:
                    numbers = allPsms.Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa), 0) == 0).Select(p => double.Parse(p.MassDiffPpm));
                    binSize = 0.1;
                    break;
                case 2:
                    numbers = allPsms.SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm));
                    binSize = 0.1;
                    break;
                case 3:
                    numbers = allPsms.Select(p => (double)(p.PrecursorCharge));
                    var results = numbers.GroupBy(p => p).OrderBy(p => p.Key).Select(p => p);
                    dict = results.ToDictionary(p => p.Key.ToString(), v => v.Count());
                    break;
                case 4:
                    numbers = allPsms.SelectMany(p => p.MatchedIons.Select(v => (double)v.Charge));
                    results = numbers.GroupBy(p => p).OrderBy(p => p.Key).Select(p => p);
                    dict = results.ToDictionary(p => p.Key.ToString(), v => v.Count());
                    break;
                case 5:
                    var psmsWithMods = allPsms.Where(p => !p.FullSequence.Contains("|") && p.FullSequence.Contains("["));
                    var mods = psmsWithMods.Select(p => new PeptideWithSetModifications(p.FullSequence, GlobalVariables.AllModsKnownDictionary)).Select(p => p.AllModsOneIsNterminus).SelectMany(p => p.Values);
                    var groupedMods = mods.GroupBy(p => p.IdWithMotif).ToList();
                    dict = groupedMods.ToDictionary(p => p.Key, v => v.Count());
                    break;
            }
            if (plotType >= 3)
            {
                ColumnSeries column = new ColumnSeries { ColumnWidth = 200, IsStacked = false, FillColor = OxyColors.Blue };
                var counter = 0;
                String[] category = new string[dict.Count];
                foreach (var d in dict)
                {
                    column.Items.Add(new ColumnItem(d.Value, counter));
                    int spaceIndex = d.Key.IndexOf(' ');
                    category[counter] = (plotType != 5 || spaceIndex < 0) ? d.Key : d.Key.Substring(0, spaceIndex) + "\n" + d.Key.Substring(spaceIndex);  // replace the first space with a newline for the PTM plot
                    counter++;
                }
                double labelAngle = plotType == 5 ? -50 : 0;
                this.privateModel.Axes.Add(new CategoryAxis
                {
                    ItemsSource = category, Angle = labelAngle
                });
                privateModel.Series.Add(column);
            }
            else
            {
                double end = roundToBin(numbers.Max(), binSize);
                double start = roundToBin(numbers.Min(), binSize);
                int numBins = (int)(((end - start) / binSize) + 1.001); // + 0.001 ensures double is just above its int value before truncating

                int[] values = new int[numBins];
                
                // put each value into the nearest bin
                foreach (var a in numbers)
                {
                    values[(int)(((roundToBin(a, binSize) - start) / binSize) + 0.001)]++;    // + 0.001 ensures double is just above its int value before truncating
                }

                int minBinLabels = 28;  // the number of labeled bins will be between minBinLabels and 2 * minBinLabels
                int skipBinLabel = numBins < minBinLabels ? 1 : numBins / minBinLabels;

                // create a column and axis label for each bin
                for (int i = 0; i < values.Length; i++)
                {
                    s1.Items.Add(new ColumnItem(values[i], i));
                    if(i % skipBinLabel == 0)
                    {
                        axes.Add(roundToBin(start + (i * binSize), binSize).ToString(CultureInfo.InvariantCulture));  // numbers need to be re-rounded so values like 0.20000000001 aren't displayed
                    }
                    else
                    {
                        axes.Add("");
                    }
                }

                privateModel.Series.Add(s1);
                privateModel.Axes.Add(new CategoryAxis
                {
                    Position = AxisPosition.Bottom,
                    ItemsSource = axes
                });
            }
        }

        private void linePlot(int plotType)
        {
            ScatterSeries series = new ScatterSeries();
            ScatterSeries variantSeries = new ScatterSeries();  // used by plot 3 for variant contianing peptides
            List<Tuple<double, double>> xy = new List<Tuple<double, double>>();
            List<Tuple<double, double>> variantxy = new List<Tuple<double, double>>();  // used by plot 3 for variant containing peptides
            var filteredList = allPsms.Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa), 0) == 0).ToList();
            var test = allPsms.SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm));
            switch (plotType)
            {
                case 1:
                    foreach (var psm in filteredList)
                    {
                        if(psm.IdentifiedSequenceVariations == null || psm.IdentifiedSequenceVariations.Equals(""))
                        {
                            xy.Add(new Tuple<double, double>(double.Parse(psm.MassDiffPpm), (double)psm.RetentionTime));
                        }
                        else
                        {
                            variantxy.Add(new Tuple<double, double>(double.Parse(psm.MassDiffPpm), (double)psm.RetentionTime));
                        }
                    }
                    break;
                case 2:
                    foreach (var psm in allPsms)
                    {
                        foreach (var ion in psm.MatchedIons)
                        {
                            xy.Add(new Tuple<double, double>((double)psm.RetentionTime, ion.MassErrorPpm));
                        }
                    }
                    break;
                case 3:
                    SSRCalc3 sSRCalc3 = new SSRCalc3("A100", SSRCalc3.Column.A100);
                    foreach (var psm in allPsms)
                    {
                        if(psm.IdentifiedSequenceVariations == null || psm.IdentifiedSequenceVariations.Equals(""))
                        {
                            xy.Add(new Tuple<double, double>(sSRCalc3.ScoreSequence(new PeptideWithSetModifications(psm.BaseSeq.Split('|')[0], null)), 
                            (double)psm.RetentionTime));
                        }
                        else
                        {
                            variantxy.Add(new Tuple<double, double>(sSRCalc3.ScoreSequence(new PeptideWithSetModifications(psm.BaseSeq.Split('|')[0], null)),
                            (double)psm.RetentionTime));
                        }
                    }
                    break;
            }
            IOrderedEnumerable<Tuple<double, double>> sorted = xy.OrderBy(x => x.Item1);
            foreach (var val in sorted)
            {
                series.Points.Add(new ScatterPoint(val.Item2, val.Item1));
            }
            series.MarkerFill = OxyColors.Blue;
            series.MarkerSize = 0.5;
            privateModel.Series.Add(series);

            // plot the variant containing peptides
            if (variantxy.Count != 0)
            {
                IOrderedEnumerable<Tuple<double, double>> variantSorted = variantxy.OrderBy(x => x.Item1);
                foreach (var val in variantSorted)
                {
                    variantSeries.Points.Add(new ScatterPoint(val.Item2, val.Item1));
                }
                variantSeries.MarkerFill = OxyColors.DarkRed;
                variantSeries.MarkerSize = 1.5;
                privateModel.Series.Add(variantSeries);
            }
        }

        // rounds a number to the nearest multiple of binsize, midpoints are rounded towards zero
        private static double roundToBin(double number, double binSize)
        {
            int sign = number < 0 ? -1 : 1;
            double d = number * sign;
            double remainder = d % binSize;
            d = remainder < 0.5 * binSize ? d - remainder : d - remainder + binSize;
            return d * sign;
        }

        //unused interface methods
        public void Update(bool updateData) { }
        public void Render(IRenderContext rc, double width, double height) { }
        public void AttachPlotView(IPlotView plotView) { }
    }
}

