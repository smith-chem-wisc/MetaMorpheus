using EngineLayer;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;

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
            "Histogram of PTM Spectral Counts"
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
                    category[counter] = d.Key;
                    counter++;
                }
                this.privateModel.Axes.Add(new CategoryAxis
                {
                    ItemsSource = category
                });
                privateModel.Series.Add(column);
            }
            else
            {
                double end = numbers.Max();
                double start = numbers.Min();
                double bins = (end - start) / binSize;
                double numbins = bins * Math.Pow(10, normalizeNumber(bins));
                int exp = (int)Math.Pow(10, normalizeNumber(end));

                long size = Convert.ToInt64(end * exp + (-1 * start * exp) + 1);
                int[,] values = new int[size, 2];

                int decimalPlaces = 0;

                if (binSize.Equals(0.1))
                {
                    decimalPlaces = 1;
                }
                foreach (var a in numbers)
                {
                    int sign = 0;
                    var current = a * exp;

                    if (a < 0)
                    {
                        current = a * -1;
                        sign = 1;
                    }
                    values[(int)Math.Round(current, decimalPlaces), sign]++;
                }

                int zeroIndex = values.Length / 2;

                //add negative value bars
                for (int i = (values.Length / 2); i > 0; i--)
                {
                    s1.Items.Add(new ColumnItem(values[i - 1, 1], zeroIndex - i));
                    axes.Add((-i).ToString());
                }

                s1.Items.Add(new ColumnItem(values[0, 0] + values[0, 1], zeroIndex));
                axes.Add(0.ToString());
                //add positive value bars
                for (int i = 1; i < values.Length / 2; i++)
                {
                    s1.Items.Add(new ColumnItem(values[i, 0], zeroIndex + i));
                    axes.Add(i.ToString());
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
            List<Tuple<double, double>> xy = new List<Tuple<double, double>>();
            var filteredList = allPsms.Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa), 0) == 0).ToList();
            var test = allPsms.SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm));
            switch (plotType)
            {
                case 1:
                    foreach (var psm in filteredList)
                    {
                        xy.Add(new Tuple<double, double>(double.Parse(psm.MassDiffPpm), (double)psm.RetentionTime));
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
            }
            IOrderedEnumerable<Tuple<double, double>> sorted = xy.OrderBy(x => x.Item1);
            foreach (var val in sorted)
            {
                series.Points.Add(new ScatterPoint(val.Item2, val.Item1));
            }
            series.MarkerFill = OxyColors.Blue;
            series.MarkerSize = 0.5;
            privateModel.Series.Add(series);
        }

        private static int normalizeNumber(double number)
        {
            string s = number.ToString("00.00E0");
            int i = Convert.ToInt32(s.Substring(s.Length - 1)) / 10;
            return i;
        }

        //unused interface methods
        public void Update(bool updateData) { }
        public void Render(IRenderContext rc, double width, double height) { }
        public void AttachPlotView(IPlotView plotView) { }
    }
}

