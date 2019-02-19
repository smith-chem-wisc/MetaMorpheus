using MassSpectrometry;
using MetaMorpheusGUI;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;

namespace ViewModels
{
    public class PlotModelStat : INotifyPropertyChanged
    {
        private const double STROKE_THICKNESS_UNANNOTATED = 0.5;
        private const double STROKE_THICKNESS_ANNOTATED = 2.0;
        private PlotModel privateModel;
        private List<MetaDrawPsm> allPSM;
        public List<string> plotNames = new List<string>{ "Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)",
                                                          "Histogram of Fragment PPM Errors",
                                                          "Histogram of Precursor Charges",
                                                          "Histogram of Fragment Charges",
                                                          "Precursor PPM Error vs. RT",
                                                          "Fragment PPM Error vs. RT",
                                                          "Histograms of Count of Different PTMs Seen at 1% FDR"};

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

        public PlotModelStat()
        {
            
        }

        public PlotModelStat(string plotName, List<MetaDrawPsm> psms)
        {
            privateModel = new PlotModel { Title = plotName, Subtitle = "using OxyPlot" };
            allPSM = psms;
            createPlot(plotName);
        }

        public List<string> plotTypes()
        {
            return plotNames;
        }

        public void createPlot(string plotType)
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
            else if (plotType.Equals("Histograms of Count of Different PTMs Seen at 1% FDR"))
            {
                histogramPlot(5);
            }
        }

        public void histogramPlot(int plotType)
        {
            double binSize = 0;
            SortedList<double,double> numCategory = new SortedList<double,double>();
            IEnumerable<double> numbers = new List<double>();
            List<string> axes = new List<string>();
            var s1 = new ColumnSeries { ColumnWidth = 200, IsStacked = false };
            switch (plotType)
            {
                case 1:
                    numbers = allPSM.Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa), 0) == 0).Select(p => double.Parse(p.MassDiffPpm));
                    binSize = 0.1;
                    break;
                case 2:
                    numbers = allPSM.SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm));
                    binSize = 0.1;
                    break;
                case 3:
                    numbers = allPSM.Select(p => (double)(p.PrecursorCharge));
                    binSize = 1;
                    break;
                case 4:
                    numbers = allPSM.SelectMany(p => p.MatchedIons.Select(v => (double)v.Charge));
                    binSize = 1;
                    break;
                case 5:
                    //dataCategory = 
                    binSize = 1;
                    break;
            }
            int[,] values = new int[numbers.Count(),2];
            double maxValue = numbers.Max();
            int decimalPlaces = 0;
            int sign = 0;
            if (binSize.Equals(0.1))
            {
                decimalPlaces = 1;
            }
            foreach(var a in numbers)
            {
                if (a == maxValue)
                {
                   values[numbers.Count()-1,0]++;
                }
                else
                {
                    var current = a;
                    if (a < 0)
                    {
                        current = a * -1;
                        sign = 1;
                    }
                    values[(int) Math.Round(current, decimalPlaces),sign]++;
                }
            }
            
            foreach(var n in values)
            {
                s1.Items.Add(new ColumnItem(values[n,0]));
                s1.Items.Add(new ColumnItem(values[n,1]));
                //var leftLimit = Math.Round(valRange * i, decimalPlaces);
                //var rightLimit = Math.Round(valRange * (i + 1), decimalPlaces);
                //axes.Add( leftLimit + " - " + rightLimit);
                axes.Add(n.ToString());
            }

            privateModel.Series.Add(s1);
            privateModel.Axes.Add(new CategoryAxis
            {
                Position = AxisPosition.Bottom,
                ItemsSource = axes
            });
        }

        public void linePlot(int plotType)
        {
            ScatterSeries series = new ScatterSeries();
            List<Tuple<double,double>> xy = new List<Tuple<double, double>>();
            var filteredList = allPSM.Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa), 0) == 0).ToList();
            switch (plotType)
            {
                case 1:
                    foreach (var psm in filteredList)
                    {
                        xy.Add(new Tuple<double, double>(double.Parse(psm.MassDiffPpm), psm.RetentionTime));
                    }
                    break;
                case 2:
                    foreach (var psm in filteredList)
                    {
                        xy.Add(new Tuple<double, double>(double.Parse(psm.MassDiffPpm), psm.RetentionTime));
                    }
                    break;
            }
            IOrderedEnumerable<Tuple<double, double>> sorted = xy.OrderBy(x => x.Item1);
            foreach(var val in sorted)
            {
                series.Points.Add(new ScatterPoint(val.Item1, val.Item2));
            }
            privateModel.Series.Add(series);
        }
    }
}
