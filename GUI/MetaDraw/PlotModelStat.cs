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

            }
            else if (plotType.Equals("Fragment PPM Error vs. RT"))
            {

            }
            else if (plotType.Equals("Histograms of Count of Different PTMs Seen at 1% FDR"))
            {
                histogramPlot(5);
            }
        }

        public void histogramPlot(int plotType)
        {
            int binRange = allPSM.Count/6;
            double dataCategory= 0;
            SortedList<double,double> numCategory = new SortedList<double,double>();

            int[] values = new int[binRange];
            List<double> numbers = new List<double>();
            List<string> axes = new List<string>();
            var s1 = new ColumnSeries { ColumnWidth = 200, IsStacked = false };
            foreach (var psm in allPSM)
            {
                switch (plotType)
                {
                    case 1:
                        dataCategory = psm.PrecursorMass;
                        break;
                    case 2:
                        //dataCategory = psm.PrecursorCharge;
                        break;
                    case 3:
                        dataCategory = psm.PrecursorCharge;
                        break;
                    case 4:
                        //dataCategory = psm.PrecursorCharge.ToString();
                        break;
                    case 5:
                        //dataCategory = psm.;
                        break;
                }
                numbers.Add(dataCategory);
                
            }
            double maxValue = numbers.Max();
            int decimalPlaces = 1;
            if (maxValue >= 100)
            {
                decimalPlaces = 0;
            }
            double valRange = maxValue / binRange;
            foreach( var a in numbers)
            {
                if (a == maxValue)
                {
                    values[binRange-1]++;
                }
                else
                {
                    values[(int)(a / valRange)]++;
                }
            }
            
            for(int i = 0; i < binRange; i++)
            {
                s1.Items.Add(new ColumnItem(values[i]));
                var leftLimit = Math.Round(valRange * i, decimalPlaces);
                var rightLimit = Math.Round(valRange * (i + 1), decimalPlaces);
                axes.Add( leftLimit + " - " + rightLimit);
            }

            privateModel.Series.Add(s1);
            privateModel.Axes.Add(new CategoryAxis
            {
                Position = AxisPosition.Bottom,
                ItemsSource = axes
            });
        }
    }
}
