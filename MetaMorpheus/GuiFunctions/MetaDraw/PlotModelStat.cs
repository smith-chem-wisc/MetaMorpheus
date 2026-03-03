using EngineLayer;
using GuiFunctions.MetaDraw;
using Omics.Fragmentation;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;

namespace GuiFunctions
{
    public class PlotModelStat
    {
        private PlotModel privateModel;
        private readonly ObservableCollection<SpectrumMatchFromTsv> allSpectralMatches;
        private readonly Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> psmsBySourceFile;
        private readonly PlotModelStatParameters parameters;

        public readonly static List<string> PlotNames = new List<string> {
            "Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)",
            "Histogram of Fragment PPM Errors",
            "Histogram of Precursor Charges",
            "Histogram of Fragment Charges",
            "Histogram of Precursor Masses",
            "Histogram of Precursor m/z",
            "Histogram of Hydrophobicity scores",
            "Precursor PPM Error vs. RT",
            "Histogram of PTM Spectral Counts",
            "Predicted RT vs. Observed RT",
            "Histogram of Missed Cleavages",
            "Histogram of Fragment Ion Types by Count",
            "Histogram of Fragment Ion Types by Intensity",
            "Histogram of Ids by Retention Time"
        };

        public PlotModel Model => privateModel;

        public PlotModelStat(string plotName, ObservableCollection<SpectrumMatchFromTsv> sms, 
            Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> smsBySourceFile, 
            PlotModelStatParameters parameters = null)
        {
            privateModel = new PlotModel { Title = plotName, DefaultFontSize = 14 };
            allSpectralMatches = sms;
            this.psmsBySourceFile = smsBySourceFile;
            this.parameters = parameters ?? new PlotModelStatParameters();
            createPlot(plotName);
            privateModel.DefaultColors = MetaDrawSettings.DataVisualizationColorOrder;
        }
        private void createPlot(string plotType)
        {
            switch (plotType)
            {
                case "Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)":
                    histogramPlot(1);
                    break;
                case "Histogram of Fragment PPM Errors": 
                    histogramPlot(2);
                    break;
                case "Histogram of Precursor Charges":
                    histogramPlot(3);
                    break;
                case "Histogram of Fragment Charges":
                    histogramPlot(4);
                    break;
                case "Histogram of Precursor Masses":
                    histogramPlot(6);
                    break;
                case "Histogram of Precursor m/z":
                    histogramPlot(7);
                    break;
                case "Histogram of Hydrophobicity scores":
                    histogramPlot(8);
                    break;
                case "Precursor PPM Error vs. RT":
                    linePlot(1);
                    break;
                case "Histogram of PTM Spectral Counts":
                    histogramPlot(5);
                    break;
                case "Predicted RT vs. Observed RT":
                    linePlot(3);
                    break;
                case "Histogram of Missed Cleavages":
                    histogramPlot(9);
                    break;
                case "Histogram of Fragment Ion Types by Count":
                    histogramPlot(10);
                    break;
                case "Histogram of Fragment Ion Types by Intensity":
                    histogramPlot(11);
                    break;
                case "Histogram of Ids by Retention Time":
                    histogramPlot(12);
                    break;
            }
        }
        private void histogramPlot(int plotType)
        {
            privateModel.LegendTitle = "Source file(s)";
            string yAxisTitle = "Count";
            string xAxisTitle = "";
            double binSize = -1;
            double labelAngle = 0;
            SortedList<double, double> numCategory = new SortedList<double, double>();
            Dictionary<string, IEnumerable<double>> numbersBySourceFile = new Dictionary<string, IEnumerable<double>>();
            Dictionary<string, Dictionary<string, int>> dictsBySourceFile = new Dictionary<string, Dictionary<string, int>>();
            
            bool isGroupingEnabled = parameters.GroupingProperty != "None";
            Dictionary<string, Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>>> groupedPsmsBySourceFile = null;

            if (isGroupingEnabled)
            {
                groupedPsmsBySourceFile = GroupPsmsByProperty(psmsBySourceFile, parameters.GroupingProperty);
            }

            switch (plotType)
            {
                case 1: // Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)
                    xAxisTitle = "Precursor error (ppm)";
                    binSize = 0.1;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        numbersBySourceFile.Add(key, psmsBySourceFile[key].Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa, CultureInfo.InvariantCulture), 0) == 0).Select(p => double.Parse(p.MassDiffPpm, CultureInfo.InvariantCulture)));
                        var results = numbersBySourceFile[key].GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 2: // Histogram of Fragment PPM Errors 
                    xAxisTitle = "Fragment error (ppm)";
                    binSize = 1;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        numbersBySourceFile.Add(key, psmsBySourceFile[key].SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm)));
                        var results = numbersBySourceFile[key].GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 3: // Histogram of Precursor Charges
                    xAxisTitle = "Precursor charge";
                    binSize = 1;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        numbersBySourceFile.Add(key, psmsBySourceFile[key].Select(p => (double)(p.PrecursorCharge)));
                        var results = numbersBySourceFile[key].GroupBy(p => p).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 4: // Histogram of Fragment Charges
                    xAxisTitle = "Fragment charge";
                    binSize = 1;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        numbersBySourceFile.Add(key, psmsBySourceFile[key].SelectMany(p => p.MatchedIons.Select(v => (double)v.Charge)));
                        var results = numbersBySourceFile[key].GroupBy(p => p).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 5: // Histogram of PTM Spectral Counts
                    xAxisTitle = "Modification";
                    labelAngle = -50;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        var psmsWithMods = psmsBySourceFile[key].Where(p => !p.FullSequence.Contains("|") && p.FullSequence.Contains("["));
                        var mods = psmsWithMods.Select(p => p.ToBioPolymerWithSetMods()).Select(p => p.AllModsOneIsNterminus).SelectMany(p => p.Values);
                        var groupedMods = mods.GroupBy(p => p.IdWithMotif).ToList();
                        dictsBySourceFile.Add(key, groupedMods.ToDictionary(p => p.Key, v => v.Count()));
                    }
                    break;
                case 6: // Histogram of Precursor mass
                    xAxisTitle = "Precursor Mass (Da)";
                    binSize = 100;
                    labelAngle = -50;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        numbersBySourceFile.Add(key, psmsBySourceFile[key].Select(p => (double)(p.PrecursorMass)));
                        var results = numbersBySourceFile[key].GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 7: // Histogram of Precursor m/z
                    xAxisTitle = "Precursor mass/charge";
                    binSize = 50;
                    labelAngle = -50;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        numbersBySourceFile.Add(key, psmsBySourceFile[key].Select(p => (double)(p.PrecursorMz)));
                        var results = numbersBySourceFile[key].GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 8: // Histogram of Hydrophobicity Scores
                    xAxisTitle = "Hydrophobicity Score (Determined by SSRCalc)";
                    binSize = 2;
                    labelAngle = -50;
                    SSRCalc3 sSRCalc3 = new SSRCalc3("A100", SSRCalc3.Column.A100);
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        var values = new List<double>();
                        foreach (var psm in psmsBySourceFile[key].Where(p => p is not OsmFromTsv))
                        {
                            values.Add(sSRCalc3.ScoreSequence(new PeptideWithSetModifications(psm.BaseSeq.Split("|")[0], null)));
                           
                        }
                        numbersBySourceFile.Add(key, values);
                        var results = numbersBySourceFile[key].GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;

                case 9: // Histogram of Missed Cleavages
                    xAxisTitle = "Missed Cleavages";
                    binSize = 1;
                    labelAngle = 0;
                    foreach (var fileName in psmsBySourceFile.Keys)
                    {
                        var values = psmsBySourceFile[fileName].Where(p => !p.MissedCleavage.Contains("|")).Select(p => double.Parse(p.MissedCleavage)).ToList();
                        numbersBySourceFile.Add(fileName, values);
                        var results = numbersBySourceFile[fileName].GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(fileName, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 10: // Histogram of Fragment Ion Types by count
                    xAxisTitle = "Fragment Types";
                    labelAngle = 0;

                    foreach (var fileName in psmsBySourceFile.Keys)
                    {
                        var allMatchedIons = psmsBySourceFile[fileName].SelectMany(p => p.MatchedIons).ToList();
                        
                        var result = allMatchedIons
                            .GroupBy(p => p.NeutralTheoreticalProduct is CustomMProduct cmp ? cmp.Annotation : p.NeutralTheoreticalProduct.ProductType.ToString())
                            .ToDictionary(p => p.Key.ToString(), p => p.Count());
                        dictsBySourceFile.Add(fileName, result);
                    }
                    break;
                case 11: // Histogram of Fragment Ion Types by intensity
                    xAxisTitle = "Fragment Types";
                    labelAngle = 0;
                    yAxisTitle = "Summed Intensity";
                    foreach (var fileName in psmsBySourceFile.Keys)
                    {
                        var allMatchedIons = psmsBySourceFile[fileName].SelectMany(p => p.MatchedIons).ToList();
                        
                        var result = allMatchedIons
                            .GroupBy(p => p.NeutralTheoreticalProduct is CustomMProduct cmp ? cmp.Annotation : p.NeutralTheoreticalProduct.ProductType.ToString())
                            .ToDictionary(p => p.Key.ToString(), p => (int)p.Sum(m => m.Intensity));
                        dictsBySourceFile.Add(fileName, result);
                    }
                    break;
                case 12: // Histogram of Ids by Retention Time
                    xAxisTitle = "Retention Time";
                    binSize = 1;
                    labelAngle = 0;

                    foreach (var fileName in psmsBySourceFile.Keys)
                    {
                        var result = psmsBySourceFile[fileName]
                            .GroupBy(p => (int)Math.Round(p.RetentionTime, 0))
                            .ToDictionary(p => p.Key.ToString(CultureInfo.InvariantCulture), p => p.Count());
                        dictsBySourceFile.Add(fileName, result);
                    }
                    break;
            }

            String[] category;
            int[] totalCounts;
            
            if (plotType == 5 || plotType == 10 || plotType == 11)  // category histogram
            {
                IEnumerable<string> allCategories = dictsBySourceFile.Values.Select(p => p.Keys).SelectMany(p => p);
                Dictionary<string, int> categoryIDs = new Dictionary<string, int>();
                int counter = 0;
                foreach (string s in allCategories)
                {
                    if (!categoryIDs.ContainsKey(s))
                    {
                        categoryIDs.Add(s, counter++);
                    }
                }
                category = new string[counter];
                totalCounts = new int[counter];

                foreach (string cat in categoryIDs.Keys)
                {
                    foreach (Dictionary<string, int> dict in dictsBySourceFile.Values)
                    {
                        totalCounts[categoryIDs[cat]] += dict.ContainsKey(cat) ? dict[cat] : 0;
                    }
                }

                int maxCount = totalCounts.Max();
                double minThreshold = maxCount * (parameters.MinRelativeCutoff / 100.0);
                double maxThreshold = maxCount * (parameters.MaxRelativeCutoff / 100.0);

                var filteredCategoryIDs = new Dictionary<string, int>();
                var filteredCategory = new List<string>();
                var filteredTotalCounts = new List<int>();
                int newIndex = 0;

                foreach (var cat in categoryIDs.Keys.OrderBy(k => k))
                {
                    int oldId = categoryIDs[cat];
                    if (totalCounts[oldId] >= minThreshold && totalCounts[oldId] <= maxThreshold)
                    {
                        filteredCategoryIDs[cat] = newIndex;
                        filteredCategory.Add(cat);
                        filteredTotalCounts.Add(totalCounts[oldId]);
                        newIndex++;
                    }
                }

                category = filteredCategory.ToArray();
                totalCounts = filteredTotalCounts.ToArray();
                categoryIDs = filteredCategoryIDs;
            
                if (isGroupingEnabled)
                {
                    var allGroupKeys = groupedPsmsBySourceFile.Values
                        .SelectMany(d => d.Keys)
                        .Distinct()
                        .OrderBy(k => k)
                        .ToList();

                    var nestedCategories = new List<string>();
                    var categoryToIndex = new Dictionary<string, int>();
                    int categoryIndex = 0;

                    foreach (var groupKey in allGroupKeys)
                    {
                        foreach (var cat in filteredCategoryIDs.Keys.OrderBy(k => k))
                        {
                            nestedCategories.Add(cat);
                            categoryToIndex[cat + "|" + groupKey] = categoryIndex++;
                        }
                    }

                    foreach (string sourceFile in psmsBySourceFile.Keys)
                    {
                        ColumnSeries column = new ColumnSeries
                        {
                            ColumnWidth = 200,
                            IsStacked = true,
                            Title = sourceFile,
                            TrackerFormatString = "Category: {bin}\n{0}: {2}\nGroup: {group}\nTotal: {total}",
                            BaseValue = parameters.UseLogScaleYAxis ? 0.1 : 0
                        };

                        foreach (string groupKey in allGroupKeys)
                        {
                            if (!groupedPsmsBySourceFile[sourceFile].ContainsKey(groupKey))
                            {
                                foreach (var cat in filteredCategoryIDs.Keys.OrderBy(k => k))
                                {
                                    string lookupKey = cat + "|" + groupKey;
                                    if (categoryToIndex.ContainsKey(lookupKey))
                                    {
                                        column.Items.Add(new HistItem(0, categoryToIndex[lookupKey], cat, totalCounts[filteredCategoryIDs[cat]], groupKey));
                                    }
                                }
                                continue;
                            }

                            var groupPsms = groupedPsmsBySourceFile[sourceFile][groupKey];
                            var allMatchedIons = groupPsms.SelectMany(p => p.MatchedIons).ToList();

                            Dictionary<string, int> groupDict;
                            if (plotType == 5)
                            {
                                var psmsWithMods = groupPsms.Where(p => !p.FullSequence.Contains("|") && p.FullSequence.Contains("["));
                                var mods = psmsWithMods.Select(p => p.ToBioPolymerWithSetMods()).Select(p => p.AllModsOneIsNterminus).SelectMany(p => p.Values);
                                var groupedMods = mods.GroupBy(p => p.IdWithMotif).ToList();
                                groupDict = groupedMods.ToDictionary(p => p.Key, v => v.Count());
                            }
                            else if (plotType == 10)
                            {
                                groupDict = allMatchedIons
                                    .GroupBy(p => p.NeutralTheoreticalProduct is CustomMProduct cmp ? cmp.Annotation : p.NeutralTheoreticalProduct.ProductType.ToString())
                                    .ToDictionary(p => p.Key.ToString(), p => p.Count());
                            }
                            else
                            {
                                groupDict = allMatchedIons
                                    .GroupBy(p => p.NeutralTheoreticalProduct is CustomMProduct cmp ? cmp.Annotation : p.NeutralTheoreticalProduct.ProductType.ToString())
                                    .ToDictionary(p => p.Key.ToString(), p => (int)p.Sum(m => m.Intensity));
                            }

                            foreach (var cat in filteredCategoryIDs.Keys.OrderBy(k => k))
                            {
                                string lookupKey = cat + "|" + groupKey;
                                if (!categoryToIndex.ContainsKey(lookupKey))
                                    continue;

                                double total = 1.0;
                                if (parameters.NormalizeHistogramToFile)
                                {
                                    total = groupDict.Values.Sum();
                                }
                                
                                int catId = filteredCategoryIDs[cat];
                                double value = groupDict.ContainsKey(cat) ? groupDict[cat] / total : 0;
                                
                                if (parameters.UseLogScaleYAxis && value > 0 && value < 0.1)
                                    value = 0.1;
                                    
                                column.Items.Add(new HistItem(value, categoryToIndex[lookupKey], cat, totalCounts[catId], groupKey));
                            }
                        }
                        privateModel.Series.Add(column);
                    }

                    category = nestedCategories.ToArray();
                    
                    // Create group axis labels with proper spacing for category histograms
                    var groupAxisLabels = new List<string>();
                    for (int i = 0; i < allGroupKeys.Count; i++)
                    {
                        for (int j = 0; j < filteredCategoryIDs.Count; j++)
                        {
                            if (j == filteredCategoryIDs.Count / 2)
                            {
                                groupAxisLabels.Add(allGroupKeys[i]);
                            }
                            else
                            {
                                groupAxisLabels.Add("");
                            }
                        }
                    }
                    
                    var groupAxis = new CategoryAxis
                    {
                        Position = AxisPosition.Bottom,
                        Title = parameters.GroupingProperty,
                        Key = "GroupAxis",
                        StartPosition = 1,
                        EndPosition = 0,
                        IsAxisVisible = true,
                        ItemsSource = groupAxisLabels.ToArray(),
                        FontSize = 12,
                        AxisTickToLabelDistance = 5,
                        MajorTickSize = 10
                    };
                    privateModel.Axes.Add(groupAxis);
                }
                else
                {
                    foreach (string key in dictsBySourceFile.Keys)
                    {
                        ColumnSeries column = new ColumnSeries { 
                            ColumnWidth = 200, 
                            IsStacked = true, 
                            Title = key, 
                            TrackerFormatString = "Bin: {bin}\n{0}: {2}\nTotal: {total}",
                            BaseValue = parameters.UseLogScaleYAxis ? 0.1 : 0
                        };
                        
                        foreach (var d in dictsBySourceFile[key])
                        {
                            if (!filteredCategoryIDs.ContainsKey(d.Key))
                                continue;

                            int id = filteredCategoryIDs[d.Key];
                            double total = 1.0;
                            if (parameters.NormalizeHistogramToFile)
                            {
                                total = dictsBySourceFile[key].Values.Sum();
                            }
                            
                            double value = d.Value / total;
                            
                            if (parameters.UseLogScaleYAxis && value < 0.1)
                                value = 0.1;
                                
                            column.Items.Add(new HistItem(value, id, d.Key, totalCounts[id]));
                        }
                        privateModel.Series.Add(column);
                    }
                }
            }
            else    // numerical histogram
            {
                IEnumerable<double> allNumbers = numbersBySourceFile.Values.SelectMany(x => x);

                int end = dictsBySourceFile.Values.Max(p => p.Max(v => int.Parse(v.Key)));
                int start = dictsBySourceFile.Values.Min(p => p.Min(v => int.Parse(v.Key)));
                int numBins = end - start + 1;
                int minBinLabels = 22;
                int skipBinLabel = numBins < minBinLabels ? 1 : numBins / minBinLabels;

                category = new string[numBins];
                totalCounts = new int[numBins];
                
                int maxCount = 0;
                foreach (Dictionary<string, int> dict in dictsBySourceFile.Values)
                {
                    foreach (var kvp in dict)
                    {
                        int idx = int.Parse(kvp.Key) - start;
                        if (idx >= 0 && idx < totalCounts.Length)
                            totalCounts[idx] += kvp.Value;
                        if (totalCounts[idx] > maxCount)
                            maxCount = totalCounts[idx];
                    }
                }

                double minThreshold = maxCount * (parameters.MinRelativeCutoff / 100.0);
                double maxThreshold = maxCount * (parameters.MaxRelativeCutoff / 100.0);

                var filteredCategories = new List<string>();
                var filteredTotalCounts = new List<int>();
                var binToFilteredIndex = new Dictionary<int, int>();
                int filteredIndex = 0;

                for (int i = start; i <= end; i++)
                {
                    int binIdx = i - start;
                    if (totalCounts[binIdx] >= minThreshold && totalCounts[binIdx] <= maxThreshold)
                    {
                        binToFilteredIndex[i] = filteredIndex;
                        if (i % skipBinLabel == 0)
                        {
                            filteredCategories.Add(Math.Round((i * binSize), 2).ToString(CultureInfo.InvariantCulture));
                        }
                        else
                        {
                            filteredCategories.Add("");
                        }
                        filteredTotalCounts.Add(totalCounts[binIdx]);
                        filteredIndex++;
                    }
                }

                category = filteredCategories.ToArray();
                totalCounts = filteredTotalCounts.ToArray();

                if (isGroupingEnabled)
                {
                    var allGroupKeys = groupedPsmsBySourceFile.Values
                        .SelectMany(d => d.Keys)
                        .Distinct()
                        .OrderBy(k => k)
                        .ToList();

                    var nestedCategories = new List<string>();
                    var categoryToIndex = new Dictionary<string, int>();
                    int categoryIndex = 0;

                    foreach (var groupKey in allGroupKeys)
                    {
                        foreach (var binKey in binToFilteredIndex.Keys.OrderBy(k => k))
                        {
                            string binLabel = Math.Round((binKey * binSize), 2).ToString(CultureInfo.InvariantCulture);
                            nestedCategories.Add(binLabel);
                            categoryToIndex[binKey.ToString() + "|" + groupKey] = categoryIndex++;
                        }
                    }

                    foreach (string sourceFile in psmsBySourceFile.Keys)
                    {
                        var column = new ColumnSeries
                        {
                            ColumnWidth = 200,
                            IsStacked = true,
                            Title = sourceFile,
                            TrackerFormatString = "Bin: {bin}\n{0}: {2}\nGroup: {group}\nTotal: {total}",
                            BaseValue = parameters.UseLogScaleYAxis ? 0.1 : 0
                        };

                        foreach (string groupKey in allGroupKeys)
                        {
                            if (!groupedPsmsBySourceFile[sourceFile].ContainsKey(groupKey))
                            {
                                foreach (var binKey in binToFilteredIndex.Keys.OrderBy(k => k))
                                {
                                    string lookupKey = binKey.ToString() + "|" + groupKey;
                                    if (categoryToIndex.ContainsKey(lookupKey))
                                    {
                                        column.Items.Add(new HistItem(0, categoryToIndex[lookupKey], 
                                            (binKey * binSize).ToString(CultureInfo.InvariantCulture), totalCounts[binToFilteredIndex[binKey]], groupKey));
                                    }
                                }
                                continue;
                            }

                            var groupPsms = groupedPsmsBySourceFile[sourceFile][groupKey];
                            var groupNumbers = GetNumbersFromPsms(groupPsms, plotType);
                            var groupResults = groupNumbers.GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key);
                            var groupDict = groupResults.ToDictionary(p => p.Key.ToString(), v => v.Count());

                            foreach (var binKey in binToFilteredIndex.Keys.OrderBy(k => k))
                            {
                                string lookupKey = binKey.ToString() + "|" + groupKey;
                                if (!categoryToIndex.ContainsKey(lookupKey))
                                    continue;

                                double total = 1.0;
                                if (parameters.NormalizeHistogramToFile)
                                {
                                    total = groupDict.Values.Sum();
                                }
                                
                                double value = groupDict.ContainsKey(binKey.ToString()) ? groupDict[binKey.ToString()] / total : 0;
                                
                                if (parameters.UseLogScaleYAxis && value > 0 && value < 0.1)
                                    value = 0.1;
                                
                                column.Items.Add(new HistItem(value, categoryToIndex[lookupKey], 
                                    (binKey * binSize).ToString(CultureInfo.InvariantCulture), totalCounts[binToFilteredIndex[binKey]], groupKey));
                            }
                        }
                        privateModel.Series.Add(column);
                    }

                    category = nestedCategories.ToArray();
                    
                    // Create group axis labels with proper spacing
                    var groupAxisLabels = new List<string>();
                    for (int i = 0; i < allGroupKeys.Count; i++)
                    {
                        for (int j = 0; j < binToFilteredIndex.Count; j++)
                        {
                            if (j == binToFilteredIndex.Count / 2)
                            {
                                groupAxisLabels.Add(allGroupKeys[i]);
                            }
                            else
                            {
                                groupAxisLabels.Add("");
                            }
                        }
                    }
                    
                    var groupAxis = new CategoryAxis
                    {
                        Position = AxisPosition.Bottom,
                        Title = parameters.GroupingProperty,
                        Key = "GroupAxis",
                        StartPosition = 1,
                        EndPosition = 0,
                        IsAxisVisible = true,
                        ItemsSource = groupAxisLabels.ToArray(),
                        FontSize = 12,
                        AxisTickToLabelDistance = 5,
                        MajorTickSize = 10
                    };
                    privateModel.Axes.Add(groupAxis);
                }
                else
                {
                    foreach (string key in dictsBySourceFile.Keys)
                    {
                        ColumnSeries column = new ColumnSeries { 
                            ColumnWidth = 200, 
                            IsStacked = true, 
                            Title = key, 
                            TrackerFormatString = "Bin: {bin}\n{0}: {2}\nTotal: {total}",
                            BaseValue = parameters.UseLogScaleYAxis ? 0.1 : 0
                        };
                        
                        foreach (var d in dictsBySourceFile[key])
                        {
                            int bin = int.Parse(d.Key);
                            
                            if (!binToFilteredIndex.ContainsKey(bin))
                                continue;

                            int filteredIdx = binToFilteredIndex[bin];
                            double total = 1.0;
                            if (parameters.NormalizeHistogramToFile)
                            {
                                total = dictsBySourceFile[key].Values.Sum();
                            }
                            
                            double value = d.Value / total;
                            
                            if (parameters.UseLogScaleYAxis && value < 0.1)
                                value = 0.1;
                        
                            column.Items.Add(new HistItem(value, filteredIdx, 
                                (bin * binSize).ToString(CultureInfo.InvariantCulture), totalCounts[filteredIdx]));
                        }
                        privateModel.Series.Add(column);
                    }
                }
            }

            // add axes
            if (parameters.NormalizeHistogramToFile)
                xAxisTitle = $"File Normalized {xAxisTitle}";
            
            var mainAxis = new CategoryAxis
            {
                Position = AxisPosition.Bottom,
                ItemsSource = category,
                Title = xAxisTitle,
                GapWidth = 0.3,
                Angle = labelAngle,
                Key = "MainAxis"
            };
            
            // Add spacing between groups when grouping is enabled
            if (isGroupingEnabled)
            {
                mainAxis.GapWidth = 0.1;  // Tighter spacing within groups
                
                // Calculate positions for group separators
                var allGroupKeys = groupedPsmsBySourceFile?.Values
                    .SelectMany(d => d.Keys)
                    .Distinct()
                    .OrderBy(k => k)
                    .ToList();
                
                if (allGroupKeys != null && allGroupKeys.Count > 1)
                {
                    int categoriesPerGroup = category.Length / allGroupKeys.Count;
                    
                    // Add extra gap positions between groups
                    mainAxis.ExtraGridlines = new double[allGroupKeys.Count - 1];
                    for (int i = 1; i < allGroupKeys.Count; i++)
                    {
                        mainAxis.ExtraGridlines[i - 1] = i * categoriesPerGroup - 0.5;
                    }
                    mainAxis.ExtraGridlineStyle = LineStyle.Solid;
                    mainAxis.ExtraGridlineColor = OxyColors.LightGray;
                    mainAxis.ExtraGridlineThickness = 2;
                }
            }
            
            privateModel.Axes.Add(mainAxis);

            if (parameters.UseLogScaleYAxis)
            {
                privateModel.Axes.Add(new LogarithmicAxis
                {
                    Title = yAxisTitle,
                    Position = AxisPosition.Left,
                    AbsoluteMinimum = 0.1,
                    Minimum = 0.1,
                    Base = 10
                });
            }
            else
            {
                privateModel.Axes.Add(new LinearAxis
                {
                    Title = yAxisTitle,
                    Position = AxisPosition.Left,
                    AbsoluteMinimum = 0,
                    Minimum = 0
                });
            }
        }

        private void linePlot(int plotType)
        {
            string yAxisTitle = "";
            string xAxisTitle = "";
            ScatterSeries series = new ScatterSeries
            {
                MarkerFill = OxyColors.Blue,
                MarkerSize = 0.5,
                TrackerFormatString = "{1}: {2:0.###}\n{3}: {4:0.###}\nFull sequence: {Tag}"
            };
            ScatterSeries variantSeries = new ScatterSeries
            {
                MarkerFill = OxyColors.DarkRed,
                MarkerSize = 1.5,
                MarkerType = MarkerType.Circle,
                TrackerFormatString = "{1}: {2:0.###}\n{3}: {4:0.###}\nFull sequence: {Tag}"
            };
            List<Tuple<double, double, string>> xy = new List<Tuple<double, double, string>>();
            List<Tuple<double, double, string>> variantxy = new List<Tuple<double, double, string>>();
            var filteredList = allSpectralMatches.Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa, CultureInfo.InvariantCulture), 0) == 0).ToList();
            var test = allSpectralMatches.SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm));
            switch (plotType)
            {
                case 1: // Precursor PPM Error vs. RT
                    yAxisTitle = "Precursor error (ppm)";
                    xAxisTitle = "Retention time";
                    foreach (var psm in filteredList)
                    {
                        if (psm.IdentifiedSequenceVariations == null || psm.IdentifiedSequenceVariations.Equals(""))
                        {
                            xy.Add(new Tuple<double, double, string>(double.Parse(psm.MassDiffPpm, CultureInfo.InvariantCulture), (double)psm.RetentionTime, psm.FullSequence));
                        }
                        else
                        {
                            variantxy.Add(new Tuple<double, double, string>(double.Parse(psm.MassDiffPpm, CultureInfo.InvariantCulture), (double)psm.RetentionTime, psm.FullSequence));
                        }
                    }
                    break;
                case 3: // Predicted RT vs. Observed RT
                    yAxisTitle = "Predicted Hydrophobicity";
                    xAxisTitle = "Observed retention time";
                    SSRCalc3 sSRCalc3 = new SSRCalc3("A100", SSRCalc3.Column.A100);
                    foreach (var psm in allSpectralMatches)
                    {
                        if (psm.IdentifiedSequenceVariations == null || psm.IdentifiedSequenceVariations.Equals(""))
                        {
                            xy.Add(new Tuple<double, double, string>(sSRCalc3.ScoreSequence(new PeptideWithSetModifications(psm.BaseSeq.Split('|')[0], null)),
                            (double)psm.RetentionTime, psm.FullSequence));
                        }
                        else
                        {
                            variantxy.Add(new Tuple<double, double, string>(sSRCalc3.ScoreSequence(new PeptideWithSetModifications(psm.BaseSeq.Split('|')[0], null)),
                            (double)psm.RetentionTime, psm.FullSequence));
                        }
                    }
                    break;
            }
            if (xy.Count != 0)
            {
                // plot each peptide
                IOrderedEnumerable<Tuple<double, double, string>> sorted = xy.OrderBy(x => x.Item1);
                foreach (var val in sorted)
                {
                    series.Points.Add(new ScatterPoint(val.Item2, val.Item1, tag: val.Item3));
                }
                privateModel.Series.Add(series);

                // add series displayed in legend, the real series will show up with a tiny dot for the symbol
                privateModel.Series.Add(new ScatterSeries { Title = "non-variant PSMs", MarkerFill = OxyColors.Blue });
            }

            if (variantxy.Count != 0)
            {
                // plot each variant peptide
                IOrderedEnumerable<Tuple<double, double, string>> variantSorted = variantxy.OrderBy(x => x.Item1);
                foreach (var val in variantSorted)
                {
                    variantSeries.Points.Add(new ScatterPoint(val.Item2, val.Item1, tag: val.Item3));
                }
                privateModel.Series.Add(variantSeries);

                // add series displayed in legend, the real series will show up with a tiny dot for the symbol
                privateModel.Series.Add(new ScatterSeries { Title = "variant PSMs", MarkerFill = OxyColors.DarkRed });
            }
            privateModel.Axes.Add(new LinearAxis { Title = xAxisTitle, Position = AxisPosition.Bottom });
            privateModel.Axes.Add(new LinearAxis { Title = yAxisTitle, Position = AxisPosition.Left });
        }

        // returns a bin index of number relative to 0, midpoints are rounded towards zero
        private static int roundToBin(double number, double binSize)
        {
            int sign = number < 0 ? -1 : 1;
            double d = number * sign;
            double remainder = d % binSize;
            int i = remainder < 0.5 * binSize ? (int)(d / binSize + 0.001) : (int)(d / binSize + 1.001);
            return i * sign;
        }

        /// <summary>
        /// Groups PSMs by the specified property for side-by-side plotting
        /// </summary>
        private Dictionary<string, Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>>> GroupPsmsByProperty(
            Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> psmsByFile, string propertyName)
        {
            var result = new Dictionary<string, Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>>>();

            foreach (var sourceFile in psmsByFile.Keys)
            {
                var groupedPsms = new Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>>();

                foreach (var psm in psmsByFile[sourceFile])
                {
                    string groupKey = GetGroupKeyFromPsm(psm, propertyName);
                    
                    if (!groupedPsms.ContainsKey(groupKey))
                        groupedPsms[groupKey] = new ObservableCollection<SpectrumMatchFromTsv>();
                    
                    groupedPsms[groupKey].Add(psm);
                }

                result[sourceFile] = groupedPsms;
            }

            return result;
        }

        /// <summary>
        /// Extracts the grouping property value from a PSM
        /// </summary>
        private string GetGroupKeyFromPsm(SpectrumMatchFromTsv psm, string propertyName)
        {
            return propertyName switch
            {
                "Notch" => psm.Notch ?? "0",
                "Precursor Charge" => psm.PrecursorCharge.ToString(),
                "File Name" => psm.FileNameWithoutExtension,
                "Ambiguity Level" => psm.AmbiguityLevel ?? "1",
                "Missed Cleavages" => psm.MissedCleavage ?? "0",
                _ => "All"
            };
        }

        /// <summary>
        /// Extracts numerical values from PSMs based on plot type
        /// </summary>
        private IEnumerable<double> GetNumbersFromPsms(ObservableCollection<SpectrumMatchFromTsv> psms, int plotType)
        {
            return plotType switch
            {
                1 => psms.Where(p => !p.MassDiffDa.Contains("|") && Math.Round(double.Parse(p.MassDiffDa, CultureInfo.InvariantCulture), 0) == 0)
                        .Select(p => double.Parse(p.MassDiffPpm, CultureInfo.InvariantCulture)),
                2 => psms.SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm)),
                3 => psms.Select(p => (double)p.PrecursorCharge),
                4 => psms.SelectMany(p => p.MatchedIons.Select(v => (double)v.Charge)),
                6 => psms.Select(p => (double)p.PrecursorMass),
                7 => psms.Select(p => (double)p.PrecursorMz),
                8 => GetHydrophobicityScores(psms),
                9 => psms.Where(p => !p.MissedCleavage.Contains("|")).Select(p => double.Parse(p.MissedCleavage)),
                12 => psms.Select(p => (double)(int)Math.Round(p.RetentionTime, 0)),
                _ => Enumerable.Empty<double>()
            };
        }

        private IEnumerable<double> GetHydrophobicityScores(ObservableCollection<SpectrumMatchFromTsv> psms)
        {
            SSRCalc3 sSRCalc3 = new SSRCalc3("A100", SSRCalc3.Column.A100);
            var values = new List<double>();
            foreach (var psm in psms.Where(p => p is not OsmFromTsv))
            {
                values.Add(sSRCalc3.ScoreSequence(new PeptideWithSetModifications(psm.BaseSeq.Split("|")[0], null)));
            }
            return values;
        }

        // used by histogram plots, gives additional properies for the tracker to display
        private class HistItem : ColumnItem
        {
            public int total { get; set; }
            public string bin { get; set; }
            public string group { get; set; }
            
            public HistItem(double value, int categoryIndex, string bin, int total, string group = null) : base(value, categoryIndex)
            {
                this.total = total;
                this.bin = bin;
                this.group = group;
            }
        }
    }
}