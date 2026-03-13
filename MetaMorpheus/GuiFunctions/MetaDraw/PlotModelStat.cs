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
using System.IO;
using System.Linq;
using System.Text;

namespace GuiFunctions
{
    public class PlotModelStat
    {
        private PlotModel privateModel;
        private readonly ObservableCollection<SpectrumMatchFromTsv> allSpectralMatches;
        private readonly Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> psmsBySourceFile;
        private readonly PlotModelStatParameters parameters;

        /// <summary>
        /// Stores the tabular data for text export. Each row is a dictionary of column name to value.
        /// </summary>
        public List<Dictionary<string, string>> PlotData { get; private set; } = new();

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

        /// <summary>
        /// Exports the plot data to a tab-separated text file.
        /// </summary>
        public void ExportToText(string filePath)
        {
            if (PlotData == null || PlotData.Count == 0)
                return;

            var columns = PlotData.SelectMany(r => r.Keys).Distinct().ToList();
            var sb = new StringBuilder();
            sb.AppendLine(string.Join("\t", columns));
            foreach (var row in PlotData)
            {
                sb.AppendLine(string.Join("\t", columns.Select(c => row.ContainsKey(c) ? row[c] : "")));
            }
            File.WriteAllText(filePath, sb.ToString());
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

            bool isGroupingEnabled = parameters.GroupingProperty != "None";
            Dictionary<string, Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>>> groupedPsmsBySourceFile = null;

            if (isGroupingEnabled)
            {
                groupedPsmsBySourceFile = GroupPsmsByProperty(psmsBySourceFile, parameters.GroupingProperty);
            }

            // Gather histogram data from PSMs
            var histData = GetHistogramData(plotType);
            string xAxisTitle = histData.XAxisTitle;
            string yAxisTitle = histData.YAxisTitle;
            double binSize = histData.BinSize;
            double labelAngle = histData.LabelAngle;
            var numbersBySourceFile = histData.NumbersBySourceFile;
            var dictsBySourceFile = histData.DictsBySourceFile;

            // Build series and categories
            string[] category;
            int[] totalCounts;
            int categoriesPerGroup = 0;
            List<string> allGroupKeys = null;
            bool isCategoryHistogram = plotType == 5 || plotType == 10 || plotType == 11;

            if (isCategoryHistogram)
            {
                (category, totalCounts, categoriesPerGroup, allGroupKeys) = BuildCategoryHistogramSeries(
                    plotType, dictsBySourceFile, groupedPsmsBySourceFile, isGroupingEnabled);
            }
            else
            {
                (category, totalCounts, categoriesPerGroup, allGroupKeys) = BuildNumericalHistogramSeries(
                    plotType, binSize, dictsBySourceFile, numbersBySourceFile, groupedPsmsBySourceFile, isGroupingEnabled);
            }

            // Configure axes
            ConfigureHistogramAxes(xAxisTitle, yAxisTitle, labelAngle, category,
                isGroupingEnabled, allGroupKeys, categoriesPerGroup);
        }

        #region Histogram Data Gathering

        private record HistogramRawData(
            string XAxisTitle, string YAxisTitle, double BinSize, double LabelAngle,
            Dictionary<string, IEnumerable<double>> NumbersBySourceFile,
            Dictionary<string, Dictionary<string, int>> DictsBySourceFile);

        private HistogramRawData GetHistogramData(int plotType)
        {
            string xAxisTitle = "";
            string yAxisTitle = "Count";
            double binSize = -1;
            double labelAngle = 0;
            var numbersBySourceFile = new Dictionary<string, IEnumerable<double>>();
            var dictsBySourceFile = new Dictionary<string, Dictionary<string, int>>();

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
                    binSize = 2;
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
                    if (dictsBySourceFile.Sum(p => p.Value.Keys.Count) >= 40)
                        labelAngle = -50;
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
                    if (dictsBySourceFile.Sum(p => p.Value.Keys.Count) >= 40)
                        labelAngle = -50;
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

            return new HistogramRawData(xAxisTitle, yAxisTitle, binSize, labelAngle, numbersBySourceFile, dictsBySourceFile);
        }

        #endregion

        #region Category Histogram Series

        private (string[] category, int[] totalCounts, int categoriesPerGroup, List<string> allGroupKeys)
            BuildCategoryHistogramSeries(int plotType,
                Dictionary<string, Dictionary<string, int>> dictsBySourceFile,
                Dictionary<string, Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>>> groupedPsmsBySourceFile,
                bool isGroupingEnabled)
        {
            IEnumerable<string> allCategories = dictsBySourceFile.Values.Select(p => p.Keys).SelectMany(p => p);
            Dictionary<string, int> categoryIDs = new Dictionary<string, int>();
            int counter = 0;
            foreach (string s in allCategories)
            {
                if (!categoryIDs.ContainsKey(s))
                    categoryIDs.Add(s, counter++);
            }
            var totalCounts = new int[counter];

            foreach (string cat in categoryIDs.Keys)
            {
                foreach (Dictionary<string, int> dict in dictsBySourceFile.Values)
                    totalCounts[categoryIDs[cat]] += dict.ContainsKey(cat) ? dict[cat] : 0;
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

            var category = filteredCategory.ToArray();
            var filteredCounts = filteredTotalCounts.ToArray();
            int categoriesPerGroup = filteredCategoryIDs.Count;
            List<string> allGroupKeys = null;

            if (isGroupingEnabled)
            {
                allGroupKeys = OrderByNaturalKey(
                    groupedPsmsBySourceFile.Values
                        .SelectMany(d => d.Keys)
                        .Distinct())
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
                                    column.Items.Add(new HistItem(0, categoryToIndex[lookupKey], cat, filteredCounts[filteredCategoryIDs[cat]], groupKey));
                            }
                            continue;
                        }

                        var groupDict = GetCategoryDictForGroup(plotType, groupedPsmsBySourceFile[sourceFile][groupKey]);

                        foreach (var cat in filteredCategoryIDs.Keys.OrderBy(k => k))
                        {
                            string lookupKey = cat + "|" + groupKey;
                            if (!categoryToIndex.ContainsKey(lookupKey))
                                continue;

                            double total = parameters.NormalizeHistogramToFile ? groupDict.Values.Sum() : 1.0;
                            int catId = filteredCategoryIDs[cat];
                            double value = groupDict.ContainsKey(cat) ? groupDict[cat] / total : 0;
                            if (parameters.UseLogScaleYAxis && value > 0 && value < 0.1)
                                value = 0.1;

                            column.Items.Add(new HistItem(value, categoryToIndex[lookupKey], cat, filteredCounts[catId], groupKey));

                            PlotData.Add(new Dictionary<string, string>
                            {
                                { "Source File", sourceFile }, { parameters.GroupingProperty, groupKey }, { "Category", cat },
                                { "Value", (groupDict.ContainsKey(cat) ? groupDict[cat] : 0).ToString() },
                                { "Total", filteredCounts[catId].ToString() }
                            });
                        }
                    }
                    privateModel.Series.Add(column);
                }

                category = nestedCategories.ToArray();
            }
            else
            {
                foreach (string key in dictsBySourceFile.Keys)
                {
                    ColumnSeries column = new ColumnSeries
                    {
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
                        double total = parameters.NormalizeHistogramToFile ? dictsBySourceFile[key].Values.Sum() : 1.0;
                        double value = d.Value / total;
                        if (parameters.UseLogScaleYAxis && value < 0.1)
                            value = 0.1;

                        column.Items.Add(new HistItem(value, id, d.Key, filteredCounts[id]));

                        PlotData.Add(new Dictionary<string, string>
                        {
                            { "Source File", key }, { "Category", d.Key },
                            { "Value", d.Value.ToString() }, { "Total", filteredCounts[id].ToString() }
                        });
                    }
                    privateModel.Series.Add(column);
                }
            }

            return (category, filteredCounts, categoriesPerGroup, allGroupKeys);
        }

        private Dictionary<string, int> GetCategoryDictForGroup(int plotType, ObservableCollection<SpectrumMatchFromTsv> groupPsms)
        {
            if (plotType == 5)
            {
                var psmsWithMods = groupPsms.Where(p => !p.FullSequence.Contains("|") && p.FullSequence.Contains("["));
                var mods = psmsWithMods.Select(p => p.ToBioPolymerWithSetMods()).Select(p => p.AllModsOneIsNterminus).SelectMany(p => p.Values);
                return mods.GroupBy(p => p.IdWithMotif).ToDictionary(p => p.Key, v => v.Count());
            }

            var allMatchedIons = groupPsms.SelectMany(p => p.MatchedIons).ToList();
            if (plotType == 10)
            {
                return allMatchedIons
                    .GroupBy(p => p.NeutralTheoreticalProduct is CustomMProduct cmp ? cmp.Annotation : p.NeutralTheoreticalProduct.ProductType.ToString())
                    .ToDictionary(p => p.Key.ToString(), p => p.Count());
            }

            return allMatchedIons
                .GroupBy(p => p.NeutralTheoreticalProduct is CustomMProduct cmp ? cmp.Annotation : p.NeutralTheoreticalProduct.ProductType.ToString())
                .ToDictionary(p => p.Key.ToString(), p => (int)p.Sum(m => m.Intensity));
        }

        #endregion

        #region Numerical Histogram Series

        private (string[] category, int[] totalCounts, int categoriesPerGroup, List<string> allGroupKeys)
            BuildNumericalHistogramSeries(int plotType, double binSize,
                Dictionary<string, Dictionary<string, int>> dictsBySourceFile,
                Dictionary<string, IEnumerable<double>> numbersBySourceFile,
                Dictionary<string, Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>>> groupedPsmsBySourceFile,
                bool isGroupingEnabled)
        {
            int end = dictsBySourceFile.Values.Max(p => p.Max(v => int.Parse(v.Key)));
            int start = dictsBySourceFile.Values.Min(p => p.Min(v => int.Parse(v.Key)));
            int numBins = end - start + 1;
            int minBinLabels = 22;
            int skipBinLabel = numBins < minBinLabels ? 1 : numBins / minBinLabels;

            var totalCounts = new int[numBins];
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
                    filteredCategories.Add(i % skipBinLabel == 0
                        ? Math.Round((i * binSize), 2).ToString(CultureInfo.InvariantCulture)
                        : "");
                    filteredTotalCounts.Add(totalCounts[binIdx]);
                    filteredIndex++;
                }
            }

            var category = filteredCategories.ToArray();
            var filteredCounts = filteredTotalCounts.ToArray();
            int categoriesPerGroup = binToFilteredIndex.Count;
            List<string> allGroupKeys = null;

            if (isGroupingEnabled)
            {
                allGroupKeys = OrderByNaturalKey(
                    groupedPsmsBySourceFile.Values
                        .SelectMany(d => d.Keys)
                        .Distinct())
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
                                    column.Items.Add(new HistItem(0, categoryToIndex[lookupKey],
                                        (binKey * binSize).ToString(CultureInfo.InvariantCulture), filteredCounts[binToFilteredIndex[binKey]], groupKey));
                            }
                            continue;
                        }

                        var groupPsms = groupedPsmsBySourceFile[sourceFile][groupKey];
                        var groupNumbers = GetNumbersFromPsms(groupPsms, plotType);
                        var groupDict = groupNumbers.GroupBy(p => roundToBin(p, binSize)).OrderBy(p => p.Key)
                            .ToDictionary(p => p.Key.ToString(), v => v.Count());

                        foreach (var binKey in binToFilteredIndex.Keys.OrderBy(k => k))
                        {
                            string lookupKey = binKey.ToString() + "|" + groupKey;
                            if (!categoryToIndex.ContainsKey(lookupKey))
                                continue;

                            double total = parameters.NormalizeHistogramToFile ? groupDict.Values.Sum() : 1.0;
                            int rawValue = groupDict.ContainsKey(binKey.ToString()) ? groupDict[binKey.ToString()] : 0;
                            double value = rawValue / total;
                            if (parameters.UseLogScaleYAxis && value > 0 && value < 0.1)
                                value = 0.1;

                            column.Items.Add(new HistItem(value, categoryToIndex[lookupKey],
                                (binKey * binSize).ToString(CultureInfo.InvariantCulture), filteredCounts[binToFilteredIndex[binKey]], groupKey));

                            PlotData.Add(new Dictionary<string, string>
                            {
                                { "Source File", sourceFile }, { parameters.GroupingProperty, groupKey },
                                { "Bin", (binKey * binSize).ToString(CultureInfo.InvariantCulture) },
                                { "Value", rawValue.ToString() }, { "Total", filteredCounts[binToFilteredIndex[binKey]].ToString() }
                            });
                        }
                    }
                    privateModel.Series.Add(column);
                }

                category = nestedCategories.ToArray();
            }
            else
            {
                foreach (string key in dictsBySourceFile.Keys)
                {
                    var column = new ColumnSeries
                    {
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
                        double total = parameters.NormalizeHistogramToFile ? dictsBySourceFile[key].Values.Sum() : 1.0;
                        double value = d.Value / total;
                        if (parameters.UseLogScaleYAxis && value < 0.1)
                            value = 0.1;

                        column.Items.Add(new HistItem(value, filteredIdx,
                            (bin * binSize).ToString(CultureInfo.InvariantCulture), filteredCounts[filteredIdx]));

                        PlotData.Add(new Dictionary<string, string>
                        {
                            { "Source File", key },
                            { "Bin", (bin * binSize).ToString(CultureInfo.InvariantCulture) },
                            { "Value", d.Value.ToString() }, { "Total", filteredCounts[filteredIdx].ToString() }
                        });
                    }
                    privateModel.Series.Add(column);
                }
            }

            return (category, filteredCounts, categoriesPerGroup, allGroupKeys);
        }

        #endregion

        #region Axis Configuration

        private void ConfigureHistogramAxes(string xAxisTitle, string yAxisTitle, double labelAngle,
            string[] category, bool isGroupingEnabled, List<string> allGroupKeys, int categoriesPerGroup)
        {
            if (parameters.NormalizeHistogramToFile)
                xAxisTitle = $"File Normalized {xAxisTitle}";

            double bottomPadding = isGroupingEnabled ? 20 : 0;
            privateModel.Padding = new OxyThickness(privateModel.Padding.Left, privateModel.Padding.Top,
                privateModel.Padding.Right, privateModel.Padding.Bottom + bottomPadding);

            var mainAxis = new CategoryAxis
            {
                Position = AxisPosition.Bottom,
                ItemsSource = category,
                Title = isGroupingEnabled ? null : xAxisTitle,
                GapWidth = 0.3,
                Angle = labelAngle,
            };

            if (isGroupingEnabled && allGroupKeys != null && allGroupKeys.Count > 0 && categoriesPerGroup > 0)
            {
                mainAxis.GapWidth = 0.1;

                if (allGroupKeys.Count > 1)
                {
                    mainAxis.ExtraGridlines = new double[allGroupKeys.Count - 1];
                    for (int i = 1; i < allGroupKeys.Count; i++)
                        mainAxis.ExtraGridlines[i - 1] = i * categoriesPerGroup - 0.5;
                    mainAxis.ExtraGridlineStyle = LineStyle.Solid;
                    mainAxis.ExtraGridlineColor = OxyColors.LightGray;
                    mainAxis.ExtraGridlineThickness = 2;
                }

                var labelledCategories = new string[category.Length];
                for (int g = 0; g < allGroupKeys.Count; g++)
                {
                    int startIdx = g * categoriesPerGroup;
                    for (int j = 0; j < categoriesPerGroup; j++)
                    {
                        int idx = startIdx + j;
                        if (idx >= category.Length) break;
                        labelledCategories[idx] = category[idx];
                    }
                }
                mainAxis.ItemsSource = labelledCategories;
                mainAxis.Title = xAxisTitle;

                var groupLabels = new string[category.Length];
                for (int g = 0; g < allGroupKeys.Count; g++)
                {
                    int startIdx = g * categoriesPerGroup;
                    int midIdx = startIdx + categoriesPerGroup / 2;
                    for (int j = 0; j < categoriesPerGroup; j++)
                    {
                        int idx = startIdx + j;
                        if (idx >= groupLabels.Length) break;
                        groupLabels[idx] = idx == midIdx ? allGroupKeys[g] : "";
                    }
                }

                var groupAxis = new CategoryAxis
                {
                    Position = AxisPosition.Bottom,
                    Key = "GroupAxis",
                    ItemsSource = groupLabels,
                    GapWidth = 0.1,
                    Angle = 0,
                    Title = parameters.GroupingProperty,
                    TitleFontWeight = FontWeights.Normal,
                    TitleFontSize = privateModel.DefaultFontSize + 2,
                    FontSize = privateModel.DefaultFontSize + 2,
                    FontWeight = FontWeights.Bold,
                    TickStyle = TickStyle.None,
                    IsAxisVisible = true,
                    AxislineStyle = LineStyle.None,
                    MajorGridlineStyle = LineStyle.None,
                    MinorGridlineStyle = LineStyle.None,
                    PositionTier = 1,
                };
                privateModel.Axes.Add(groupAxis);
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
                    Base = 10,
                    Key = "Primary"
                });
            }
            else
            {
                privateModel.Axes.Add(new LinearAxis
                {
                    Title = yAxisTitle,
                    Position = AxisPosition.Left,
                    AbsoluteMinimum = 0,
                    Minimum = 0,
                    Key = "Primary"
                });
            }
        }

        #endregion

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

                        PlotData.Add(new Dictionary<string, string>
                        {
                            { "Retention Time", psm.RetentionTime.ToString(CultureInfo.InvariantCulture) },
                            { "Precursor Error Ppm", psm.MassDiffPpm },
                            { "Full Sequence", psm.FullSequence },
                            { "Is Variant", (psm.IdentifiedSequenceVariations != null && !psm.IdentifiedSequenceVariations.Equals("")).ToString() }
                        });
                    }
                    break;
                case 3: // Predicted RT vs. Observed RT
                    yAxisTitle = "Predicted Hydrophobicity";
                    xAxisTitle = "Observed retention time";
                    SSRCalc3 sSRCalc3 = new SSRCalc3("A100", SSRCalc3.Column.A100);
                    foreach (var psm in allSpectralMatches)
                    {
                        double predicted = sSRCalc3.ScoreSequence(new PeptideWithSetModifications(psm.BaseSeq.Split('|')[0], null));
                        if (psm.IdentifiedSequenceVariations == null || psm.IdentifiedSequenceVariations.Equals(""))
                        {
                            xy.Add(new Tuple<double, double, string>(predicted, (double)psm.RetentionTime, psm.FullSequence));
                        }
                        else
                        {
                            variantxy.Add(new Tuple<double, double, string>(predicted, (double)psm.RetentionTime, psm.FullSequence));
                        }

                        PlotData.Add(new Dictionary<string, string>
                        {
                            { "Retention Time", psm.RetentionTime.ToString(CultureInfo.InvariantCulture) },
                            { "Predicted Hydrophobicity", predicted.ToString(CultureInfo.InvariantCulture) },
                            { "Full Sequence", psm.FullSequence },
                            { "Is Variant", (psm.IdentifiedSequenceVariations != null && !psm.IdentifiedSequenceVariations.Equals("")).ToString() }
                        });
                    }
                    break;
            }
            if (xy.Count != 0)
            {
                IOrderedEnumerable<Tuple<double, double, string>> sorted = xy.OrderBy(x => x.Item1);
                foreach (var val in sorted)
                {
                    series.Points.Add(new ScatterPoint(val.Item2, val.Item1, tag: val.Item3));
                }
                privateModel.Series.Add(series);
                privateModel.Series.Add(new ScatterSeries { Title = "non-variant PSMs", MarkerFill = OxyColors.Blue });
            }

            if (variantxy.Count != 0)
            {
                IOrderedEnumerable<Tuple<double, double, string>> variantSorted = variantxy.OrderBy(x => x.Item1);
                foreach (var val in variantSorted)
                {
                    variantSeries.Points.Add(new ScatterPoint(val.Item2, val.Item1, tag: val.Item3));
                }
                privateModel.Series.Add(variantSeries);
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
        /// Groups PSMs by the specified property for side-by-side plotting.
        /// Excludes ambiguous PSMs (those with "|" in FullSequence) since they cannot be
        /// cleanly assigned to a single group.
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

                    if (groupKey.Contains('|'))
                    {
                        if (MetaDrawSettings.AmbiguityFilter == "1")
                            continue;

                        if (!parameters.AllowAmbiguousGroups)
                            continue;

                        groupKey = NormalizeAmbiguousGroupKey(groupKey);
                    }

                    if (!groupedPsms.ContainsKey(groupKey))
                        groupedPsms[groupKey] = new ObservableCollection<SpectrumMatchFromTsv>();

                    groupedPsms[groupKey].Add(psm);
                }

                result[sourceFile] = groupedPsms;
            }

            return result;
        }

        /// <summary>
        /// Normalizes an ambiguous group key by sorting its pipe-delimited parts into a canonical order.
        /// This ensures that keys like "0|2|1" and "1|0|2" map to the same group ("0|1|2").
        /// Uses natural (numeric) sorting when all parts are numeric.
        /// </summary>
        private static string NormalizeAmbiguousGroupKey(string groupKey)
        {
            var parts = groupKey.Split('|');
            return string.Join("|", OrderByNaturalKey(parts));
        }

        private string GetGroupKeyFromPsm(SpectrumMatchFromTsv psm, string propertyName)
        {
            return propertyName switch
            {
                "Notch" => psm.Notch ?? "0",
                "Precursor Charge" => psm.PrecursorCharge.ToString(),
                "File Name" => psm.FileNameWithoutExtension,
                "Ambiguity Level" => psm.AmbiguityLevel ?? "1",
                "Missed Cleavages" => psm.MissedCleavage ?? "0",
                "SpectralAngle" => psm.SpectralAngle?.ToString() ?? "N/A",
                "PrecursorIntensity" => psm.PrecursorIntensity?.ToString() ?? "0",
                "QValue" => psm.QValue.ToString(),
                "PEP_QValue" => psm.PEP_QValue.ToString(),
                "TotalIonCurrent" => psm.TotalIonCurrent?.ToString() ?? "0",
                "OrganismName" => psm.OrganismName,
                "DecoyContamTarget" => psm.DecoyContamTarget,

                _ => "All"
            };
        }

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

        /// <summary>
        /// Orders strings numerically if all values are numeric, otherwise alphabetically.
        /// </summary>
        private static IOrderedEnumerable<string> OrderByNaturalKey(IEnumerable<string> keys)
        {
            bool allNumeric = keys.All(k => double.TryParse(k, NumberStyles.Any, CultureInfo.InvariantCulture, out _));
            return allNumeric
                ? keys.OrderBy(k => double.Parse(k, CultureInfo.InvariantCulture))
                : keys.OrderBy(k => k);
        }
    }
}