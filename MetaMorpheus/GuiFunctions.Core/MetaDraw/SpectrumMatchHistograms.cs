using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using Readers;

namespace GuiFunctions
{
    /// <summary>
    /// The binning behind MetaDraw's statistics plots, as plain data rather than an OxyPlot model.
    /// Extracted from PlotModelStat so a second front end can render the same numbers instead of
    /// recomputing them - the two cannot then disagree about a count.
    /// </summary>
    public static class SpectrumMatchHistograms
    {

        public record HistogramRawData(
            string XAxisTitle, string YAxisTitle, double BinSize, double LabelAngle,
            Dictionary<string, IEnumerable<double>> NumbersBySourceFile,
            Dictionary<string, Dictionary<string, int>> DictsBySourceFile);

        public static HistogramRawData GetHistogramData(int plotType,
            Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> psmsBySourceFile)
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
                        var results = numbersBySourceFile[key].GroupBy(p => RoundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(key, results.ToDictionary(p => p.Key.ToString(), v => v.Count()));
                    }
                    break;
                case 2: // Histogram of Fragment PPM Errors 
                    xAxisTitle = "Fragment error (ppm)";
                    binSize = 2;
                    foreach (string key in psmsBySourceFile.Keys)
                    {
                        numbersBySourceFile.Add(key, psmsBySourceFile[key].SelectMany(p => p.MatchedIons.Select(v => v.MassErrorPpm)));
                        var results = numbersBySourceFile[key].GroupBy(p => RoundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
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
                        var results = numbersBySourceFile[key].GroupBy(p => RoundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
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
                        var results = numbersBySourceFile[key].GroupBy(p => RoundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
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
                        var results = numbersBySourceFile[key].GroupBy(p => RoundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
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
                        var results = numbersBySourceFile[fileName].GroupBy(p => RoundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
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
                case 13: // Histogram of Spectral Match Ambiguity Levels
                    xAxisTitle = "Ambiguity Level";
                    labelAngle = 0;
                    foreach (var fileName in psmsBySourceFile.Keys)
                    {
                        var result = psmsBySourceFile[fileName]
                            .GroupBy(p => NormalizeAmbiguityLevel(p.AmbiguityLevel))
                            .ToDictionary(p => p.Key, p => p.Count());
                        dictsBySourceFile.Add(fileName, result);
                    }
                    break;
                case 14: // Histogram of Notch (ambiguous PSMs contribute to each of their pipe-delimited notches)
                    xAxisTitle = "Notch";
                    binSize = 1;
                    labelAngle = 0;
                    foreach (var fileName in psmsBySourceFile.Keys)
                    {
                        var values = new List<double>();
                        foreach (var psm in psmsBySourceFile[fileName])
                        {
                            foreach (double notch in ParseAmbiguousNotchValues(psm.Notch))
                                values.Add(notch);
                        }
                        numbersBySourceFile.Add(fileName, values);
                        var results = values.GroupBy(p => RoundToBin(p, binSize)).OrderBy(p => p.Key).Select(p => p);
                        dictsBySourceFile.Add(fileName, results.ToDictionary(p => p.Key.ToString(CultureInfo.InvariantCulture), v => v.Count()));
                    }
                    break;
            }

            return new HistogramRawData(xAxisTitle, yAxisTitle, binSize, labelAngle, numbersBySourceFile, dictsBySourceFile);
        }


        public static string NormalizeAmbiguityLevel(string ambiguityLevel)
        {
            if (string.IsNullOrWhiteSpace(ambiguityLevel))
                return "1";

            var first = ambiguityLevel.Split('|')[0].Trim();
            return string.IsNullOrEmpty(first) ? "1" : first;
        }

        public static int RoundToBin(double number, double binSize)
        {
            int sign = number < 0 ? -1 : 1;
            double d = number * sign;
            double remainder = d % binSize;
            int i = remainder < 0.5 * binSize ? (int)(d / binSize + 0.001) : (int)(d / binSize + 1.001);
            return i * sign;
        }

        private static IEnumerable<double> ParseAmbiguousNotchValues(string notch)
        {
            foreach (var part in SplitAmbiguousNotch(notch))
            {
                if (double.TryParse(part, NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
                    yield return value;
                else
                    yield return 0.0;
            }
        }

        private static IEnumerable<string> SplitAmbiguousNotch(string notch)
        {
            if (string.IsNullOrWhiteSpace(notch))
                return new[] { "0" };
            return notch.Split('|', StringSplitOptions.RemoveEmptyEntries).Select(p => p.Trim());
        }
    }
}
