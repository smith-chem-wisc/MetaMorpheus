using EngineLayer;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace TaskLayer
{
    public static class GlobalTaskLevelSettings
    {

        #region Public Constructors

        static GlobalTaskLevelSettings()
        {
            SearchModesKnown = LoadSearchModesFromFile().ToList();
        }

        #endregion Public Constructors

        #region Public Properties

        public static List<SearchMode> SearchModesKnown { get; private set; }

        #endregion Public Properties

        #region Private Methods

        private static IEnumerable<SearchMode> LoadSearchModesFromFile()
        {
            yield return new SinglePpmAroundZeroSearchMode(10);
            yield return new SinglePpmAroundZeroSearchMode(5);
            yield return new SinglePpmAroundZeroSearchMode(20);
            yield return new SingleAbsoluteAroundZeroSearchMode(0.05);
            yield return new DotSearchMode("3mm", new double[] { 0, 1.003, 2.006, 3.009 }, new Tolerance(ToleranceUnit.PPM, 5));
            yield return new IntervalSearchMode("2.1aroundZero", new List<DoubleRange>() { new DoubleRange(-2.1, 2.1) });
            yield return new OpenSearchMode();
            yield return new IntervalSearchMode("ZeroAndSodium", new List<DoubleRange> { new DoubleRange(-0.005, 0.005), new DoubleRange(21.981943 - 0.005, 21.981943 + 0.005) });
            yield return new IntervalSearchMode("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });
            foreach (var sm in GetResidueInclusionExclusionSearchModes(new DoubleRange(-187, double.PositiveInfinity), 0.0075))
                yield return sm;
        }

        /// <summary>
        /// Ideally v is less than 0.00168565165, so no overlaps happen
        /// </summary>
        /// <param name="doubleRange"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static IEnumerable<SearchMode> GetResidueInclusionExclusionSearchModes(DoubleRange doubleRange, double v)
        {
            List<double> massesToExclude = new List<double>();
            for (char c = 'A'; c <= 'Z'; c++)
            {
                Residue residue;
                if (Residue.TryGetResidue(c, out residue))
                {
                    massesToExclude.Add(residue.MonoisotopicMass);
                    massesToExclude.Add(-residue.MonoisotopicMass);
                    for (char cc = 'A'; cc <= 'Z'; cc++)
                    {
                        Residue residueCC;
                        if (Residue.TryGetResidue(cc, out residueCC))
                        {
                            massesToExclude.Add(residue.MonoisotopicMass + residueCC.MonoisotopicMass);
                            massesToExclude.Add(residue.MonoisotopicMass - residueCC.MonoisotopicMass);
                            massesToExclude.Add(-residue.MonoisotopicMass + residueCC.MonoisotopicMass);
                            massesToExclude.Add(-residue.MonoisotopicMass - residueCC.MonoisotopicMass);
                        }
                    }
                }
            }
            List<double> filteredMasses = massesToExclude.GroupBy(b => Math.Round(b, 6)).Select(b => b.FirstOrDefault()).OrderBy(b => b).ToList();

            yield return new DotSearchMode("OnlyAAs", filteredMasses, new Tolerance(ToleranceUnit.Absolute, v));

            List<DoubleRange> doubleRanges = new List<DoubleRange>();

            var prevGoodMin = double.NegativeInfinity;

            for (int i = 0; i < filteredMasses.Count; i++)
            {
                if (Math.Round(filteredMasses[i], 6) == 0)
                    continue;
                doubleRanges.Add(new DoubleRange(prevGoodMin, filteredMasses[i] - v));
                prevGoodMin = filteredMasses[i] + v;
            }
            doubleRanges.Add(new DoubleRange(prevGoodMin, double.PositiveInfinity));

            doubleRanges = doubleRanges.Where(b => b.Minimum <= doubleRange.Maximum && b.Maximum >= doubleRange.Minimum).Select(b => new DoubleRange(Math.Max(doubleRange.Minimum, b.Minimum), Math.Min(doubleRange.Maximum, b.Maximum))).ToList();

            yield return new IntervalSearchMode("ExcludeAAs", doubleRanges);
        }

        #endregion Private Methods

    }
}