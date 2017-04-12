﻿using EngineLayer;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public static class GlobalTaskLevelSettings
    {

        #region Public Constructors

        static GlobalTaskLevelSettings()
        {
            SearchModesKnown = LoadSearchModesFromFile().ToList();
            ProteaseDictionary = LoadProteaseDictionary();
            AllModsKnown = new List<Modification>();
        }

        #endregion Public Constructors

        #region Public Properties

        public static Dictionary<string, Protease> ProteaseDictionary { get; }
        public static List<SearchMode> SearchModesKnown { get; set; }
        public static List<Modification> AllModsKnown { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static void AddMods(IEnumerable<ModificationWithLocation> enumerable)
        {
            foreach (var ye in enumerable)
            {
                if (string.IsNullOrEmpty(ye.modificationType) || string.IsNullOrEmpty(ye.id))
                    throw new Exception(ye.ToString() + Environment.NewLine + " has null or empty modification type");
                if (AllModsKnown.Any(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType) && !b.Equals(ye)))
                    throw new Exception(ye.ToString() + Environment.NewLine + " has same and id and modification type as " + Environment.NewLine + AllModsKnown.First(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType)));
                else if (AllModsKnown.Any(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType)))
                    continue;
                else
                    AllModsKnown.Add(ye);
            }
        }

        #endregion Public Methods

        #region Private Methods

        private static Dictionary<string, Protease> LoadProteaseDictionary()
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();
            using (StreamReader proteases = new StreamReader(Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Data", "proteases.tsv")))
            {
                proteases.ReadLine();

                while (proteases.Peek() != -1)
                {
                    string line = proteases.ReadLine();
                    string[] fields = line.Split('\t');

                    string name = fields[0];
                    string[] sequences_inducing_cleavage = fields[1].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    string[] sequences_preventing_cleavage = fields[2].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    var cleavage_terminus = (TerminusType)Enum.Parse(typeof(TerminusType), fields[3], true);
                    var cleavage_specificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true);
                    string psi_ms_accession_number = fields[5];
                    string psi_ms_name = fields[6];
                    string site_regexp = fields[7];
                    var protease = new Protease(name, sequences_inducing_cleavage, sequences_preventing_cleavage, cleavage_terminus, cleavage_specificity, psi_ms_accession_number, psi_ms_name, site_regexp);
                    dict.Add(protease.Name, protease);
                }
            }
            return dict;
        }

        private static IEnumerable<SearchMode> LoadSearchModesFromFile()
        {
            yield return new SinglePpmAroundZeroSearchMode(10);
            yield return new SinglePpmAroundZeroSearchMode(5);
            yield return new SinglePpmAroundZeroSearchMode(20);
            yield return new SingleAbsoluteAroundZeroSearchMode(0.05);
            yield return new DotSearchMode("3mm", new double[] { 0, 1.003, 2.006, 3.009 }, new Tolerance(ToleranceUnit.PPM, 5));
            yield return new IntervalSearchMode("3.5aroundZero", new List<DoubleRange>() { new DoubleRange(-3.5, 3.5) });
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