// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace Development.Dia
{
    /// <summary>
    /// Parsed entry from a .msp spectral library file.
    /// </summary>
    public readonly struct MspEntry
    {
        public string Sequence { get; }
        public int Charge { get; }
        public double PrecursorMz { get; }
        public double? RetentionTime { get; }
        public bool IsDecoy { get; }
        public float[] FragmentMzs { get; }
        public float[] FragmentIntensities { get; }

        public MspEntry(string sequence, int charge, double precursorMz, double? rt,
            bool isDecoy, float[] fragmentMzs, float[] fragmentIntensities)
        {
            Sequence = sequence;
            Charge = charge;
            PrecursorMz = precursorMz;
            RetentionTime = rt;
            IsDecoy = isDecoy;
            FragmentMzs = fragmentMzs;
            FragmentIntensities = fragmentIntensities;
        }
    }

    /// <summary>
    /// Parses .msp spectral library files into MspEntry structs.
    ///
    /// Expected format per entry:
    ///   Name: SEQUENCE/CHARGE
    ///   MW: precursorMz
    ///   Comment: Parent=precursorMz RT=retentionTime
    ///   Num peaks: N
    ///   mz\tintensity\t"annotation"   (× N lines)
    ///
    /// This is a standalone parser so DiaDevBenchmarks has no dependency on the
    /// Readers.SpectralLibrary stack (which requires full mzLib Readers initialization).
    /// </summary>
    public static class MspParser
    {
        public static List<MspEntry> Parse(string path, bool isDecoy)
        {
            var entries = new List<MspEntry>();

            string currentSequence = null;
            int currentCharge = 0;
            double currentPrecursorMz = 0;
            double? currentRt = null;
            var mzList = new List<float>(64);
            var intList = new List<float>(64);

            foreach (string rawLine in File.ReadLines(path))
            {
                string line = rawLine.TrimEnd('\r', '\n');

                if (line.StartsWith("Name: ", StringComparison.Ordinal))
                {
                    // Flush previous entry
                    if (currentSequence != null && mzList.Count > 0)
                    {
                        entries.Add(new MspEntry(
                            currentSequence, currentCharge, currentPrecursorMz,
                            currentRt, isDecoy, mzList.ToArray(), intList.ToArray()));
                    }

                    string nameValue = line.Substring(6).Trim();
                    int slashIdx = nameValue.LastIndexOf('/');
                    if (slashIdx > 0 && int.TryParse(nameValue.AsSpan(slashIdx + 1), out int z))
                    {
                        currentSequence = nameValue.Substring(0, slashIdx);
                        currentCharge = z;
                    }
                    else
                    {
                        currentSequence = nameValue;
                        currentCharge = 2;
                    }

                    currentPrecursorMz = 0;
                    currentRt = null;
                    mzList.Clear();
                    intList.Clear();
                }
                else if (line.StartsWith("MW: ", StringComparison.Ordinal) ||
                         line.StartsWith("PRECURSORMZ: ", StringComparison.Ordinal))
                {
                    int colonIdx = line.IndexOf(':');
                    if (colonIdx >= 0)
                        double.TryParse(line.AsSpan(colonIdx + 1).Trim(stackalloc char[] { ' ', '\t' }),
                            NumberStyles.Float, CultureInfo.InvariantCulture, out currentPrecursorMz);
                }
                else if (line.StartsWith("Comment:", StringComparison.Ordinal))
                {
                    ReadOnlySpan<char> comment = line.AsSpan(8).Trim();

                    int parentIdx = comment.IndexOf("Parent=".AsSpan(), StringComparison.Ordinal);
                    if (parentIdx >= 0)
                    {
                        var afterParent = comment.Slice(parentIdx + 7);
                        int endIdx = afterParent.IndexOfAny(' ', '\t');
                        var parentStr = endIdx >= 0 ? afterParent.Slice(0, endIdx) : afterParent;
                        if (double.TryParse(parentStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double p))
                            currentPrecursorMz = p;
                    }

                    int rtIdx = comment.IndexOf("RT=".AsSpan(), StringComparison.Ordinal);
                    if (rtIdx >= 0)
                    {
                        var afterRt = comment.Slice(rtIdx + 3);
                        int endIdx = afterRt.IndexOfAny(' ', '\t');
                        var rtStr = endIdx >= 0 ? afterRt.Slice(0, endIdx) : afterRt;
                        if (double.TryParse(rtStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double rt))
                            currentRt = rt;
                    }
                }
                else if (line.Length > 0 && char.IsDigit(line[0]))
                {
                    // Peak line: mz\tintensity\t...
                    int firstTab = line.IndexOf('\t');
                    if (firstTab > 0)
                    {
                        int secondTab = line.IndexOf('\t', firstTab + 1);
                        int intEnd = secondTab > 0 ? secondTab : line.Length;

                        if (float.TryParse(line.AsSpan(0, firstTab),
                                NumberStyles.Float, CultureInfo.InvariantCulture, out float mz) &&
                            float.TryParse(line.AsSpan(firstTab + 1, intEnd - firstTab - 1),
                                NumberStyles.Float, CultureInfo.InvariantCulture, out float intensity))
                        {
                            mzList.Add(mz);
                            intList.Add(intensity);
                        }
                    }
                }
            }

            // Flush last entry
            if (currentSequence != null && mzList.Count > 0)
            {
                entries.Add(new MspEntry(
                    currentSequence, currentCharge, currentPrecursorMz,
                    currentRt, isDecoy, mzList.ToArray(), intList.ToArray()));
            }

            return entries;
        }
    }
}