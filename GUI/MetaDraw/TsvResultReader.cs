using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text.RegularExpressions;
using System.Windows;
using MassSpectrometry;
using Chemistry;

namespace MetaMorpheusGUI
{
    public class TsvResultReader
    {
        private const string FULL_SEQUENCE_HEADER = "Full Sequence";
        private const string SCAN_NUMBER_HEADER = "Scan Number";
        private const string FILENAME_HEADER = "File Name";
        private const string MATCHED_MZ_HEADER = "Matched Ion Mass-To-Charge Ratios";
        private static Regex ionParser = new Regex(@"([a-zA-Z]+)(\d+)");
        private static char[] split = new char[] { '\t' };
        private static char[] mzSplit = new char[] { '[', ',', ']', ';' };

        public static List<MetaDrawPsm> ReadTsv(string filePath)
        {
            List<MetaDrawPsm> psms = new List<MetaDrawPsm>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filePath);
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not read the file: " + e.Message);
                return psms;
            }

            int lineCount = 0;
            int fullSeqIndex = -1;
            int scanNumberIndex = -1;
            int fileNameIndex = -1;
            int matchedMzIndex = -1;
            string line;
            while (reader.Peek() > 0)
            {
                lineCount++;

                line = reader.ReadLine();
                var spl = line.Split(split);

                if (lineCount == 1)
                {
                    fullSeqIndex = Array.IndexOf(spl, FULL_SEQUENCE_HEADER);
                    scanNumberIndex = Array.IndexOf(spl, SCAN_NUMBER_HEADER);
                    fileNameIndex = Array.IndexOf(spl, FILENAME_HEADER);
                    matchedMzIndex = Array.IndexOf(spl, MATCHED_MZ_HEADER);
                }

                if (fullSeqIndex < 0 || scanNumberIndex < 0 || fileNameIndex < 0 || matchedMzIndex < 0)
                {
                    MessageBox.Show("Could not read PSMs file. Is it from an older version of MetaMorpheus?");
                    return psms;
                }

                try
                {
                    int oneBasedScanNumber = int.Parse(spl[scanNumberIndex]);
                    string peptideSequence = spl[fullSeqIndex];
                    string fileName = spl[fileNameIndex];
                    string matchedPeaks = spl[matchedMzIndex];
                    List<TheoreticalFragmentIon> peaks = ReadFragmentIonsFromString(matchedPeaks);

                    psms.Add(new MetaDrawPsm(oneBasedScanNumber, fileName, peptideSequence, peaks));
                }
                catch (Exception)
                {
                    // TODO: write some kind of warning here?
                }
            }

            reader.Close();

            if ((lineCount - 1) != psms.Count)
            {
                MessageBox.Show("Warning: " + ((lineCount - 1) - psms.Count) + " PSMs were not read.");
            }

            return psms;
        }

        private static List<TheoreticalFragmentIon> ReadFragmentIonsFromString(string matchedMzString)
        {
            var peaks = matchedMzString.Split(mzSplit, StringSplitOptions.RemoveEmptyEntries).Select(v => v.Trim());
            List<TheoreticalFragmentIon> matchedIons = new List<TheoreticalFragmentIon>();

            foreach (var peak in peaks)
            {
                var split = peak.Split(new char[] { '+', ':' });

                string ionTypeAndNumber = split[0];
                Match result = ionParser.Match(ionTypeAndNumber);
                string ionType = result.Groups[1].Value;
                int ionNumber = int.Parse(result.Groups[2].Value);

                ProductType p = ProductType.None;
                if (ionType.Contains("b"))
                    p = ProductType.B;
                else if (ionType.Contains("y"))
                    p = ProductType.Y;
                else if (ionType.Contains("c"))
                    p = ProductType.C;
                else if (ionType.Contains("z"))
                    p = ProductType.Zdot;

                int z = int.Parse(split[1]);
                double mz = double.Parse(split[2]);

                matchedIons.Add(new TheoreticalFragmentIon(mz.ToMass(z), double.NaN, z, p, ionNumber));
            }

            return matchedIons;
        }
    }
}
