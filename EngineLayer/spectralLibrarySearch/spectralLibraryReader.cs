using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Globalization;

namespace EngineLayer
{
    public static class SpectralLibraryReader
    {
        public static Dictionary<string, LibrarySpectrum> ReadSpectralLibrary(string filePath)
        {
            Dictionary<string, LibrarySpectrum> spectralLibraryDictionary = new Dictionary<string, LibrarySpectrum>();
            string[] lines;

            lines = File.ReadAllLines(filePath);

            //find the lines which contain "name"
            var nameLines = new List<int>();
            for (int i = 0; i < lines.Length; i++)
            {
                if (lines[i].Contains("name", StringComparison.OrdinalIgnoreCase))
                {
                    nameLines.Add(i);
                }
            }
            nameLines.Add(lines.Length);// for the convenience to separate the file to different parts

            for (int i = 0; i < nameLines.Count - 1; i++)
            {
                string sequence = "";
                int z = 1;
                double precursorMz = 0;
                double rt = 0;
                List<MatchedFragmentIon> matchedFragmentIons = new List<MatchedFragmentIon>();

                for (int j = nameLines[i]; j < nameLines[i + 1] - 1; j++)
                {
                    // get name of each spectrum
                    if (lines[j].Contains("name", StringComparison.OrdinalIgnoreCase))
                    {
                        string[] name = lines[j].Split(new char[] { ':', '=' }, 2).Select(b => b.Trim()).ToArray();
                        string sequenceWithCharge = name[1];
                        string[] sequenceAndCharge = sequenceWithCharge.Split(new char[] { '/' }, 2).Select(b => b.Trim()).ToArray();
                        sequence = sequenceAndCharge[0];
                        if (sequenceAndCharge.Length > 1)
                        {
                            z = int.Parse(sequenceAndCharge[1]);
                        }
                    }

                    // get m/z of the peptide. "MW" is not the molecular weight but is the m/z.
                    else if ((lines[j].Contains("MW", StringComparison.OrdinalIgnoreCase) ||
                        lines[j].Contains("Monoisotopic Mass", StringComparison.OrdinalIgnoreCase))
                        && !lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                    {
                        double mw = double.Parse(lines[j].Split(":", 2).Select(b => b.Trim()).ToArray()[1], CultureInfo.InvariantCulture);
                    }

                    // get information from comment
                    if (lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                    {
                        string[] comment = lines[j].Split(" ").Select(b => b.Trim()).ToArray();
                        for (int l = 0; l < comment.Length; l++)
                        {
                            if (comment[l].Contains("parent", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("precursor", StringComparison.OrdinalIgnoreCase))
                            {
                                precursorMz = double.Parse(comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray()[1], CultureInfo.InvariantCulture);
                            }

                            if (comment[l].Contains("iRT", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("retention time", StringComparison.OrdinalIgnoreCase))
                            {
                                rt = double.Parse(comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray()[1], CultureInfo.InvariantCulture);
                            }
                        }
                    }

                    else if (lines[j].Contains("peaks", StringComparison.OrdinalIgnoreCase))
                    {
                        int numberOfPeaks = int.Parse(lines[j].Split(":").Select(b => b.Trim()).ToArray()[1]);

                        // read each peak 
                        for (int k = j + 1; k < j + 1 + numberOfPeaks; k++)
                        {
                            var peakLine = lines[k];

                            string[] peak = peakLine.Split("\t").Select(b => b.Trim()).ToArray();
                            var experMz = double.Parse(peak[0], CultureInfo.InvariantCulture);
                            var experIntensity = double.Parse(peak[1], CultureInfo.InvariantCulture);
                            string[] ionInfo = peak[2].Split(new char[] { '/', '\"', ')', '(' }, StringSplitOptions.RemoveEmptyEntries).Select(b => b.Trim()).ToArray();
                            ionInfo[1] = ionInfo[1].Replace("ppm", "", ignoreCase: true, CultureInfo.InvariantCulture);

                            //TODO: figure out a more robust way to do this
                            var spectrumPeakProductType = ionInfo[0].ToCharArray()[0].ToString();

                            int fragmentNumber = int.Parse(new string(ionInfo[0].Split(new char[] { '^' })[0].Where(Char.IsDigit).ToArray()));

                            int ionCharge = 1;
                            if (ionInfo[0].Contains('^'))
                            {
                                ionCharge = int.Parse(ionInfo[0].Split('^')[1]);
                            }

                            ProductType peakProductType = (ProductType)Enum.Parse(typeof(ProductType), spectrumPeakProductType, true);

                            //TODO: figure out terminus
                            FragmentationTerminus terminus = (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus), "None", true);

                            //TODO: figure out amino acid position
                            var product = new Product(peakProductType, terminus, experMz, fragmentNumber, 0, 0);

                            matchedFragmentIons.Add(new MatchedFragmentIon(ref product, experMz, experIntensity, ionCharge));
                        }
                    }
                }

                spectralLibraryDictionary.Add(sequence + z, new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt));
            }

            return spectralLibraryDictionary;
        }
    }
}