using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.SpectralLibrarySearch
{
    public static class SpectralLibraryReader
    {

        public static Dictionary<string, Book> ReadSpectralLibrary(string filePath)
        {
            Dictionary<string, Book> spectralLibraryDictionary = new Dictionary<string, Book>();
            string[] lines;

            lines = File.ReadAllLines(filePath);

            //find the lines which contain "name"
            var nameLine = new List<int>();
            for (int i = 0; i < lines.Length; i++)
            {
                if (lines[i].Contains("name", StringComparison.OrdinalIgnoreCase))
                {
                    nameLine.Add(i);
                }
            }
            nameLine.Add(lines.Length);// for the convenience to separate the file to different parts

            //for each spectrum 
            for (int i = 0; i < nameLine.Count - 1; i++)
            {
                string sequence = "";
                int z = 1;
                double precursorMz = 0;
                double rt = 0;
                List<MatchedFragmentIon> matchedFragmentIons = new List<MatchedFragmentIon>();
                //for each line
                for (int j = nameLine[i]; j < nameLine[i + 1] - 1; j++)
                {
                    //get name of each spectrum
                    if (lines[j].Contains("name", StringComparison.OrdinalIgnoreCase))
                    {
                        string[] name = lines[j].Split(new char[] { ':', '=' }, 2).Select(b => b.Trim()).ToArray();
                        string sequenceWithCharge = name[1];
                        string[] sequenceAndCharge = sequenceWithCharge.Split(new char[] { '/' }, 2).Select(b => b.Trim()).ToArray();
                        sequence = sequenceAndCharge[0];
                        if (sequenceAndCharge.Length > 1)
                        {
                            z = Convert.ToInt32(sequenceAndCharge[1]);
                        }
                    }

                    //get MW of each spectrum
                    else if ((lines[j].Contains("MW", StringComparison.OrdinalIgnoreCase) ||
                        lines[j].Contains("Monoisotopic Mass", StringComparison.OrdinalIgnoreCase))
                        && !lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                    {
                        double mw = Convert.ToDouble(lines[j].Split(":", 2).Select(b => b.Trim()).ToArray()[1]);
                    }

                    // get information from comment
                    if (lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                    {
                        string[] comment = lines[j].Split(" ").Select(b => b.Trim()).ToArray();
                        for (int l = 0; l < comment.Length; l++)
                        {
                            if (comment[l].Contains("parent", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("precursor", StringComparison.OrdinalIgnoreCase))
                            {
                                precursorMz = Convert.ToDouble(comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray()[1]);
                            }

                            if (comment[l].Contains("iRT", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("retention time", StringComparison.OrdinalIgnoreCase))
                            {
                                rt = Convert.ToDouble(comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray()[1]);
                            }
                        }
                    }

                    else if (lines[j].Contains("peaks", StringComparison.OrdinalIgnoreCase))
                    {
                        int numberOfPeaks = Convert.ToInt32(lines[j].Split(":").Select(b => b.Trim()).ToArray()[1]);

                        //load each peak 
                        for (int k = j + 1; k < j + 1 + numberOfPeaks; k++)
                        {
                            string[] eachPeak = lines[k].Split("\t").Select(b => b.Trim()).ToArray();
                            var experMz = double.Parse(eachPeak[0]);
                            var experIntensity = double.Parse(eachPeak[1]);
                            string[] ions = eachPeak[2].Split(new char[] { '/', '\"', 'p' }, StringSplitOptions.RemoveEmptyEntries).Select(b => b.Trim()).ToArray();
                            var spectrumPeakProductType = ions[0].ToCharArray()[0].ToString();
                            var fragmentNumber = (int)char.GetNumericValue(ions[0].ToCharArray()[1]);
                            int ionCharge = 1;
                            if (ions[0].ToCharArray().Length > 3)
                            {
                                ionCharge = (int)char.GetNumericValue(ions[0].ToCharArray()[3]);
                            }                        

                            ProductType peakProductType = (ProductType)Enum.Parse(typeof(ProductType), spectrumPeakProductType, true);
                            FragmentationTerminus terminus = (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus), "None", true);
                            var product = new Product(peakProductType, terminus, experMz, fragmentNumber, 0, 0);
                            matchedFragmentIons.Add(new MatchedFragmentIon(ref product, experMz, experIntensity, ionCharge));
                        }
                    }
                }
                spectralLibraryDictionary.Add(sequence + z, new Book(sequence, precursorMz, z, matchedFragmentIons));
            }

            return spectralLibraryDictionary;
        }
    }
}