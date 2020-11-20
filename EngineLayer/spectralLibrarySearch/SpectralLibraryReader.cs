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
        private static Dictionary<string, string> PrositToMetaMorpheusModDictionary = new Dictionary<string, string>
        {
            { "Oxidation","[Common Variable:Oxidation on M]" },
            {"Carbamidomethyl", "[Common Fixed:Carbamidomethyl on C]" }
        };

        public static Dictionary<string, LibrarySpectrum> ReadSpectralLibrary(string filePath)
        {
            Dictionary<string, LibrarySpectrum> spectralLibraryDictionary = new Dictionary<string, LibrarySpectrum>();
            string[] lines;

            lines = File.ReadAllLines(filePath);

            //set up variables: 
            string sequence = "";
            int z = 0;
            double precursorMz = 0;
            double rt = 0;
            List<MatchedFragmentIon> matchedFragmentIons = new List<MatchedFragmentIon>();

            for (int lineIndex = 0; lineIndex < lines.Length; lineIndex++)
            {
                string line = lines[lineIndex];
                if (line.Substring(0, 4).Equals("Name"))
                {
                    //this is a new entry, save the old entry
                    if (lineIndex != 0)
                    {
                        LibrarySpectrum librarySpectrum = new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
                        spectralLibraryDictionary.Add(librarySpectrum.Name, librarySpectrum);
                    }

                    //start making new entry
                    string[] name = line.Substring(6).Split('/');
                    sequence = name[0];
                    z = Convert.ToInt32(name[1]);
                    matchedFragmentIons = new List<MatchedFragmentIon>();

                    //move to the next line for the MW
                    line = lines[++lineIndex];
                    precursorMz = Convert.ToDouble(line.Substring(4));

                    //move to the next line for the comment
                    line = lines[++lineIndex];
                    //separate by spaces
                    string[] comments = line.Split(' ');
                    //see if there's a Mods entry and parse it
                    foreach (string comment in comments)
                    {
                        string[] splitComment = comment.Split('=');
                        if (splitComment[0].Equals("Mods"))
                        {
                            string[] mods = splitComment[1].Split('/');
                            for (int i = mods.Length - 1; i > 0; i--)
                            {
                                string[] modInfo = mods[i].Split(',');
                                int index = Convert.ToInt32(modInfo[0]);
                                string mod = modInfo[2];
                                string metaMorpheusMod = PrositToMetaMorpheusModDictionary[mod];
                                //add the mod into the sequence
                                string leftSeq = sequence.Substring(0, index + 1);
                                string rightSeq = sequence.Substring(index + 1);
                                sequence = leftSeq + metaMorpheusMod + rightSeq;
                            }
                        }
                        else if (splitComment[0].Equals("iRT"))
                        {
                            rt = Convert.ToDouble(splitComment[1]);
                        }
                    }
                    //skip first line of fragment ions
                    lineIndex++;
                }
                else //we're going throught the fragment ions
                {
                    string[] fragmentInfo = line.Split("\t").Select(b => b.Trim()).ToArray();
                    var experMz = double.Parse(fragmentInfo[0], CultureInfo.InvariantCulture);
                    var experIntensity = double.Parse(fragmentInfo[1], CultureInfo.InvariantCulture);
                    string[] ionInfo = fragmentInfo[2].Split(new char[] { '/', '\"', ')', '(' }, StringSplitOptions.RemoveEmptyEntries).Select(b => b.Trim()).ToArray();
                    ionInfo[1] = ionInfo[1].Replace("ppm", "", ignoreCase: true, CultureInfo.InvariantCulture);

                    //TODO: figure out a more robust way to do this
                    string spectrumPeakProductType = ionInfo[0].ToCharArray()[0].ToString();

                    int fragmentNumber = int.Parse(new string(ionInfo[0].Split(new char[] { '^' })[0].Where(char.IsDigit).ToArray()));

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
            //add last entry
            if (lines.Length != 0)
            {
                LibrarySpectrum librarySpectrum = new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
                spectralLibraryDictionary.Add(librarySpectrum.Name, librarySpectrum);
            }

            return spectralLibraryDictionary;
        }
    }
}
