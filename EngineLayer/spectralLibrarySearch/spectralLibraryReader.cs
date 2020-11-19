using Easy.Common.Extensions;
using Proteomics.Fragmentation;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.spectralLibrarySearch
{
    public class SpectralLibraryReader
    {
        public Dictionary<String, Spectrum> SpectralLibraryDictionary { get; set; }
        public List<Spectrum> targetSpectralLibrary { get; set; }
        public List<Spectrum> decoySpectralLibrary { get; set; }

        //int i = 1;

        public SpectralLibraryReader(string filePath)
        {
            targetSpectralLibrary = new List<Spectrum>();
            decoySpectralLibrary = new List<Spectrum>();
            SpectralLibraryDictionary = new Dictionary<string, Spectrum>();
            string[] lines;

            try
            {
                lines = File.ReadAllLines(filePath);
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Could not read file: " + e.Message);
            }

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
           
            //load separate spectrums 
            for (int i = 0; i < nameLine.Count - 1; i++)
                {
                    //load each spectrum
                    var singleSpectrum = new Spectrum();
                    for (int j = nameLine[i]; j < nameLine[i + 1]-1; j++)
                    {
                        //get name of each spectrum
                        if (lines[j].Contains("name", StringComparison.OrdinalIgnoreCase))
                        {
                            try
                            {
                                string[] name = lines[j].Split(new char[] { ':', '=' }, 2).Select(b => b.Trim()).ToArray();
                                singleSpectrum.SequenceWithCharge = name[1];
                                string[] sequenceAndCharge = singleSpectrum.SequenceWithCharge.Split(new char[] { '/'}, 2).Select(b => b.Trim()).ToArray();
                                singleSpectrum.Sequence = sequenceAndCharge[0];
                                if (sequenceAndCharge.Length > 1)
                                {
                                    singleSpectrum.Charge_state = Int32.Parse(sequenceAndCharge[1]);
                                }
                                else
                                {
                                    singleSpectrum.Charge_state = 1;
                                }
                               
                            }
                            catch (Exception e)
                            {
                                throw new MetaMorpheusException("Could not find the name : " + e.Message);
                            }
                        }

                        //get MW of each spectrum
                        else if ((lines[j].Contains("MW", StringComparison.OrdinalIgnoreCase) || lines[j].Contains("Monoisotopic Mass", StringComparison.OrdinalIgnoreCase))&& !lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                        {
                        string[] mw = lines[j].Split(":", 2).Select(b => b.Trim()).ToArray();
                            try
                            {
                                var x = mw[0] + mw[1];
                                singleSpectrum.MW = double.Parse(mw[1]);
                            }
                            catch (Exception e)
                            {
                                throw new MetaMorpheusException("Could not find the MW : " + e.Message);
                            }
                        }

                        // get information from comment
                        if (lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                        {

                            string[] comment = lines[j].Split(" ").Select(b => b.Trim()).ToArray();
                            for (int l = 0; l < comment.Length; l++)
                            {
                                if (comment[l].Contains("parent", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("precursor", StringComparison.OrdinalIgnoreCase))
                                {
                                    try
                                    {
                                        string[] precursorMz = comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray();
                                        singleSpectrum.PrecursorMz = double.Parse(precursorMz[1]);
                                    }
                                    catch (Exception e)
                                    {
                                        throw new MetaMorpheusException("Could not find the mz of precursor : " + e.Message);
                                    }
                                }

                                if (comment[l].Contains("iRT", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("retention time", StringComparison.OrdinalIgnoreCase))
                                {
                                    
                                    try
                                    {
                                        string[] rententionTime = comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray();
                                        singleSpectrum.rententionTime = double.Parse(rententionTime[1]);
                                    }
                                    catch 
                                    {
                                    }
                                }

                            }
                        }

                        if (lines[j].Contains("peaks", StringComparison.OrdinalIgnoreCase))
                        {
                        var numberOfPeaks = 0;
                        string[] numPeaks = lines[j].Split(":").Select(b => b.Trim()).ToArray();
                        try
                        {
                            numberOfPeaks = int.Parse(numPeaks[1]);
                        }
                        catch
                        {
                            
                        }
                          //load each peak 
                            var peaks = new List<MatchedFragmentIon>();
                           
                            for (int k = j + 1; k < j + 1 + numberOfPeaks; k++)
                            {
                                string[] eachPeak = lines[k].Split("\t").Select(b => b.Trim()).ToArray();
                                var experMz = double.Parse(eachPeak[0]);
                                var experIntensity = double.Parse(eachPeak[1]);
                                string[] ions = eachPeak[2].Split(new char[] { '/', '\"', 'p' }, StringSplitOptions.RemoveEmptyEntries).Select(b => b.Trim()).ToArray();
                                var spectrumPeakProductType = ions[0].ToCharArray()[0].ToString();
                                var fragmentNumber = (int)Char.GetNumericValue(ions[0].ToCharArray()[1]);
                                int ionCharge;
                                if (ions[0].ToCharArray().Length>3)
                               {
                                    ionCharge = (int)Char.GetNumericValue(ions[0].ToCharArray()[3]);
                               }
                               else
                               {
                                    ionCharge = 1;
                                }
            
                                ProductType peakProductType = (ProductType)Enum.Parse(typeof(ProductType), spectrumPeakProductType, true);
                                FragmentationTerminus terminus = (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus), "None", true);
                                ProductType.a.ToString();
                                var product = new Product(peakProductType, terminus, experMz, fragmentNumber, 0, 0);
                                peaks.Add(new MatchedFragmentIon(ref product, experMz, experIntensity, ionCharge));
                            }
                            singleSpectrum.MatchedFragmentIons = peaks;
                        }
                    }
                    if (singleSpectrum.SequenceWithCharge != null &&  singleSpectrum.MatchedFragmentIons.Count != 0)
                    {
                        targetSpectralLibrary.Add( singleSpectrum);
                        SpectralLibraryDictionary.Add(singleSpectrum.SequenceWithCharge, singleSpectrum);
                }
                }
            if (targetSpectralLibrary.Count == 0)
            {
                Console.WriteLine("the library doesn't contain any spectrums!");
            }
           
        }

        public List<Spectrum> SpectralLibrary { get; set; }

        public List<Spectrum> SortSpectralLibrary(List<Spectrum> unsortedSpectralLibrary)
        {
            List<Spectrum> SortedSpectralLibrary = unsortedSpectralLibrary.OrderByDescending(o => o.PrecursorMz).ToList();
            return SortedSpectralLibrary;
        }

        

    }
}
