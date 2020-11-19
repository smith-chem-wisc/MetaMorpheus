using System;
using System.Collections.Generic;
using System.Text;
using Proteomics.Fragmentation;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.spectralLibrarySearch
{
    public class DecoyLibrary
    {

        public List<Spectrum> GenerateDecoyLibrary(List<Spectrum> unsortedSpectralLibrary, LibraryDecoyType decoyType)
        {
            var sortedSpectralLibrary = SortSpectralLibrary(unsortedSpectralLibrary);
            if (decoyType == LibraryDecoyType.None)
            {
                return new List<Spectrum>();
            }
            else if (decoyType == LibraryDecoyType.Reverse)
            {
                return PrecursorSwapDecoyLibraryGeneration(sortedSpectralLibrary);
            }
            else
            {
                throw new ArgumentException("Decoy type " + decoyType.ToString() + " is not implemented.");
            }
        }
        public List<Spectrum> SpectralLibrary { get; set; }

        public List<Spectrum> SortSpectralLibrary(List<Spectrum> unsortedSpectralLibrary)
        {
            List<Spectrum> SortedSpectralLibrary = unsortedSpectralLibrary.OrderByDescending(o => o.PrecursorMz).ToList();
            return SortedSpectralLibrary;
        }

        public List<Spectrum> PrecursorSwapDecoyLibraryGeneration(List<Spectrum> sortedSpectralLibrary)
        {
            var decoySpectralLibrary = new List<Spectrum>();
            double libraryHalfLenth = sortedSpectralLibrary.Count / 2;
            int libraryRoundUpHalfLenth = (int)Math.Ceiling(libraryHalfLenth);
            int libraryRoundDownHalfLenth = (int)libraryHalfLenth;
            for (int i=0; i< libraryRoundUpHalfLenth; i++)
            {
                int indexOfSwap = i + libraryRoundDownHalfLenth;
                while(indexOfSwap< sortedSpectralLibrary.Count)
                {
                    if (Math.Abs(sortedSpectralLibrary[i].PrecursorMz - sortedSpectralLibrary[indexOfSwap].PrecursorMz) < 8)
                    {
                        decoySpectralLibrary.Add(new Spectrum(sortedSpectralLibrary[i].SequenceWithCharge, sortedSpectralLibrary[i].Sequence, sortedSpectralLibrary[indexOfSwap].PrecursorMz, sortedSpectralLibrary[indexOfSwap].Charge_state, sortedSpectralLibrary[i].MatchedFragmentIons,true));
                        break;
                    }
                    else
                    {
                        indexOfSwap = indexOfSwap + 1;
                    }
                }
                
            }
            return decoySpectralLibrary;
        }



    }
}
