using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class SpectrumWithIndexAndNotchInfo
    {
        public LibrarySpectrum TheSpectrum;
        public int Notch;
        public int SpectrumIndex;

        public SpectrumWithIndexAndNotchInfo(LibrarySpectrum theSpectrum, int notch, int spectrumIndex)
        {
            TheSpectrum = theSpectrum;
            Notch = notch;
            SpectrumIndex = spectrumIndex;
        }
    }
}
