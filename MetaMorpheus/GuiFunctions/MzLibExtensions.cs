using System;
using EngineLayer;
using MassSpectrometry;
using Omics;
using Proteomics.ProteolyticDigestion;
using Readers;

namespace GuiFunctions
{
    /// <summary>
    /// Provides additional functionality for the MzLib library.
    /// </summary>
    public static class MzLibExtensions
    {

        /// <summary>
        /// Converts the given <see cref="DeconvolutionParameters"/> to a <see cref="DeconParamsViewModel"/>.
        /// </summary>
        /// <param name="parameters">The deconvolution parameters to convert.</param>
        /// <returns>A <see cref="DeconParamsViewModel"/> representing the given parameters.</returns>
        /// <exception cref="NotImplementedException">
        /// Thrown when the type of <paramref name="parameters"/> is not supported.
        /// </exception>
        public static DeconParamsViewModel ToViewModel(this DeconvolutionParameters parameters)
        {
            if (parameters is ClassicDeconvolutionParameters classicParams)
            {
                return new ClassicDeconParamsViewModel(classicParams);
            }
            else if (parameters is IsoDecDeconvolutionParameters isoParams)
            {
                return new IsoDecDeconParamsViewModel(isoParams);
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        public static bool IsCrossLinkedPeptide(this SpectrumMatchFromTsv sm)
        {
            return sm is PsmFromTsv { BetaPeptideBaseSequence: not null };
        }

        public static bool IsPeptide(this SpectrumMatchFromTsv sm)
        {
            //if (sm is OsmFromTsv)
            //    return false;
            return true;
        }

        public static IBioPolymerWithSetMods ToBioPolymerWithSetMods(this SpectrumMatchFromTsv sm, string fullSequence = null)
        {
            //if (sm.IsPeptide())
                return new PeptideWithSetModifications(fullSequence ?? sm.FullSequence, GlobalVariables.AllModsKnownDictionary);
            //else
            //    return new OligoWithSetMods(fullSequence ?? sm.FullSequence, GlobalVariables.AllRnaModsKnownDictionary);
        }

        public static SpectrumMatchFromTsv ReplaceFullSequence(this SpectrumMatchFromTsv sm, string fullSequence, string baseSequence = "")
        {
            //if (sm.IsPeptide())
                return new PsmFromTsv(sm as PsmFromTsv, fullSequence, baseSequence: baseSequence);
            //else
            //    return new OsmFromTsv(sm as OsmFromTsv, fullSequence, baseSequence: baseSequence);
        }
    }
}
