using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using Transcriptomics;
using Transcriptomics.Digestion;

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
            if (sm is OsmFromTsv)
                return false;
            return true;
        }

        public static string GetDigestionProductLabel(this SpectrumMatchFromTsv sm)
        {
            if (sm.IsPeptide())
                return "Peptide";
            else
                return "Oligo";
        }

        public static IBioPolymerWithSetMods ToBioPolymerWithSetMods(this SpectrumMatchFromTsv sm, string fullSequence = null)
        {
            int startResidue = 0, endResidue = sm.BaseSeq.Length;

            if (sm.StartAndEndResiduesInParentSequence != null)
            {
                var splitStartAndEnd = sm.StartAndEndResiduesInParentSequence
                    .Replace('[', ' ')
                    .Replace(']', ' ').Split("to");
                startResidue = int.Parse(splitStartAndEnd[0].Trim());
                endResidue = int.Parse(splitStartAndEnd[1].Split('|')[0].Trim());
            }

            if (sm.IsPeptide())
                return new PeptideWithSetModifications(fullSequence ?? sm.FullSequence, GlobalVariables.AllModsKnownDictionary, oneBasedStartResidueInProtein: startResidue, oneBasedEndResidueInProtein: endResidue);
            else
            {
                // Determine termini based on position in parent sequence, future versions would read it from the SM. 
                IHasChemicalFormula fivePrimeTerminus = startResidue == 1 ? NucleicAcid.DefaultFivePrimeTerminus : Rnase.DefaultFivePrimeTerminus;
                IHasChemicalFormula threePrimeTerminus = endResidue == sm.BaseSeq.Length ? NucleicAcid.DefaultThreePrimeTerminus : Rnase.DefaultThreePrimeTerminus;

                return new OligoWithSetMods(fullSequence ?? sm.FullSequence, GlobalVariables.AllRnaModsKnownDictionary,
                    oneBaseStartResidue: startResidue, oneBasedEndResidue: endResidue, fivePrimeTerminus: fivePrimeTerminus, threePrimeTerminus: threePrimeTerminus);
            }
        }

        public static SpectrumMatchFromTsv ReplaceFullSequence(this SpectrumMatchFromTsv sm, string fullSequence, string baseSequence = "")
        {
            if (sm.IsPeptide())
                return new PsmFromTsv(sm as PsmFromTsv, fullSequence, baseSequence: baseSequence);
            else
                return new OsmFromTsv(sm as OsmFromTsv, fullSequence, baseSequence: baseSequence);
        }

        public static Dictionary<DissociationType, List<ProductType>> ProductsFromDissociationType(this SpectrumMatchFromTsv sm)
        {
            if (sm.IsPeptide())
                return Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType;
            else
                return Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType;
        }

        public static Dictionary<string, Modification> AllModsKnownModificationDictionary(this SpectrumMatchFromTsv sm)
        {
            if (sm.IsPeptide())
                return GlobalVariables.AllModsKnownDictionary;
            else
                return GlobalVariables.AllRnaModsKnownDictionary;
        }

        /// <summary>
        /// Determines if a majority of values are within a range
        /// </summary>
        /// <param name="range"></param>
        /// <param name="values"></param>
        /// <returns></returns>
        public static bool MajorityWithin(this MzRange range, IEnumerable<double> values)
        {
            int within = values.Count(p => p >= range.Minimum && p <= range.Maximum);
            return within > values.Count() / 2;
        }

        public static IEnumerable<(int Start, int End)> GetStartAndEndPosition(this SpectrumMatchFromTsv sm)
        {
            foreach (var ambigSplit in sm.StartAndEndResiduesInParentSequence.Split('|'))
            {
                var split = ambigSplit.Replace("[", "").Replace("]", "").Split("to");
                yield return (int.Parse(split[0]), int.Parse(split[1]));
            }
        }
    }
}
