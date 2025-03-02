using System.Collections.Generic;
using System.Linq;
using EngineLayer.SpectrumMatch;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;

namespace EngineLayer
{
    public class PeptideSpectralMatch : SpectralMatch
    {
        public PeptideSpectralMatch(IBioPolymerWithSetMods peptide, int notch, double score, int scanIndex,
            Ms2ScanWithSpecificMass scan, CommonParameters commonParameters,
            List<MatchedFragmentIon> matchedFragmentIons, double xcorr = 0) : base(peptide, notch, score, scanIndex,
            scan, commonParameters, matchedFragmentIons, xcorr)
        {

        }
        #region Silac
            

        /// <summary>
        /// This method changes the base and full sequences to reflect heavy silac labels
        /// translates SILAC sequence into the proper peptide sequence ("PEPTIDEa" into "PEPTIDEK(+8.014)")
        /// </summary>
        public void ResolveHeavySilacLabel(List<SilacLabel> labels, IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            //FullSequence
            FullSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.WithSetMods.FullSequence)).ResolvedString; //string, not value
            FullSequence = SilacConversions.GetAmbiguousLightSequence(FullSequence, labels, false);

            //BaseSequence
            BaseSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.WithSetMods.BaseSequence)).ResolvedString; //string, not value
            BaseSequence = SilacConversions.GetAmbiguousLightSequence(BaseSequence, labels, true);

            //EssentialSequence
            EssentialSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.WithSetMods.EssentialSequence(modsToWritePruned))).ResolvedString; //string, not value
            EssentialSequence = SilacConversions.GetAmbiguousLightSequence(EssentialSequence, labels, false);
        }

        /// <summary>
        /// This method is used by SILAC quantification to add heavy/light psms
        /// Don't have access to the scans at that point, so a new contructor is needed
        /// </summary>
        public PeptideSpectralMatch Clone(List<TentativeSpectralMatch> bestMatchingPeptides) => new PeptideSpectralMatch(this, bestMatchingPeptides);
        
        protected PeptideSpectralMatch(SpectralMatch psm, List<TentativeSpectralMatch> bestMatchingPeptides) 
            : base(psm, bestMatchingPeptides)
        {
        }
        #endregion


    }
}
