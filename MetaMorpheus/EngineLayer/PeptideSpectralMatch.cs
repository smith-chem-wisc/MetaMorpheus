using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Transcriptomics;

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

        public new DigestionParams DigestionParams => base.DigestionParams as DigestionParams;

        #region Silac

        /// <summary>
        /// This method changes the base and full sequences to reflect heavy silac labels
        /// translates SILAC sequence into the proper peptide sequence ("PEPTIDEa" into "PEPTIDEK(+8.014)")
        /// </summary>
        public void ResolveHeavySilacLabel(List<SilacLabel> labels, IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            //FullSequence
            FullSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.FullSequence)).ResolvedString; //string, not value
            FullSequence = SilacConversions.GetAmbiguousLightSequence(FullSequence, labels, false);

            //BaseSequence
            BaseSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.BaseSequence)).ResolvedString; //string, not value
            BaseSequence = SilacConversions.GetAmbiguousLightSequence(BaseSequence, labels, true);

            //EssentialSequence
            EssentialSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.EssentialSequence(modsToWritePruned))).ResolvedString; //string, not value
            EssentialSequence = SilacConversions.GetAmbiguousLightSequence(EssentialSequence, labels, false);
        }

        /// <summary>
        /// This method is used by SILAC quantification to add heavy/light psms
        /// Don't have access to the scans at that point, so a new contructor is needed
        /// </summary>
        public PeptideSpectralMatch Clone(List<(int Notch, IBioPolymerWithSetMods Peptide)> bestMatchingPeptides) => new PeptideSpectralMatch(this, bestMatchingPeptides);
        
        protected PeptideSpectralMatch(SpectralMatch psm, List<(int Notch, IBioPolymerWithSetMods Peptide)> bestMatchingPeptides) 
            : base(psm, bestMatchingPeptides)
        {
        }

        #endregion

    }

    public class OligoSpectralMatch : SpectralMatch
    {
        public OligoSpectralMatch(IBioPolymerWithSetMods peptide, int notch, double score, int scanIndex,
            Ms2ScanWithSpecificMass scan, CommonParameters commonParameters,
            List<MatchedFragmentIon> matchedFragmentIons, double xcorr = 0) : base(peptide, notch, score, scanIndex,
            scan, commonParameters, matchedFragmentIons, xcorr)
        {

        }

        protected OligoSpectralMatch(SpectralMatch psm, List<(int Notch, IBioPolymerWithSetMods Peptide)> bestMatchingPeptides)
            : base(psm, bestMatchingPeptides)
        {
        }
    }
}
