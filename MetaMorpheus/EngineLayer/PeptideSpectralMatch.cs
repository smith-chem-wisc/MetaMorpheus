using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;

namespace EngineLayer
{
    public class PeptideSpectralMatch : IComparable<PeptideSpectralMatch>
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


        /// <summary>
        /// There are a few key locations in MetaMorpheus where we want to have psms sorted in a consistent manner.
        /// These are for q-value determination and for when we write the psms to psmtsv. 
        /// </summary>
        /// <param name="otherPsm"></param>
        /// <returns></returns>
        public int CompareTo(PeptideSpectralMatch otherPsm)
        {
            if (Math.Abs(this.Score - otherPsm.Score) > ToleranceForScoreDifferentiation)
            {
                return this.Score.CompareTo(otherPsm.Score);
            }
            else if (Math.Abs(this.DeltaScore - otherPsm.DeltaScore) > ToleranceForScoreDifferentiation)
            {
                return this.RunnerUpScore.CompareTo(otherPsm.RunnerUpScore);
            }
            else if (otherPsm.PrecursorMassErrorPpm != null && (Math.Abs(otherPsm.PrecursorMassErrorPpm.First() - this.PrecursorMassErrorPpm.First()) > 0.01))
            {
                return Math.Abs(otherPsm.PrecursorMassErrorPpm.First()).CompareTo(Math.Abs(this.PrecursorMassErrorPpm.First())); //precursor mass errors defined for both otherPsms. Reverse the comparision so that lower ppm error comes first
            }
            return otherPsm.ScanNumber.CompareTo(this.ScanNumber); //reverse the comparision so that the lower scan number comes first.
        }
    }
}
