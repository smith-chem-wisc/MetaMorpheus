using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Localization
{
    /// <summary>
    /// The mass difference between an experimental precursor and the theoretical mass of the assigned peptide is determined. The LocalizationEngine attempts
    /// to localize this mass to one of the residues. It does this by adding the mass difference to each theoretical ion mass and looking for additional matches
    /// in the experimental spectrum. This engine should only be run for open, notch or custom searches. It should not be run for exact mass or missed
    /// monoisopic searches.
    /// </summary>
    public class LocalizationEngine : MetaMorpheusEngine
    {
        private readonly IEnumerable<PeptideSpectralMatch> AllResultingIdentifications;
        private readonly MsDataFile MyMsDataFile;

        public LocalizationEngine(IEnumerable<PeptideSpectralMatch> allResultingIdentifications, MsDataFile myMsDataFile, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            AllResultingIdentifications = allResultingIdentifications;
            MyMsDataFile = myMsDataFile;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // don't try to localize mass differences for ambiguous peptides
            PeptideSpectralMatch[] unambiguousPsms = AllResultingIdentifications.Where(b => b.FullSequence != null).ToArray();

            double psmsSearched = 0;
            int oldPercentProgress = 0;

            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Localizing mass-differences... ", NestedIds));

            Parallel.ForEach(Partitioner.Create(0, unambiguousPsms.Length),
                new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                (range, loopState) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        PeptideSpectralMatch psm = unambiguousPsms[i];

                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { break; }

                        MsDataScan scan = MyMsDataFile.GetOneBasedScan(psm.ScanNumber);
                        Ms2ScanWithSpecificMass scanWithSpecificMass = new Ms2ScanWithSpecificMass(scan, psm.ScanPrecursorMonoisotopicPeakMz, psm.ScanPrecursorCharge, psm.FullFilePath, CommonParameters);
                        PeptideWithSetModifications peptide = psm.BestMatchingPeptides.First().Peptide;
                        double massDifference = psm.ScanPrecursorMass - peptide.MonoisotopicMass;

                        // this section will iterate through all residues of the peptide and try to localize the mass-diff at each residue and report a score for each residue
                        var localizedScores = new List<double>();
                        for (int r = 0; r < peptide.Length; r++)
                        {
                            // create new PeptideWithSetMods with unidentified mass difference at the given residue
                            PeptideWithSetModifications peptideWithLocalizedMassDiff = peptide.Localize(r, massDifference);

                            // this is the list of theoretical products for this peptide with mass-difference on this residue
                            List<Product> productsWithLocalizedMassDiff = peptideWithLocalizedMassDiff.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus).ToList();

                            var matchedIons = MatchFragmentIons(scanWithSpecificMass, productsWithLocalizedMassDiff, CommonParameters);

                            // score when the mass-diff is on this residue
                            double localizedScore = CalculatePeptideScore(scan, matchedIons);

                            localizedScores.Add(localizedScore);
                        }

                        psm.LocalizedScores = localizedScores;

                        // report search progress (PSM mass differences localized so far out of total)
                        psmsSearched++;
                        var percentProgress = (int)((psmsSearched / unambiguousPsms.Length) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Localizing mass-differences... ", NestedIds));
                        }
                    }
                });

            return new LocalizationEngineResults(this);
        }
    }
}