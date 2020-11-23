using Chemistry;
using EngineLayer.FdrAnalysis;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class PeptideSpectralMatch
    {
        public const double ToleranceForScoreDifferentiation = 1e-9;
        private List<(int Notch, PeptideWithSetModifications Pwsm)> _BestMatchingPeptides;

        public PeptideSpectralMatch(PeptideWithSetModifications peptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters, List<MatchedFragmentIon> matchedFragmentIons, double xcorr = 0)
        {
            _BestMatchingPeptides = new List<(int, PeptideWithSetModifications)>();
            ScanIndex = scanIndex;
            FullFilePath = scan.FullFilePath;
            ScanNumber = scan.OneBasedScanNumber;
            PrecursorScanNumber = scan.OneBasedPrecursorScanNumber;
            ScanRetentionTime = scan.RetentionTime;
            ScanExperimentalPeaks = scan.NumPeaks;
            TotalIonCurrent = scan.TotalIonCurrent;
            ScanPrecursorCharge = scan.PrecursorCharge;
            ScanPrecursorMonoisotopicPeakMz = scan.PrecursorMonoisotopicPeakMz;
            ScanPrecursorMass = scan.PrecursorMass;
            DigestionParams = commonParameters.DigestionParams;
            PeptidesToMatchingFragments = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();
            Xcorr = xcorr;
            NativeId = scan.NativeId;
            RunnerUpScore = commonParameters.ScoreCutoff;

            AddOrReplace(peptide, score, notch, true, matchedFragmentIons, xcorr);
        }

        public ChemicalFormula ModsChemicalFormula { get; private set; } // these fields will be null if they are ambiguous
        public string FullSequence { get; private set; }
        public string EssentialSequence { get; private set; }
        public int? Notch { get; private set; }
        public string BaseSequence { get; private set; }
        public int? PeptideLength { get; private set; }
        public int? OneBasedStartResidueInProtein { get; private set; }
        public int? OneBasedEndResidueInProtein { get; private set; }
        public double? PeptideMonisotopicMass { get; private set; }
        public int? ProteinLength { get; private set; }
        public string ProteinAccession { get; private set; }
        public string Organism { get; private set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; protected set; }
        public int PsmCount { get; internal set; }
        public Dictionary<string, int> ModsIdentified { get; private set; } // these should never be null under normal circumstances
        public List<double> LocalizedScores { get; internal set; }
        public int ScanNumber { get; }
        public int? PrecursorScanNumber { get; }
        public double ScanRetentionTime { get; }
        public int ScanExperimentalPeaks { get; }
        public double TotalIonCurrent { get; }
        public int ScanPrecursorCharge { get; }
        public double ScanPrecursorMonoisotopicPeakMz { get; }
        public double ScanPrecursorMass { get; }
        public string FullFilePath { get; private set; }
        public int ScanIndex { get; }
        public int NumDifferentMatchingPeptides { get { return _BestMatchingPeptides.Count; } }
        public FdrInfo FdrInfo { get; private set; }
        public PsmData PsmData_forPEPandPercolator { get; set; }

        public double Score { get; private set; }
        public double Xcorr;
        public double SpectralAngle { get; set; }
        public string NativeId; // this is a property of the scan. used for mzID writing

        public double DeltaScore { get { return (Score - RunnerUpScore); } }

        public double RunnerUpScore { get; set; }
        public bool IsDecoy { get; private set; }
        public bool IsContaminant { get; private set; }

        public DigestionParams DigestionParams { get; }
        public Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> PeptidesToMatchingFragments { get; private set; }
        
        public IEnumerable<(int Notch, PeptideWithSetModifications Peptide)> BestMatchingPeptides
        {
            get
            {
                return _BestMatchingPeptides.OrderBy(p => p.Pwsm.FullSequence)
                    .ThenBy(p => p.Pwsm.Protein.Accession)
                    .ThenBy(p => p.Pwsm.OneBasedStartResidueInProtein);
            }
        }

        public static string GetTabSeparatedHeader()
        {
            return string.Join("\t", DataDictionary(null, null).Keys);
        }

        public void AddOrReplace(PeptideWithSetModifications pwsm, double newScore, int notch, bool reportAllAmbiguity, List<MatchedFragmentIon> matchedFragmentIons, double newXcorr)
        {
            if (newScore - Score > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                _BestMatchingPeptides.Clear();
                _BestMatchingPeptides.Add((notch, pwsm));

                if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
                {
                    RunnerUpScore = Score;
                }

                Score = newScore;
                Xcorr = newXcorr;

                PeptidesToMatchingFragments.Clear();
                PeptidesToMatchingFragments.Add(pwsm, matchedFragmentIons);
            }
            else if (newScore - Score > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                _BestMatchingPeptides.Add((notch, pwsm));

                if (!PeptidesToMatchingFragments.ContainsKey(pwsm))
                {
                    PeptidesToMatchingFragments.Add(pwsm, matchedFragmentIons);
                }
            }
            else if (newScore - RunnerUpScore > ToleranceForScoreDifferentiation)
            {
                RunnerUpScore = newScore;
            }
        }

        //PEP-Value analysis identifies ambiguous peptides with lower probability. These are removed from the bestmatchingpeptides dictionary, which lowers ambiguity.
        public void RemoveThisAmbiguousPeptide(int notch, PeptideWithSetModifications pwsm)
        {
            _BestMatchingPeptides.Remove((notch, pwsm));
            this.ResolveAllAmbiguities();
        }

        public override string ToString()
        {
            return ToString(new Dictionary<string, int>());
        }

        public string ToString(IReadOnlyDictionary<string, int> ModstoWritePruned)
        {
            return string.Join("\t", DataDictionary(this, ModstoWritePruned).Values);
        }

        public static Dictionary<string, string> DataDictionary(PeptideSpectralMatch psm, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddBasicMatchData(s, psm);
            PsmTsvWriter.AddPeptideSequenceData(s, psm, ModsToWritePruned);
            PsmTsvWriter.AddMatchedIonsData(s, psm?.MatchedFragmentIons);
            PsmTsvWriter.AddMatchScoreData(s, psm);
            return s;
        }

        public void SetFdrValues(double cumulativeTarget, double cumulativeDecoy, double qValue, double cumulativeTargetNotch, double cumulativeDecoyNotch, double qValueNotch, double pep, double pepQValue)
        {
            FdrInfo = new FdrInfo
            {
                CumulativeTarget = cumulativeTarget,
                CumulativeDecoy = cumulativeDecoy,
                QValue = qValue,
                CumulativeTargetNotch = cumulativeTargetNotch,
                CumulativeDecoyNotch = cumulativeDecoyNotch,
                QValueNotch = qValueNotch,
                PEP = pep,
                PEP_QValue = pepQValue
            };
        }

        /// <summary>
        /// This method saves properties of this PSM for internal use. It is NOT used for any output.
        /// These resolved fields are (usually) null if there is more than one option.
        /// e.g., if this PSM can be explained by more than one base sequence, the BaseSequence property will be null
        /// </summary>
        public void ResolveAllAmbiguities()
        {
            IsDecoy = _BestMatchingPeptides.Any(p => p.Pwsm.Protein.IsDecoy);
            IsContaminant = _BestMatchingPeptides.Any(p => p.Pwsm.Protein.IsContaminant);
            FullSequence = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.FullSequence)).ResolvedValue;
            BaseSequence = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.BaseSequence)).ResolvedValue;
            PeptideLength = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.Length)).ResolvedValue;
            OneBasedStartResidueInProtein = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.OneBasedStartResidueInProtein)).ResolvedValue;
            OneBasedEndResidueInProtein = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.OneBasedEndResidueInProtein)).ResolvedValue;
            ProteinLength = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.Protein.Length)).ResolvedValue;
            PeptideMonisotopicMass = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.MonoisotopicMass)).ResolvedValue;
            ProteinAccession = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.Protein.Accession)).ResolvedValue;
            Organism = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.Protein.Organism)).ResolvedValue;
            ModsIdentified = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.AllModsOneIsNterminus)).ResolvedValue;
            ModsChemicalFormula = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.AllModsOneIsNterminus.Select(c => (c.Value)))).ResolvedValue;
            Notch = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Notch)).ResolvedValue;

            // if the PSM matches a target and a decoy and they are the SAME SEQUENCE, remove the decoy
            if (IsDecoy)
            {
                bool removedPeptides = false;
                var hits = _BestMatchingPeptides.GroupBy(p => p.Pwsm.FullSequence);

                foreach (var hit in hits)
                {
                    if (hit.Any(p => p.Pwsm.Protein.IsDecoy) && hit.Any(p => !p.Pwsm.Protein.IsDecoy))
                    {
                        // at least one peptide with this sequence is a target and at least one is a decoy
                        // remove the decoys with this sequence
                        var pwsmToRemove = _BestMatchingPeptides.Where(p => p.Pwsm.FullSequence == hit.Key && p.Pwsm.Protein.IsDecoy).ToList();
                        _BestMatchingPeptides.RemoveAll(p => p.Pwsm.FullSequence == hit.Key && p.Pwsm.Protein.IsDecoy);
                        foreach ((int, PeptideWithSetModifications) pwsm in pwsmToRemove)
                        {
                            PeptidesToMatchingFragments.Remove(pwsm.Item2);
                        }

                        removedPeptides = true;
                    }
                }

                if (removedPeptides)
                {
                    ResolveAllAmbiguities();
                }
            }

            // TODO: technically, different peptide options for this PSM can have different matched ions
            // we can write a Resolve method for this if we want...
            MatchedFragmentIons = PeptidesToMatchingFragments.First().Value;
        }

        public static int GetLongestIonSeriesBidirectional(Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> PeptidesToMatchingFragments, PeptideWithSetModifications peptide)
        {
            List<int> maxDiffs = new List<int> { 1 };
            if (PeptidesToMatchingFragments != null && PeptidesToMatchingFragments.TryGetValue(peptide, out var matchedFragments) && matchedFragments != null && matchedFragments.Any())
            {
                var jointSeries = matchedFragments.Select(p => p.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

                if (jointSeries.Count > 0)
                {
                    jointSeries.Sort();

                    List<int> aminoAcidPostionsThatCouldBeObserved = Enumerable.Range(1, peptide.BaseSequence.Length).ToList();

                    List<int> missing = aminoAcidPostionsThatCouldBeObserved.Except(jointSeries).ToList();

                    int localMaxDiff = 0;
                    for (int i = 0; i < aminoAcidPostionsThatCouldBeObserved.Count; i++)
                    {
                        if (!missing.Contains(aminoAcidPostionsThatCouldBeObserved[i]))
                        {
                            localMaxDiff++;
                        }
                        else
                        {
                            maxDiffs.Add(localMaxDiff);
                            localMaxDiff = 0;
                        }
                    }
                    maxDiffs.Add(localMaxDiff);
                }
            }

            return maxDiffs.Max();
        }

        public static int GetCountComplementaryIons(Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> PeptidesToMatchingFragments, PeptideWithSetModifications peptide)
        {
            if (PeptidesToMatchingFragments != null && PeptidesToMatchingFragments.TryGetValue(peptide, out var matchedFragments) && matchedFragments != null && matchedFragments.Any())
            {
                List<int> nIons = matchedFragments.Where(f => f.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N).Select(f => f.NeutralTheoreticalProduct.FragmentNumber).ToList();
                List<int> cIons = matchedFragments.Where(f => f.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C).Select(f => (peptide.BaseSequence.Length - f.NeutralTheoreticalProduct.FragmentNumber)).ToList();
                if (nIons.Any() && cIons.Any())
                {
                    return nIons.Intersect(cIons).Count();
                }
                else
                {
                    return 0;
                }
            }
            else
            {
                return 0;
            }
        }

        /// <summary>
        /// This method changes the base and full sequences to reflect heavy silac labels
        /// translates SILAC sequence into the proper peptide sequence ("PEPTIDEa" into "PEPTIDEK(+8.014)")
        /// </summary>
        public void ResolveHeavySilacLabel(List<SilacLabel> labels, IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            //FullSequence
            FullSequence = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.FullSequence)).ResolvedString; //string, not value
            FullSequence = SilacConversions.GetAmbiguousLightSequence(FullSequence, labels, false);

            //BaseSequence
            BaseSequence = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.BaseSequence)).ResolvedString; //string, not value
            BaseSequence = SilacConversions.GetAmbiguousLightSequence(BaseSequence, labels, true);

            //EssentialSequence
            EssentialSequence = PsmTsvWriter.Resolve(_BestMatchingPeptides.Select(b => b.Pwsm.EssentialSequence(modsToWritePruned))).ResolvedString; //string, not value
            EssentialSequence = SilacConversions.GetAmbiguousLightSequence(EssentialSequence, labels, false);
        }

        /// <summary>
        /// This method is used by protein parsimony to remove PeptideWithSetModifications objects that have non-parsimonious protein associations
        /// </summary>
        public void TrimProteinMatches(HashSet<Protein> parsimoniousProteins)
        {
            if (IsDecoy)
            {
                if (_BestMatchingPeptides.Any(p => parsimoniousProteins.Contains(p.Pwsm.Protein) && p.Pwsm.Protein.IsDecoy))
                {
                    _BestMatchingPeptides.RemoveAll(p => !parsimoniousProteins.Contains(p.Pwsm.Protein));
                }
                // else do nothing
            }
            else
            {
                _BestMatchingPeptides.RemoveAll(p => !parsimoniousProteins.Contains(p.Pwsm.Protein));
            }

            ResolveAllAmbiguities();
        }

        /// <summary>
        /// This method is used by protein parsimony to add PeptideWithSetModifications objects for modification-agnostic parsimony
        /// </summary>
        public void AddProteinMatch((int, PeptideWithSetModifications) peptideWithNotch, List<MatchedFragmentIon> mfi)
        {
            if (!_BestMatchingPeptides.Select(p => p.Pwsm).Contains(peptideWithNotch.Item2))
            {
                _BestMatchingPeptides.Add(peptideWithNotch);
                if (!PeptidesToMatchingFragments.ContainsKey(peptideWithNotch.Item2))
                {
                    PeptidesToMatchingFragments.Add(peptideWithNotch.Item2, mfi);
                }
                ResolveAllAmbiguities();
            }
        }

        /// <summary>
        /// This method is used by SILAC quantification to add heavy/light psms
        /// Don't have access to the scans at that point, so a new contructor is needed
        /// </summary>
        public PeptideSpectralMatch Clone(List<(int Notch, PeptideWithSetModifications Peptide)> bestMatchingPeptides)
        {
            return new PeptideSpectralMatch(this, bestMatchingPeptides);
        }

        private PeptideSpectralMatch(PeptideSpectralMatch psm, List<(int Notch, PeptideWithSetModifications Peptide)> bestMatchingPeptides)
        {
            _BestMatchingPeptides = bestMatchingPeptides;
            BaseSequence = PsmTsvWriter.Resolve(bestMatchingPeptides.Select(b => b.Peptide.BaseSequence)).ResolvedValue;
            FullSequence = PsmTsvWriter.Resolve(bestMatchingPeptides.Select(b => b.Peptide.FullSequence)).ResolvedValue;

            ModsChemicalFormula = psm.ModsChemicalFormula;
            Notch = psm.Notch;
            PeptideLength = psm.PeptideLength;
            OneBasedStartResidueInProtein = psm.OneBasedStartResidueInProtein;
            OneBasedEndResidueInProtein = psm.OneBasedEndResidueInProtein;
            PeptideMonisotopicMass = psm.PeptideMonisotopicMass;
            ProteinLength = psm.ProteinLength;
            ProteinAccession = psm.ProteinAccession;
            Organism = psm.Organism;
            MatchedFragmentIons = psm.MatchedFragmentIons;
            PsmCount = psm.PsmCount;
            ModsIdentified = psm.ModsIdentified;
            LocalizedScores = psm.LocalizedScores;
            ScanNumber = psm.ScanNumber;
            PrecursorScanNumber = psm.PrecursorScanNumber;
            ScanRetentionTime = psm.ScanRetentionTime;
            ScanExperimentalPeaks = psm.ScanExperimentalPeaks;
            TotalIonCurrent = psm.TotalIonCurrent;
            ScanPrecursorCharge = psm.ScanPrecursorCharge;
            ScanPrecursorMonoisotopicPeakMz = psm.ScanPrecursorMonoisotopicPeakMz;
            ScanPrecursorMass = psm.ScanPrecursorMass;
            FullFilePath = psm.FullFilePath;
            ScanIndex = psm.ScanIndex;
            FdrInfo = psm.FdrInfo;
            Score = psm.Score;
            Xcorr = psm.Xcorr;
            RunnerUpScore = psm.RunnerUpScore;
            IsDecoy = psm.IsDecoy;
            IsContaminant = psm.IsContaminant;
            DigestionParams = psm.DigestionParams;
            PeptidesToMatchingFragments = psm.PeptidesToMatchingFragments;
            SpectralAngle = psm.SpectralAngle;
        }
    }
}