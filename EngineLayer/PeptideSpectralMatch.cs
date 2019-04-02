using Chemistry;
using EngineLayer.FdrAnalysis;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class PeptideSpectralMatch
    {
        public const double ToleranceForScoreDifferentiation = 1e-9;
        private List<(int Notch, PeptideWithSetModifications Pwsm)> _BestMatchingPeptides;

        public PeptideSpectralMatch(PeptideWithSetModifications peptide, int notch, double score, int scanIndex, IScan scan, DigestionParams digestionParams, List<MatchedFragmentIon> matchedFragmentIons, double xcorr = 0)
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
            AllScores = new List<double>();
            DigestionParams = digestionParams;
            PeptidesToMatchingFragments = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();

            Xcorr = xcorr;
            

            AddOrReplace(peptide, score, notch, true, matchedFragmentIons, xcorr);
        }

        public ChemicalFormula ModsChemicalFormula { get; private set; } // these fields will be null if they are ambiguous
        public string FullSequence { get; private set; }
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
        public string FullFilePath { get; }
        public int ScanIndex { get; }
        public int NumDifferentMatchingPeptides { get { return _BestMatchingPeptides.Count; } }
        public FdrInfo FdrInfo { get; private set; }
        public double Score { get; private set; }
        public double Xcorr;
        public double DeltaScore { get; private set; }
        public double RunnerUpScore { get; set; }
        public bool IsDecoy { get; private set; }
        public bool IsContaminant { get; private set; }
        public DigestionParams DigestionParams { get; }
        public List<double> AllScores { get; internal set; }
        public Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> PeptidesToMatchingFragments { get; private set; }

        public IEnumerable<(int Notch, PeptideWithSetModifications Peptide)> BestMatchingPeptides
        {
            get
            {
                return _BestMatchingPeptides.OrderBy(p => p.Item2.FullSequence)
                    .ThenBy(p => p.Item2.Protein.Accession)
                    .ThenBy(p => p.Item2.OneBasedStartResidueInProtein);
            }
        }

        /// <summary>
        /// Used for Percolator output
        /// </summary>
        public double[] Features
        {
            get
            {
                return new[] { Math.Round(Score), Score - Math.Round(Score) };
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
            else if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
            {
                RunnerUpScore = newScore;
            }
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
            PsmTsvWriter.AddMatchedIonsData(s, psm == null ? null : psm.MatchedFragmentIons);
            PsmTsvWriter.AddMatchScoreData(s, psm);
            return s;
        }

        public void CalculateDeltaScore(double scoreCutoff)
        {
            DeltaScore = Score - Math.Max(RunnerUpScore, scoreCutoff);
        }

        public void SetFdrValues(double cumulativeTarget, double cumulativeDecoy, double qValue, double cumulativeTargetNotch, double cumulativeDecoyNotch, double qValueNotch, double maximumLikelihood, double eValue, double eScore, bool calculateEValue)
        {
            FdrInfo = new FdrInfo
            {
                CumulativeTarget = cumulativeTarget,
                CumulativeDecoy = cumulativeDecoy,
                QValue = qValue,
                CumulativeTargetNotch = cumulativeTargetNotch,
                CumulativeDecoyNotch = cumulativeDecoyNotch,
                QValueNotch = qValueNotch,
                MaximumLikelihood = maximumLikelihood,
                EScore = eScore,
                EValue = eValue,
                CalculateEValue = calculateEValue
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
                        _BestMatchingPeptides.RemoveAll(p => p.Pwsm.FullSequence == hit.Key && p.Pwsm.Protein.IsDecoy);
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

        /// <summary>
        /// This method is used by protein parsimony to remove PeptideWithSetModifications objects that have non-parsimonious protein associations
        /// </summary>
        public void TrimProteinMatches(HashSet<Protein> parsimoniousProteins)
        {
            if (IsDecoy)
            {
                if (_BestMatchingPeptides.Any(p => parsimoniousProteins.Contains(p.Pwsm.Protein) && p.Pwsm.Protein.IsDecoy))
                {
                    _BestMatchingPeptides.RemoveAll(p => !parsimoniousProteins.Contains(p.Item2.Protein));
                }
                // else do nothing
            }
            else
            {
                _BestMatchingPeptides.RemoveAll(p => !parsimoniousProteins.Contains(p.Item2.Protein));
            }

            ResolveAllAmbiguities();
        }

        /// <summary>
        /// This method is used by protein parsimony to add PeptideWithSetModifications objects for modification-agnostic parsimony
        /// </summary>
        public void AddProteinMatch((int, PeptideWithSetModifications) peptideWithNotch)
        {
            _BestMatchingPeptides.Add(peptideWithNotch);
            ResolveAllAmbiguities();
        }
    }
}