using Chemistry;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using Omics;
using System;
using Omics.Digestion;
using EngineLayer.CrosslinkSearch;
using EngineLayer.SpectrumMatch;

namespace EngineLayer
{
    public abstract class SpectralMatch : IComparable<SpectralMatch>
    {
        public const double ToleranceForScoreDifferentiation = 1e-9;

        protected SpectralMatch(IBioPolymerWithSetMods peptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters, List<MatchedFragmentIon> matchedFragmentIons)
        {
            _BestMatchingBioPolymersWithSetMods = new List<SpectralMatchHypothesis>();
            ScanIndex = scanIndex;
            FullFilePath = scan.FullFilePath;
            ScanNumber = scan.OneBasedScanNumber;
            PrecursorScanNumber = scan.OneBasedPrecursorScanNumber;
            ScanRetentionTime = scan.RetentionTime;
            ScanExperimentalPeaks = scan.NumPeaks;
            PrecursorScanIntensity = scan.PrecursorIntensity;
            TotalIonCurrent = scan.TotalIonCurrent;
            ScanPrecursorCharge = scan.PrecursorCharge;
            ScanPrecursorMonoisotopicPeakMz = scan.PrecursorMonoisotopicPeakMz;
            ScanPrecursorMass = scan.PrecursorMass;
            PrecursorScanEnvelopePeakCount = scan.PrecursorEnvelopePeakCount;
            PrecursorFractionalIntensity = scan.PrecursorFractionalIntensity;
            DigestionParams = commonParameters.DigestionParams;
            NativeId = scan.NativeId;
            RunnerUpScore = commonParameters.ScoreCutoff;
            SpectralAngle = -1;
            IsobaricMassTagReporterIonIntensities = scan.IsobaricMassTagReporterIonIntensities;

            if (scan.TheScan is TimsDataScan timsScan)
            {
                ScanOneOverK0 = timsScan.OneOverK0;
            }
            else
            {
                ScanOneOverK0 = null; // this is only used for ion mobility data, so it can be null
            }

            AddOrReplace(peptide, score, notch, true, matchedFragmentIons);
        }

        public ChemicalFormula ModsChemicalFormula { get; private set; } // these fields will be null if they are ambiguous
        public string FullSequence { get; protected set; }
        public string EssentialSequence { get; protected set; }
        public int? Notch { get; private set; }
        public string BaseSequence { get; protected set; }
        public int? BioPolymerWithSetModsLength { get; private set; }
        public int? OneBasedStartResidue { get; private set; }
        public int? OneBasedEndResidue { get; private set; }
        public double? BioPolymerWithSetModsMonoisotopicMass { get; private set; }
        public int? ParentLength { get; private set; }
        public string Accession { get; private set; }
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
        public double PrecursorScanIntensity { get; }
        public int PrecursorScanEnvelopePeakCount { get; }
        public double PrecursorFractionalIntensity { get; }
        public double ScanPrecursorMass { get; }
        public double? ScanOneOverK0 { get; set; } // this is only used for ion mobility data, so it can be null
        public string FullFilePath { get; private set; }
        /// <summary>
        /// Refers to the index of the Ms2ScanWithSpecificMass in an array of Ms2ScansWithSpecificMass that is sorted by precursor mass
        /// </summary>
        public int ScanIndex { get; }
        public int NumDifferentMatchingPeptides { get { return _BestMatchingBioPolymersWithSetMods.Count; } }

        public FdrInfo FdrInfo
        {
            get => PsmFdrInfo;
            set => PsmFdrInfo = value;
        }
        public FdrInfo PsmFdrInfo { get; set; }
        public FdrInfo PeptideFdrInfo { get; set; }
        public FdrInfo GetFdrInfo(bool peptideLevel)
        {
            return peptideLevel ? PeptideFdrInfo : PsmFdrInfo;
        }

        public PsmData PsmData_forPEPandPercolator { get; set; }

        public double Score { get; private set; }
        public double SpectralAngle { get; set; }
        public string NativeId; // this is a property of the scan. used for mzID writing

        public double DeltaScore { get { return (Score - RunnerUpScore); } }

        public double RunnerUpScore { get; set; }
        public bool IsDecoy { get; private set; }
        public bool IsContaminant { get; private set; }

        //One-based positions in peptide that are covered by fragments on both sides of amino acids
        public List<int> FragmentCoveragePositionInPeptide { get; private set; }

        public List<double> PrecursorMassErrorDa
        {
            get
            {
                return this._BestMatchingBioPolymersWithSetMods.Select(p => Math.Round(this.ScanPrecursorMass - p.SpecificBioPolymer.MonoisotopicMass, 5))
                    .ToList();
            }
        }

        public List<double> PrecursorMassErrorPpm
        {
            get
            {
                return this._BestMatchingBioPolymersWithSetMods.Select(p => Math.Round((this.ScanPrecursorMass - p.SpecificBioPolymer.MonoisotopicMass) / p.SpecificBioPolymer.MonoisotopicMass * 1e6, 2)).ToList();
            }
        }

        #region GPTMD

        /// <summary>
        /// This property is only used by the GPTMD Task
        /// </summary>
        public MsDataScan Ms2Scan { get; private set; }

        /// <summary>
        /// This should only be called within the GPTMD Task
        /// </summary>
        /// <param name="scan"></param>
        public void SetMs2Scan(MsDataScan scan)
        {
            Ms2Scan = scan;
        }

        #endregion

        #region Search

        public IDigestionParams DigestionParams { get; }

        public static BioPolymerNotchFragmentIonComparer BioPolymerNotchFragmentIonComparer = new();

        protected List<SpectralMatchHypothesis> _BestMatchingBioPolymersWithSetMods;

        /// <summary>
        /// An array containing the intensities of the reporter ions for isobaric mass tags. 
        /// If multiplex quantification wasn't performed, this will be null
        /// </summary>
        public double[]? IsobaricMassTagReporterIonIntensities { get; private set; }

        public IEnumerable<SpectralMatchHypothesis> BestMatchingBioPolymersWithSetMods
        {
            get
            {
                // This property gets called frequently
                // It might be worth considering stashing the sorted list in a field instead of sorting every time

                // Order high (better matches) to low (worse matches)
                return _BestMatchingBioPolymersWithSetMods.OrderByDescending(p => p, BioPolymerNotchFragmentIonComparer);
            }
        }

        public void AddOrReplace(IBioPolymerWithSetMods pwsm, double newScore, int notch, bool reportAllAmbiguity, List<MatchedFragmentIon> matchedFragmentIons)
        {
            if (newScore - Score > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                _BestMatchingBioPolymersWithSetMods.Clear();
                _BestMatchingBioPolymersWithSetMods.Add(new(notch, pwsm, matchedFragmentIons, newScore));

                if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
                {
                    RunnerUpScore = Score;
                }
                Score = newScore;
                Notch = notch;
            }
            else if (newScore - Score > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                _BestMatchingBioPolymersWithSetMods.Add(new(notch, pwsm, matchedFragmentIons, newScore));

                if (Notch != notch)
                    Notch = null;
            }
            else if (newScore - RunnerUpScore > ToleranceForScoreDifferentiation)
            {
                RunnerUpScore = newScore;
            }
        }

        //PEP-Value analysis identifies ambiguous peptides with lower probability. These are removed from the bestmatchingpeptides dictionary, which lowers ambiguity.
        public void RemoveThisAmbiguousPeptide(SpectralMatchHypothesis tentativeSpectralMatch)
        {
            _BestMatchingBioPolymersWithSetMods.Remove(tentativeSpectralMatch);
            this.ResolveAllAmbiguities();
        }

        /// <summary>
        /// This method saves properties of this PSM for internal use. It is NOT used for any output.
        /// These resolved fields are (usually) null if there is more than one option.
        /// e.g., if this PSM can be explained by more than one base sequence, the BaseSequence property will be null
        /// </summary>
        public void ResolveAllAmbiguities()
        {
            // Order the BPWSM list for stability
            var bestMatches = BestMatchingBioPolymersWithSetMods.ToList();

            IsDecoy = bestMatches.Any(p => p.SpecificBioPolymer.Parent.IsDecoy);
            IsContaminant = bestMatches.Any(p => p.SpecificBioPolymer.Parent.IsContaminant);
            FullSequence = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.FullSequence)).ResolvedValue;
            BaseSequence = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.BaseSequence)).ResolvedValue;
            BioPolymerWithSetModsLength = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.Length)).ResolvedValue;
            OneBasedStartResidue = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.OneBasedStartResidue)).ResolvedValue;
            OneBasedEndResidue = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.OneBasedEndResidue)).ResolvedValue;
            ParentLength = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.Parent.Length)).ResolvedValue;
            BioPolymerWithSetModsMonoisotopicMass = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.MonoisotopicMass)).ResolvedValue;
            Accession = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.Parent.Accession)).ResolvedValue;
            Organism = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.Parent.Organism)).ResolvedValue;
            ModsIdentified = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.AllModsOneIsNterminus)).ResolvedValue;
            ModsChemicalFormula = PsmTsvWriter.Resolve(bestMatches.Select(b => b.SpecificBioPolymer.AllModsOneIsNterminus.Select(c => (c.Value)))).ResolvedValue;
            Notch = PsmTsvWriter.Resolve(bestMatches.Select(b => b.Notch)).ResolvedValue;
            if (Notch == null) // notch ambiguous resolve q-value notch fields if possible. 
            {
                if (bestMatches.First().QValueNotch != null) // Psm Level FDR has run
                {
                    var qNotch = PsmTsvWriter.Resolve(bestMatches.Select(p => p.QValueNotch!.Value)).ResolvedValue;
                    if (qNotch != null) // Q-value notches are identical even if the notches are not
                    {
                        PsmFdrInfo.QValueNotch = qNotch.Value;
                        PsmFdrInfo.CumulativeTargetNotch = PsmTsvWriter.Resolve(bestMatches.Select(p => p.CumulativeTargetNotch!.Value)).ResolvedValue ?? 0;
                        PsmFdrInfo.CumulativeDecoyNotch = PsmTsvWriter.Resolve(bestMatches.Select(p => p.CumulativeDecoyNotch!.Value)).ResolvedValue ?? 0;
                    }
                }

                if (bestMatches.First().PeptideQValueNotch != null) // Peptide Level FDR has run
                {
                    var pepQNotch = PsmTsvWriter.Resolve(bestMatches.Select(p => p.PeptideQValueNotch!.Value)).ResolvedValue;
                    if (pepQNotch != null) // Q-value notches are identical even if the notches are not
                    {
                        PeptideFdrInfo.QValueNotch = pepQNotch.Value;
                        PeptideFdrInfo.CumulativeTargetNotch = PsmTsvWriter.Resolve(bestMatches.Select(p => p.PeptideCumulativeTargetNotch!.Value)).ResolvedValue ?? 0;
                        PeptideFdrInfo.CumulativeDecoyNotch = PsmTsvWriter.Resolve(bestMatches.Select(p => p.PeptideCumulativeDecoyNotch!.Value)).ResolvedValue ?? 0;
                    }
                }
            }

            //if the PSM matches a target and a decoy and they are the SAME SEQUENCE, remove the decoy
            if (IsDecoy)
            {
                bool removedPeptides = false;
                var hits = _BestMatchingBioPolymersWithSetMods.GroupBy(p => p.FullSequence);

                foreach (var hit in hits)
                {
                    if (hit.Any(p => p.IsDecoy) && hit.Any(p => !p.IsDecoy))
                    {
                        // at least one peptide with this sequence is a target and at least one is a decoy
                        // remove the decoys with this sequence
                        _BestMatchingBioPolymersWithSetMods.RemoveAll(p => p.FullSequence == hit.Key && p.IsDecoy);
                        removedPeptides = true;
                    }
                }

                if (removedPeptides)
                {
                    ResolveAllAmbiguities();
                }
            }

            // Technically, different peptide options for this PSM can have different matched ions
            // However, writing out all the matched ions for all the peptide options would break excel
            // Instead, we set MatchedFragmentIons as the ions matched to the best peptide option
            if (this is CrosslinkSpectralMatch) // CrosslinkSpectralMatch has its own way of handling this, however, this method of retrieving the "First" item in a dictionary is problematic and should be revisted at some point
                MatchedFragmentIons = bestMatches.First().MatchedIons;
            else 
                MatchedFragmentIons = bestMatches.First().MatchedIons;
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

        #endregion

        #region IO

        public static string GetTabSeparatedHeader(bool includeOneOverK0Column = false)
        {
            return string.Join("\t", DataDictionary(null, null, includeOneOverK0Column: includeOneOverK0Column).Keys);
        }

        public override string ToString()
        {
            return ToString(new Dictionary<string, int>());
        }

        public string ToString(IReadOnlyDictionary<string, int> ModstoWritePruned, bool writePeptideLevelFdr = false, bool includeOneOverK0Column = false)
        {
            return string.Join("\t", DataDictionary(this, ModstoWritePruned, writePeptideLevelFdr, includeOneOverK0Column).Values.Select(v => v?.Trim() ?? string.Empty));
        }

        public static Dictionary<string, string> DataDictionary(SpectralMatch psm, IReadOnlyDictionary<string, int> ModsToWritePruned, bool writePeptideLevelFdr = false, bool includeOneOverK0Column = false)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddBasicMatchData(s, psm, includeOneOverK0Column);
            PsmTsvWriter.AddPeptideSequenceData(s, psm, ModsToWritePruned);
            PsmTsvWriter.AddMatchedIonsData(s, psm?.MatchedFragmentIons);
            PsmTsvWriter.AddMatchScoreData(s, psm, writePeptideLevelFdr);
            return s;
        }

        #endregion

        #region Parsimony

        /// <summary>
        /// This method is used by protein parsimony to remove PeptideWithSetModifications objects that have non-parsimonious protein associations
        /// </summary>
        public void TrimProteinMatches(HashSet<IBioPolymer> parsimoniousProteins)
        {
            if (IsDecoy)
            {
                if (_BestMatchingBioPolymersWithSetMods.Any(p => parsimoniousProteins.Contains(p.SpecificBioPolymer.Parent) && p.SpecificBioPolymer.Parent.IsDecoy))
                {
                    _BestMatchingBioPolymersWithSetMods.RemoveAll(p => !parsimoniousProteins.Contains(p.SpecificBioPolymer.Parent));
                }
                // else do nothing
            }
            else
            {
                _BestMatchingBioPolymersWithSetMods.RemoveAll(p => !parsimoniousProteins.Contains(p.SpecificBioPolymer.Parent));
            }

            ResolveAllAmbiguities();
        }

        /// <summary>
        /// This method is used by protein parsimony to add PeptideWithSetModifications objects for modification-agnostic parsimony
        /// </summary>
        public void AddProteinMatch(SpectralMatchHypothesis tentativeSpectralMatch)
        {
            if (!_BestMatchingBioPolymersWithSetMods.Contains(tentativeSpectralMatch))
            {
                _BestMatchingBioPolymersWithSetMods.Add(tentativeSpectralMatch); 
                ResolveAllAmbiguities();
            }
        }

        #endregion

        #region Silac

        protected SpectralMatch(SpectralMatch psm, List<SpectralMatchHypothesis> bestMatchingPeptides)
        {
            _BestMatchingBioPolymersWithSetMods = bestMatchingPeptides;
            BaseSequence = PsmTsvWriter.Resolve(bestMatchingPeptides.Select(b => b.SpecificBioPolymer.BaseSequence)).ResolvedValue;
            FullSequence = PsmTsvWriter.Resolve(bestMatchingPeptides.Select(b => b.SpecificBioPolymer.FullSequence)).ResolvedValue;

            ModsChemicalFormula = psm.ModsChemicalFormula;
            Notch = psm.Notch;
            BioPolymerWithSetModsLength = psm.BioPolymerWithSetModsLength;
            OneBasedStartResidue = psm.OneBasedStartResidue;
            OneBasedEndResidue = psm.OneBasedEndResidue;
            BioPolymerWithSetModsMonoisotopicMass = psm.BioPolymerWithSetModsMonoisotopicMass;
            ParentLength = psm.ParentLength;
            Accession = psm.Accession;
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
            RunnerUpScore = psm.RunnerUpScore;
            IsDecoy = psm.IsDecoy;
            IsContaminant = psm.IsContaminant;
            DigestionParams = psm.DigestionParams;
            SpectralAngle = psm.SpectralAngle;
        }

        #endregion

        #region FDR

        private string _chimeraIdString;
        public string ChimeraIdString => _chimeraIdString ??= $"{ScanNumber}{FullFilePath}{PrecursorScanNumber}";

        /// <summary>
        /// Returns an integer representing the longest continuous number of residues in the match covered on both sides by fragment ions
        /// </summary>
        public static int GetLongestIonSeriesBidirectional(List<MatchedFragmentIon> matchedFragments, IBioPolymerWithSetMods peptide)
        {
            List<int> maxDiffs = new List<int> { 1 };
            if (matchedFragments != null && matchedFragments.Count != 0)
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

        /// <summary>
        /// Returns an integer representing the longest continuous number of residues in the match covered on both sides by fragment ions
        /// </summary>
        public static int GetLongestIonSeriesBidirectional(SpectralMatchHypothesis hypothesis) => GetLongestIonSeriesBidirectional(hypothesis.MatchedIons, hypothesis.SpecificBioPolymer);

        /// <summary>
        /// Determine the Fragment Coverage the PSM
        /// Assigns fragment coverage indices for the PSM and the protein based on Amino Acid Position in Matched Ion Fragments
        /// </summary>
        public void GetAminoAcidCoverage()
        {
            if (string.IsNullOrEmpty(this.BaseSequence) ||
                !this.MatchedFragmentIons.Any()) return;
            //Pull C terminal and N terminal Fragments and amino acid numbers
            var nTermFragmentAAPositions = this.MatchedFragmentIons.Where(p =>
                    p.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime)
                .Select(j => j.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

            var cTermFragmentAAPositions = this.MatchedFragmentIons.Where(p =>
                    p.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime)
                .Select(j => j.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

            //Create a hashset to store the covered amino acid positions
            HashSet<int> fragmentCoveredAminoAcids = new();

            //Check N term frags first
            if (nTermFragmentAAPositions.Any())
            {
                nTermFragmentAAPositions.Sort();

                //if the final NFragment is present, last AA is covered
                if (nTermFragmentAAPositions.Contains(this.BaseSequence.Length - 1))
                {
                    fragmentCoveredAminoAcids.Add(this.BaseSequence.Length);
                }

                // if the first NFragment is present, first AA is covered
                if (nTermFragmentAAPositions.Contains(1))
                {
                    fragmentCoveredAminoAcids.Add(1);
                }

                //Check all amino acids except for the last one in the list
                for (int i = 0; i < nTermFragmentAAPositions.Count - 1; i++)
                {
                    //sequential AA, second one is covered
                    if (nTermFragmentAAPositions[i + 1] - nTermFragmentAAPositions[i] == 1)
                    {
                        fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1]);
                    }

                    //check to see if the position is covered from both directions, inclusive
                    if (cTermFragmentAAPositions.Contains(nTermFragmentAAPositions[i + 1]))
                    {
                        fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1]);
                    }

                    //check to see if the position is covered from both directions, exclusive
                    if (cTermFragmentAAPositions.Contains(nTermFragmentAAPositions[i + 1] + 2))
                    {
                        fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1] + 1);
                    }
                }

            }

            //Check C term frags
            if (cTermFragmentAAPositions.Any())
            {
                cTermFragmentAAPositions.Sort();

                //if the second AA is present, the first AA is covered
                if (cTermFragmentAAPositions.Contains(2))
                {
                    fragmentCoveredAminoAcids.Add(1);
                }

                //if the last AA is present, the final AA is covered
                if (cTermFragmentAAPositions.Contains(this.BaseSequence.Length))
                {
                    fragmentCoveredAminoAcids.Add(this.BaseSequence.Length);
                }

                //check all amino acids except for the last one in the list
                for (int i = 0; i < cTermFragmentAAPositions.Count - 1; i++)
                {
                    //sequential AA, the first one is covered
                    if (cTermFragmentAAPositions[i + 1] - cTermFragmentAAPositions[i] == 1)
                    {
                        fragmentCoveredAminoAcids.Add(cTermFragmentAAPositions[i]);
                    }
                }
            }

            //store in PSM
            var fragmentCoveredAminoAcidsList = fragmentCoveredAminoAcids.ToList();
            fragmentCoveredAminoAcidsList.Sort();
            this.FragmentCoveragePositionInPeptide = fragmentCoveredAminoAcidsList;
        }

        /// <summary>
        /// Gets the count of complementary ions in the PSM, or those that have been observed on both the N and C termini for the same residue
        /// </summary>
        public static int GetCountComplementaryIons(List<MatchedFragmentIon> matchedFragments, IBioPolymerWithSetMods peptide)
        {
            if (matchedFragments != null && matchedFragments.Count != 0)
            {
                var nIons = matchedFragments.Where(f => f.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime).Select(f => f.NeutralTheoreticalProduct.FragmentNumber);
                var cIons = matchedFragments.Where(f => f.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime).Select(f => peptide.BaseSequence.Length - f.NeutralTheoreticalProduct.FragmentNumber);

                return nIons.Intersect(cIons).Count();
            }

            return 0;
        }

        /// <summary>
        /// Gets the count of complementary ions in the PSM, or those that have been observed on both the N and C termini for the same residue
        /// </summary>
        public static int GetCountComplementaryIons(SpectralMatchHypothesis hypothesis) => GetCountComplementaryIons(hypothesis.MatchedIons, hypothesis.SpecificBioPolymer);

        #endregion

        /// <summary>
        /// There are a few key locations in MetaMorpheus where we want to have psms sorted in a consistent manner.
        /// These are for q-value determination and for when we write the psms to psmtsv. 
        /// </summary>
        /// <param name="otherPsm"></param>
        /// <returns></returns>
        public int CompareTo(SpectralMatch otherPsm)
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