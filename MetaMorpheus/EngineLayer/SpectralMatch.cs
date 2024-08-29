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

namespace EngineLayer
{
    public abstract class SpectralMatch : IComparable<SpectralMatch>
    {
        public const double ToleranceForScoreDifferentiation = 1e-9;

        protected SpectralMatch(IBioPolymerWithSetMods peptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters, List<MatchedFragmentIon> matchedFragmentIons, double xcorr = 0)
        {
            _BestMatchingBioPolymersWithSetMods = new List<(int, IBioPolymerWithSetMods)>();
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
            BioPolymersWithSetModsToMatchingFragments = new Dictionary<IBioPolymerWithSetMods, List<MatchedFragmentIon>>();
            Xcorr = xcorr;
            NativeId = scan.NativeId;
            RunnerUpScore = commonParameters.ScoreCutoff;
            MsDataScan = scan.TheScan;
            SpectralAngle = -1;

            AddOrReplace(peptide, score, notch, true, matchedFragmentIons, xcorr);
        }

        public MsDataScan MsDataScan { get; set; }
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
        public string FullFilePath { get; private set; }
        public int ScanIndex { get; }
        public int NumDifferentMatchingPeptides { get { return _BestMatchingBioPolymersWithSetMods.Count; } }

        public FdrInfo FdrInfo
        {
            get => PsmFdrInfo;
            set => PsmFdrInfo = value;

        }
        public FdrInfo PsmFdrInfo { get;  set; }
        public FdrInfo PeptideFdrInfo { get;  set; }
        public FdrInfo GetFdrInfo(bool peptideLevel)
        {
            return peptideLevel ? PeptideFdrInfo : PsmFdrInfo;
        }

        public PsmData PsmData_forPEPandPercolator { get; set; }

        public double Score { get; private set; }
        public double Xcorr;
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
                return this._BestMatchingBioPolymersWithSetMods.Select(p => Math.Round(this.ScanPrecursorMass - p.Pwsm.MonoisotopicMass, 5))
                    .ToList();
            }
        }

        public List<double> PrecursorMassErrorPpm
        {
            get
            {
                return this._BestMatchingBioPolymersWithSetMods.Select(p => Math.Round((this.ScanPrecursorMass - p.Pwsm.MonoisotopicMass) / p.Pwsm.MonoisotopicMass * 1e6, 2)).ToList();
            }
        }

        #region Search
        public DigestionParams DigestionParams { get; }
        public Dictionary<IBioPolymerWithSetMods, List<MatchedFragmentIon>> BioPolymersWithSetModsToMatchingFragments { get; private set; }

        protected List<(int Notch, IBioPolymerWithSetMods Pwsm)> _BestMatchingBioPolymersWithSetMods;
        public IEnumerable<(int Notch, IBioPolymerWithSetMods Peptide)> BestMatchingBioPolymersWithSetMods
        {
            get
            {
                return _BestMatchingBioPolymersWithSetMods.OrderBy(p => p.Pwsm.FullSequence)
                    .ThenBy(p => p.Pwsm.Parent.Accession)
                    .ThenBy(p => p.Pwsm.OneBasedStartResidue);
            }
        }

        public void AddOrReplace(IBioPolymerWithSetMods pwsm, double newScore, int notch, bool reportAllAmbiguity, List<MatchedFragmentIon> matchedFragmentIons, double newXcorr)
        {
            if (newScore - Score > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                _BestMatchingBioPolymersWithSetMods.Clear();
                _BestMatchingBioPolymersWithSetMods.Add((notch, pwsm));

                if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
                {
                    RunnerUpScore = Score;
                }

                Score = newScore;
                Xcorr = newXcorr;

                BioPolymersWithSetModsToMatchingFragments.Clear();
                BioPolymersWithSetModsToMatchingFragments.Add(pwsm, matchedFragmentIons);
            }
            else if (newScore - Score > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                _BestMatchingBioPolymersWithSetMods.Add((notch, pwsm));
                BioPolymersWithSetModsToMatchingFragments.TryAdd(pwsm, matchedFragmentIons);
            }
            else if (newScore - RunnerUpScore > ToleranceForScoreDifferentiation)
            {
                RunnerUpScore = newScore;
            }
        }

        //PEP-Value analysis identifies ambiguous peptides with lower probability. These are removed from the bestmatchingpeptides dictionary, which lowers ambiguity.
        public void RemoveThisAmbiguousPeptide(int notch, IBioPolymerWithSetMods pwsm)
        {
            _BestMatchingBioPolymersWithSetMods.Remove((notch, pwsm));
            if (!_BestMatchingBioPolymersWithSetMods.Any(x => x.Pwsm.Equals(pwsm)))
            {
                BioPolymersWithSetModsToMatchingFragments.Remove(pwsm);
            }
            this.ResolveAllAmbiguities();
        }

        /// <summary>
        /// This method saves properties of this PSM for internal use. It is NOT used for any output.
        /// These resolved fields are (usually) null if there is more than one option.
        /// e.g., if this PSM can be explained by more than one base sequence, the BaseSequence property will be null
        /// </summary>
        public void ResolveAllAmbiguities()
        {
            IsDecoy = _BestMatchingBioPolymersWithSetMods.Any(p => p.Pwsm.Parent.IsDecoy);
            IsContaminant = _BestMatchingBioPolymersWithSetMods.Any(p => p.Pwsm.Parent.IsContaminant);
            FullSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.FullSequence)).ResolvedValue;
            BaseSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.BaseSequence)).ResolvedValue;
            BioPolymerWithSetModsLength = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.Length)).ResolvedValue;
            OneBasedStartResidue = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.OneBasedStartResidue)).ResolvedValue;
            OneBasedEndResidue = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.OneBasedEndResidue)).ResolvedValue;
            ParentLength = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.Parent.Length)).ResolvedValue;
            BioPolymerWithSetModsMonoisotopicMass = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.MonoisotopicMass)).ResolvedValue;
            Accession = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.Parent.Accession)).ResolvedValue;
            Organism = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.Parent.Organism)).ResolvedValue;
            ModsIdentified = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.AllModsOneIsNterminus)).ResolvedValue;
            ModsChemicalFormula = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.AllModsOneIsNterminus.Select(c => (c.Value)))).ResolvedValue;
            Notch = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Notch)).ResolvedValue;

            //if the PSM matches a target and a decoy and they are the SAME SEQUENCE, remove the decoy
            if (IsDecoy)
            {
                bool removedPeptides = false;
                var hits = _BestMatchingBioPolymersWithSetMods.GroupBy(p => p.Pwsm.FullSequence);

                foreach (var hit in hits)
                {
                    if (hit.Any(p => p.Pwsm.Parent.IsDecoy) && hit.Any(p => !p.Pwsm.Parent.IsDecoy))
                    {
                        // at least one peptide with this sequence is a target and at least one is a decoy
                        // remove the decoys with this sequence
                        var pwsmToRemove = _BestMatchingBioPolymersWithSetMods.Where(p => p.Pwsm.FullSequence == hit.Key && p.Pwsm.Parent.IsDecoy).ToList();
                        _BestMatchingBioPolymersWithSetMods.RemoveAll(p => p.Pwsm.FullSequence == hit.Key && p.Pwsm.Parent.IsDecoy);
                        foreach ((int, IBioPolymerWithSetMods) pwsm in pwsmToRemove)
                        {
                            BioPolymersWithSetModsToMatchingFragments.Remove(pwsm.Item2);
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
            MatchedFragmentIons = BioPolymersWithSetModsToMatchingFragments.First().Value;
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
        /// This method is used to compute qValue etc for the inverted set of psms
        /// We neither compute nor calculated cumulativeTarget, cumulativeDecoy, etc for the inverted set.
        /// CumulativeTarget, CumulativeDecoy, CumulativeTargetNotch, CumulativeDecoyNotch, were computed
        /// for the non-inverted set.   We don't want to use them for the inverted set.
        /// </summary>
        /// <param name="qValue"></param>
        /// <param name="qValueNotch"></param>
        /// <param name="pep"></param>
        /// <param name="pepQValue"></param>
        public void SetQandPEPvalues(double qValue, double qValueNotch, double pep, double pepQValue)
        {
            FdrInfo.QValue = qValue;
            FdrInfo.QValueNotch = qValueNotch;
            FdrInfo.PEP = pep;
            FdrInfo.PEP_QValue = pepQValue;
        }



        #endregion

        #region IO

        public static string GetTabSeparatedHeader()
        {
            return string.Join("\t", DataDictionary(null, null).Keys);
        }

        public override string ToString()
        {
            return ToString(new Dictionary<string, int>());
        }

        public string ToString(IReadOnlyDictionary<string, int> ModstoWritePruned, bool writePeptideLevelFdr = false)
        {
            return string.Join("\t", DataDictionary(this, ModstoWritePruned, writePeptideLevelFdr).Values);
        }

        public static Dictionary<string, string> DataDictionary(SpectralMatch psm, IReadOnlyDictionary<string, int> ModsToWritePruned, bool writePeptideLevelFdr = false)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddBasicMatchData(s, psm);
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
        public void TrimProteinMatches(HashSet<Protein> parsimoniousProteins)
        {
            if (IsDecoy)
            {
                if (_BestMatchingBioPolymersWithSetMods.Any(p => parsimoniousProteins.Contains(p.Pwsm.Parent) && p.Pwsm.Parent.IsDecoy))
                {
                    _BestMatchingBioPolymersWithSetMods.RemoveAll(p => !parsimoniousProteins.Contains(p.Pwsm.Parent));
                }
                // else do nothing
            }
            else
            {
                _BestMatchingBioPolymersWithSetMods.RemoveAll(p => !parsimoniousProteins.Contains(p.Pwsm.Parent));
            }

            ResolveAllAmbiguities();
        }

        /// <summary>
        /// This method is used by protein parsimony to add PeptideWithSetModifications objects for modification-agnostic parsimony
        /// </summary>
        public void AddProteinMatch((int, IBioPolymerWithSetMods) peptideWithNotch, List<MatchedFragmentIon> mfi)
        {
            if (!_BestMatchingBioPolymersWithSetMods.Select(p => p.Pwsm).Contains(peptideWithNotch.Item2))
            {
                _BestMatchingBioPolymersWithSetMods.Add(peptideWithNotch);
                if (!BioPolymersWithSetModsToMatchingFragments.ContainsKey(peptideWithNotch.Item2))
                {
                    BioPolymersWithSetModsToMatchingFragments.Add(peptideWithNotch.Item2, mfi);
                }
                ResolveAllAmbiguities();
            }
        }

        #endregion

        #region Silac

        protected SpectralMatch(SpectralMatch psm, List<(int Notch, IBioPolymerWithSetMods Peptide)> bestMatchingPeptides)
        {
            _BestMatchingBioPolymersWithSetMods = bestMatchingPeptides;
            BaseSequence = PsmTsvWriter.Resolve(bestMatchingPeptides.Select(b => b.Peptide.BaseSequence)).ResolvedValue;
            FullSequence = PsmTsvWriter.Resolve(bestMatchingPeptides.Select(b => b.Peptide.FullSequence)).ResolvedValue;

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
            Xcorr = psm.Xcorr;
            RunnerUpScore = psm.RunnerUpScore;
            IsDecoy = psm.IsDecoy;
            IsContaminant = psm.IsContaminant;
            DigestionParams = psm.DigestionParams;
            BioPolymersWithSetModsToMatchingFragments = psm.BioPolymersWithSetModsToMatchingFragments;
            SpectralAngle = psm.SpectralAngle;
        }

        #endregion

        #region FDR

        private string _chimeraIdString;
        public string ChimeraIdString => _chimeraIdString ??= $"{ScanNumber}{FullFilePath}{PrecursorScanNumber}";

        /// <summary>
        /// Returns an integer representing the longest continuous number of residues in the match covered on both sides by fragment ions
        /// </summary>
        /// <param name="PeptidesToMatchingFragments"></param>
        /// <param name="peptide"></param>
        /// <returns></returns>
        public static int GetLongestIonSeriesBidirectional(Dictionary<IBioPolymerWithSetMods, List<MatchedFragmentIon>> PeptidesToMatchingFragments, IBioPolymerWithSetMods peptide)
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
                    p.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N)
                .Select(j => j.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

            var cTermFragmentAAPositions = this.MatchedFragmentIons.Where(p =>
                    p.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C)
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

        public static int GetCountComplementaryIons(Dictionary<IBioPolymerWithSetMods, List<MatchedFragmentIon>> PeptidesToMatchingFragments, IBioPolymerWithSetMods peptide)
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