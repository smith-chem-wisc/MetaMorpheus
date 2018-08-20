using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ClassicSearch
{
    public class ClassicSearchEngine : MetaMorpheusEngine
    {
        private readonly MassDiffAcceptor SearchMode;
        private readonly List<Protein> Proteins;
        private readonly List<ModificationWithMass> FixedModifications;
        private readonly List<ModificationWithMass> VariableModifications;
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches; // PeptideSpectralMatch[Scan]
        protected readonly PeptideSpectralMatch[][] DecoyPeptideSpectralMatches; // PeptideSpectralMatch[Database][Scan]
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly double[] MyScanPrecursorMasses;
        private readonly List<ProductType> ProductTypes;
        private readonly List<DissociationType> DissociationTypes;

        public ClassicSearchEngine(PeptideSpectralMatch[] globalPsms, PeptideSpectralMatch[][] decoyPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, List<ProductType> lp, MassDiffAcceptor searchMode, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            DecoyPeptideSpectralMatches = decoyPsms;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            Proteins = proteinList;
            SearchMode = searchMode;
            ProductTypes = lp;
            DissociationTypes = DetermineDissociationType(lp);
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Getting ms2 scans...");

            // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
            var scanLocks = new object[PeptideSpectralMatches.Length];
            for (int i = 0; i < scanLocks.Length; i++)
            {
                scanLocks[i] = new object();
            }

            string message = "Performing classic search...";
            Status(message);
            if (Proteins.Any())
            {
                Search(message, scanLocks);
            }

            if (commonParameters.DoGenerateShuffledDecoys)
            {
                message = "Performing shuffled decoy classic search...";
                Status(message);
                HashSet<string> forbiddenShuffledSequences = new HashSet<string>(PeptideSpectralMatches
                    .Where(x => x != null && x.Score > commonParameters.ScoreCutoff)
                    .Select(x => x.BaseSequence));

                for (int i = 0; i < commonParameters.NumDecoyDatabases; i++)
                {
                    message = "Performing shuffled decoy classic search " + (i + 1).ToString() + " of " + commonParameters.NumDecoyDatabases.ToString() + "...";
                    Status(message);
                    DecoyPeptideSpectralMatches[i] = new PeptideSpectralMatch[PeptideSpectralMatches.Length];
                    SearchShuffledDecoys(message, scanLocks, i, forbiddenShuffledSequences, new Random(i));
                }
            }

            // remove peptides below the score cutoff that were stored to calculate expectation values
            if (commonParameters.CalculateEValue)
            {
                for (int i = 0; i < PeptideSpectralMatches.Length; i++)
                {
                    if (PeptideSpectralMatches[i] != null && PeptideSpectralMatches[i].Score < commonParameters.ScoreCutoff)
                    {
                        PeptideSpectralMatches[i] = null;
                    }
                }
            }

            return new MetaMorpheusEngineResults(this);
        }

        private void Search(string message, object[] scanLocks)
        {
            double proteinsSearched = 0;
            int oldPercentProgress = 0;
            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(ProductTypes);

            Parallel.ForEach(Partitioner.Create(0, Proteins.Count), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (partitionRange, loopState) =>
            {
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                    foreach (var peptide in Proteins[i].Digest(commonParameters.DigestionParams, FixedModifications, VariableModifications))
                    {
                        var peptideTheorIons = peptide.GetTheoreticalFragments(ProductTypes);
                        var compactPeptide = peptide.CompactPeptide(terminusType);

                        foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(compactPeptide.MonoisotopicMassIncludingFixedMods, SearchMode))
                        {
                            Score(scan, peptideTheorIons, scanLocks, compactPeptide, peptide, false, 0);
                        }
                    }

                    // report search progress (proteins searched so far out of total proteins in database)
                    proteinsSearched++;
                    var percentProgress = (int)((proteinsSearched / Proteins.Count) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, message, nestedIds));
                    }
                }
            });
        }

        private void SearchShuffledDecoys(string message, object[] scanLocks, int decoyDatabase, HashSet<string> forbiddenShuffledSequences, Random random)
        {
            double psmsSearched = 0;
            int oldPercentProgress = 0;
            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(ProductTypes);

            Parallel.ForEach(Partitioner.Create(0, PeptideSpectralMatches.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (partitionRange, loopState) =>
            {
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    var psm = PeptideSpectralMatches[i];
                    if (psm == null) { continue; }
                    var decoyPeptide = ShufflePeptide(psm.Peptide, forbiddenShuffledSequences, random);
                    Score(GetScan(psm), psm.Peptide.GetTheoreticalFragments(ProductTypes), scanLocks, psm.Peptide.CompactPeptide(terminusType), psm.Peptide, true, decoyDatabase);

                    // report search progress (proteins searched so far out of total proteins in database)
                    psmsSearched++;
                    var percentProgress = (int)((psmsSearched / PeptideSpectralMatches.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, message, nestedIds));
                    }
                }
            });
        }

        private void Score(ScanWithIndexAndNotchInfo scan, List<TheoreticalFragmentIon> peptideTheorIons, object[] scanLocks,
            CompactPeptide compactPeptide, PeptideWithSetModifications peptide, bool doGenerateShuffledDecoys, int decoyDatabase)
        {
            var matchedIons = MatchFragmentIons(scan.TheScan.TheScan.MassSpectrum, peptideTheorIons, commonParameters);
            double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons, 0);
            bool meetsScoreCutoff = thisScore >= commonParameters.ScoreCutoff;

            // this is thread-safe because even if the score improves from another thread writing to this PSM,
            // the lock combined with AddOrReplace method will ensure thread safety
            if (meetsScoreCutoff || commonParameters.CalculateEValue)
            {
                // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
                lock (scanLocks[scan.ScanIndex])
                {
                    PeptideSpectralMatch[] psmArray = doGenerateShuffledDecoys ? DecoyPeptideSpectralMatches[decoyDatabase] : PeptideSpectralMatches;
                    PeptideSpectralMatch psm = psmArray[scan.ScanIndex];
                    bool scoreImprovement = psm == null || thisScore - psm.RunnerUpScore > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                    if (scoreImprovement)
                    {
                        if (psm == null)
                        {
                            psmArray[scan.ScanIndex] = new PeptideSpectralMatch(compactPeptide, peptide, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, matchedIons, commonParameters.DigestionParams);
                        }
                        else
                        {
                            psmArray[scan.ScanIndex].AddOrReplace(compactPeptide, peptide, thisScore, scan.Notch, matchedIons, commonParameters.ReportAllAmbiguity);
                        }
                    }

                    if (commonParameters.CalculateEValue)
                        psm.AllScores.Add(thisScore);
                }
            }
        }

        private ScanWithIndexAndNotchInfo GetScan(PeptideSpectralMatch psm)
        {
            return new ScanWithIndexAndNotchInfo(ArrayOfSortedMS2Scans[psm.ScanIndex], 0, psm.ScanIndex);
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor searchMode)
        {
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass).ToList())
            {
                DoubleRange allowedInterval = allowedIntervalWithNotch.AllowedInterval;
                int ScanIndex = GetFirstScanWithMassOverOrEqual(allowedInterval.Minimum);
                if (ScanIndex < ArrayOfSortedMS2Scans.Length)
                {
                    var scanMass = MyScanPrecursorMasses[ScanIndex];
                    while (scanMass <= allowedInterval.Maximum)
                    {
                        var TheScan = ArrayOfSortedMS2Scans[ScanIndex];
                        yield return new ScanWithIndexAndNotchInfo(TheScan, allowedIntervalWithNotch.Notch, ScanIndex);
                        ScanIndex++;
                        if (ScanIndex == ArrayOfSortedMS2Scans.Length)
                            break;
                        scanMass = MyScanPrecursorMasses[ScanIndex];
                    }
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum)
        {
            int index = Array.BinarySearch(MyScanPrecursorMasses, minimum);
            if (index < 0)
                index = ~index;

            // index of the first element that is larger than value
            return index;
        }

        /// <summary>
        /// Shuffles a peptide sequence, maintaining initiating methionine and cleavage sites
        /// Continues to shuffle until sequence is not in forbidden set (or until it's not a valid length and will be filtered later)
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="forbiddenSequences"></param>
        /// <param name="random"></param>
        /// <returns></returns>
        private PeptideWithSetModifications ShufflePeptide(PeptideWithSetModifications peptide, HashSet<string> forbiddenSequences, Random random)
        {
            while (true)
            {
                Dictionary<int, ModificationWithMass> decoyOneBasedModifications = new Dictionary<int, ModificationWithMass>();

                // Shuffle the one based indices and correct the start / end bases if important
                string seq = peptide.BaseSequence;
                bool keepFirstBase = peptide.OneBasedStartResidueInProtein == 1 && seq[0] == 'M'
                    || peptide.DigestionParams.Protease.SequencesInducingCleavage.Any(x => x.Item2 == TerminusType.C && seq.StartsWith(x.Item1));
                bool keepLastBase = peptide.OneBasedEndResidueInProtein == peptide.Protein.Length
                    && peptide.DigestionParams.Protease.SequencesInducingCleavage.Any(x => x.Item2 == TerminusType.N && seq.EndsWith(x.Item1));
                List<int> unshuffledOneBasedIndices = Enumerable.Range(peptide.OneBasedStartResidueInProtein, seq.Length).ToList();
                List<int> shuffledOneBasedProteinIndices = unshuffledOneBasedIndices.OrderBy(x => random.Next()).ToList();
                Dictionary<int, int> oneBasedIdxLookup = Enumerable.Range(0, seq.Length).ToDictionary(x => unshuffledOneBasedIndices[x], x => shuffledOneBasedProteinIndices[x]);
                if (keepFirstBase)
                {
                    shuffledOneBasedProteinIndices[shuffledOneBasedProteinIndices.IndexOf(peptide.OneBasedStartResidueInProtein)] = shuffledOneBasedProteinIndices[0];
                    shuffledOneBasedProteinIndices[0] = peptide.OneBasedStartResidueInProtein;
                }
                if (keepLastBase)
                {
                    shuffledOneBasedProteinIndices[shuffledOneBasedProteinIndices.IndexOf(peptide.OneBasedEndResidueInProtein)] = shuffledOneBasedProteinIndices[seq.Length - 1];
                    shuffledOneBasedProteinIndices[seq.Length - 1] = peptide.OneBasedEndResidueInProtein;
                }

                // Construct the decoy sequence and modification set
                // Append the decoy peptide if already in the database (and a valid length by digestion parameters)
                string decoyPeptideSeq = new string(Enumerable.Range(0, seq.Length)
                    .Select(x => seq[shuffledOneBasedProteinIndices[x] - peptide.OneBasedStartResidueInProtein]).ToArray());
                bool validTargetLength = decoyPeptideSeq.Length >= peptide.DigestionParams.MinPeptideLength && peptide.DigestionParams.MaxPeptideLength <= decoyPeptideSeq.Length;
                if (validTargetLength && !forbiddenSequences.Contains(decoyPeptideSeq.ToString()) || !validTargetLength)
                {
                    for (int j = 0; j < seq.Length; j++)
                    {
                        // if there is a modification at this position in the regular sequence, add it to the appropriate position in the shuffled one
                        if (peptide.AllModsOneIsNterminus.TryGetValue(shuffledOneBasedProteinIndices[j], out var mod))
                        {
                            decoyOneBasedModifications[j] = mod;
                        }
                    }
                    return new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, "DECOY_" + peptide.PeptideDescription, peptide.MissedCleavages, decoyOneBasedModifications, peptide.NumFixedMods);
                }
            }
        }
    }
}