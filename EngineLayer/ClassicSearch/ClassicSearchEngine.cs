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
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly double[] MyScanPrecursorMasses;
        private readonly List<ProductType> ProductTypes;
        private readonly List<DissociationType> DissociationTypes;

        public ClassicSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<ModificationWithMass> variableModifications,
            List<ModificationWithMass> fixedModifications, List<Protein> proteinList, List<ProductType> productType, MassDiffAcceptor searchMode,
            CommonParameters commonParameters, List<string> nestedIds)
            : base(commonParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            Proteins = proteinList;
            SearchMode = searchMode;
            ProductTypes = productType;
            DissociationTypes = DetermineDissociationType(productType);
        }

        public static int x = 0;

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Getting ms2 scans...");

            double proteinsSearched = 0;
            int oldPercentProgress = 0;
            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(ProductTypes);

            // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
            var myLocks = new object[PeptideSpectralMatches.Length];
            for (int i = 0; i < myLocks.Length; i++)
            {
                myLocks[i] = new object();
            }

            Status("Performing classic search...");

            if (Proteins.Any())
            {
                Parallel.ForEach(Partitioner.Create(0, Proteins.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                    (partitionRange, loopState) =>
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
                        foreach (var peptide in Proteins[i].Digest(CommonParameters.DigestionParams, FixedModifications, VariableModifications))
                        {
                            var peptideTheorIons = peptide.GetTheoreticalFragments(ProductTypes);
                            var compactPeptide = peptide.CompactPeptide(terminusType);

                            foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(compactPeptide.MonoisotopicMassIncludingFixedMods, SearchMode))
                            {
                                var matchedIons = MatchFragmentIons(scan.TheScan.TheScan.MassSpectrum, peptideTheorIons, CommonParameters);

                                if (CommonParameters.AddCompIons)
                                {
                                    foreach (var dissociationType in DissociationTypes)
                                    {
                                        MzSpectrum complementarySpectrum = GenerateComplementarySpectrum(scan.TheScan.TheScan.MassSpectrum, scan.TheScan.PrecursorMass, dissociationType);
                                        matchedIons.AddRange(MatchFragmentIons(complementarySpectrum, peptideTheorIons, CommonParameters));
                                    }
                                }

                                double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons, 0);

                                bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;

                                // this is thread-safe because even if the score improves from another thread writing to this PSM,
                                // the lock combined with AddOrReplace method will ensure thread safety
                                if (meetsScoreCutoff || CommonParameters.CalculateEValue)
                                {
                                    // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
                                    lock (myLocks[scan.ScanIndex])
                                    {
                                        bool scoreImprovement = PeptideSpectralMatches[scan.ScanIndex] == null || (thisScore - PeptideSpectralMatches[scan.ScanIndex].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                                        if (scoreImprovement)
                                        {
                                            if (PeptideSpectralMatches[scan.ScanIndex] == null)
                                            {
                                                PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(compactPeptide, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, CommonParameters.DigestionParams);
                                            }
                                            else
                                            {
                                                PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(compactPeptide, thisScore, scan.Notch, CommonParameters.ReportAllAmbiguity);
                                            }

                                            //TODO: move this into the PeptideSpectralMatch constructor
                                            PeptideSpectralMatches[scan.ScanIndex].SetMatchedFragments(matchedIons);
                                        }

                                        if (CommonParameters.CalculateEValue)
                                        {
                                            PeptideSpectralMatches[scan.ScanIndex].AllScores.Add(thisScore);
                                        }
                                    }
                                }
                            }
                        }

                        // report search progress (proteins searched so far out of total proteins in database)
                        proteinsSearched++;
                        var percentProgress = (int)((proteinsSearched / Proteins.Count) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Performing classic search... ", NestedIds));
                        }
                    }
                });
            }

            // remove peptides below the score cutoff that were stored to calculate expectation values
            if (CommonParameters.CalculateEValue)
            {
                for (int i = 0; i < PeptideSpectralMatches.Length; i++)
                {
                    if (PeptideSpectralMatches[i] != null && PeptideSpectralMatches[i].Score < CommonParameters.ScoreCutoff)
                    {
                        PeptideSpectralMatches[i] = null;
                    }
                }
            }

            return new MetaMorpheusEngineResults(this);
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor searchMode)
        {
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass).ToList())
            {
                DoubleRange allowedInterval = allowedIntervalWithNotch.AllowedInterval;
                int scanIndex = GetFirstScanWithMassOverOrEqual(allowedInterval.Minimum);
                if (scanIndex < ArrayOfSortedMS2Scans.Length)
                {
                    var scanMass = MyScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedInterval.Maximum)
                    {
                        var theScan = ArrayOfSortedMS2Scans[scanIndex];
                        yield return new ScanWithIndexAndNotchInfo(theScan, allowedIntervalWithNotch.Notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == ArrayOfSortedMS2Scans.Length)
                        {
                            break;
                        }
                        scanMass = MyScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum)
        {
            int index = Array.BinarySearch(MyScanPrecursorMasses, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            // index of the first element that is larger than value
            return index;
        }
    }
}