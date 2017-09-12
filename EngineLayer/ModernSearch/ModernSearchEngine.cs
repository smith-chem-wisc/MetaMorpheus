using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ModernSearch
{
    public class ModernSearchEngine : MetaMorpheusEngine
    {
        #region Protected Fields

        protected readonly List<int>[] fragmentIndex;
        protected readonly float[] keys;
        protected readonly Psm[] globalPsms;
        protected readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;
        protected readonly List<CompactPeptide> peptideIndex;
        protected readonly List<ProductType> lp;
        protected readonly int currentPartition;
        protected readonly CommonParameters CommonParameters;
        protected readonly bool addCompIons;
        protected readonly MassDiffAcceptor massDiffAcceptors;
        protected readonly List<DissociationType> dissociationTypes;

        #endregion Protected Fields

        #region Private Fields

        private const double tolForScoreImprovement = 1e-9;

        #endregion Private Fields

        #region Public Constructors

        public ModernSearchEngine(Psm[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters CommonParameters, bool addCompIons, MassDiffAcceptor massDiffAcceptors, List<string> nestedIds) : base(nestedIds)
        {
            this.globalPsms = globalPsms;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;
            this.lp = lp;
            this.currentPartition = currentPartition + 1;
            this.CommonParameters = CommonParameters;
            this.addCompIons = addCompIons;
            this.massDiffAcceptors = massDiffAcceptors;
            this.dissociationTypes = DetermineDissociationType(lp);
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In modern search engine..." + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds);

            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;

            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                List<CompactPeptide> bestPeptides;
                double bestScores;
                List<int> bestNotches;
                double[] fullPeptideScores = new double[peptideIndexCount];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = listOfSortedms2Scans[i];
                    var thisScanprecursorMass = thisScan.PrecursorMass;
                    Array.Clear(fullPeptideScores, 0, peptideIndexCount);
                    double thePrecursorMass = thisScan.PrecursorMass;
                    CalculatePeptideScores(thisScan.TheScan, fullPeptideScores, thePrecursorMass);

                    bestPeptides = null;
                    bestScores = 0;
                    bestNotches = null;

                    for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                    {
                        var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                        if (consideredScore > CommonParameters.ScoreCutoff) //intentionally high. 99.9% of 4-mers are present in a given UniProt database. This saves considerable time
                        {
                            CompactPeptide candidatePeptide = peptideIndex[possibleWinningPeptideIndex];

                            // Check if makes sense to add due to peptidescore!
                            var searchMode = massDiffAcceptors;
                            double currentBestScore = bestScores;
                            if (currentBestScore > 1)
                            {
                                // Existed! Need to compare with old match
                                if ((Math.Abs(currentBestScore - consideredScore) < tolForScoreImprovement) && (CommonParameters.ReportAllAmbiguity || bestPeptides.Count == 0))
                                {
                                    // Score is same, need to see if accepts and if prefer the new one
                                    int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                    if (notch >= 0)
                                    {
                                        bestPeptides.Add(candidatePeptide);
                                        bestNotches.Add(notch);
                                    }
                                }
                                else if (currentBestScore < consideredScore)
                                {
                                    // Score is better, only make sure it is acceptable
                                    int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                    if (notch >= 0)
                                    {
                                        bestPeptides = new List<CompactPeptide> { candidatePeptide };
                                        bestScores = consideredScore;
                                        bestNotches = new List<int> { notch };
                                    }
                                }
                            }
                            // Did not exist! Only make sure that it is acceptable
                            else
                            {
                                int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                if (notch >= 0)
                                {
                                    bestPeptides = new List<CompactPeptide> { candidatePeptide };
                                    bestScores = consideredScore;
                                    bestNotches = new List<int> { notch };
                                }
                            }
                        }
                    }
                    if (bestPeptides != null)
                    {
                        int startIndex = 0;

                        if (globalPsms[i] == null)
                        {
                            globalPsms[i] = new Psm(bestPeptides[0], bestNotches[0], bestScores, i, thisScan, CommonParameters.ExcelCompatible);
                            startIndex = 1;
                        }

                        for (int k = startIndex; k < bestPeptides.Count; k++)
                        {
                            globalPsms[i].AddOrReplace(bestPeptides[k], bestScores, bestNotches[k], CommonParameters.ReportAllAmbiguity);
                        }
                    }
                }
                lock (outputObject)
                {
                    scansSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)scansSeen / (listOfSortedms2ScansLength) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In modern search loop" + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
                        old_progress = new_progress;
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        protected void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores, double thePrecursorMass)
        {
            //create previous variables to determine if peaks can be sequestered
            double previousTheAdd = 1 + spectrum.MassSpectrum.YArray[0] / spectrum.TotalIonCurrent;
            double previousExperimentalPeakInDaltons = spectrum.MassSpectrum.XArray[0] - Constants.protonMass;
            double previousMinRange = CommonParameters.ProductMassTolerance.GetMinimumValue(previousExperimentalPeakInDaltons);
            double previousMaxRange;
            //search observed peaks
            for (int i = 1; i < spectrum.MassSpectrum.Size; i++)
            {
                double experimentalPeakInDaltons = spectrum.MassSpectrum.XArray[i] - Constants.protonMass;
                if (CommonParameters.ProductMassTolerance.Within(previousExperimentalPeakInDaltons, experimentalPeakInDaltons))
                {
                    previousTheAdd += spectrum.MassSpectrum.YArray[i] / spectrum.TotalIonCurrent; //open to debate, currently sum intensities of all peaks within tolerance like it was low res
                }
                else
                {
                    previousMaxRange = CommonParameters.ProductMassTolerance.GetMaximumValue(previousExperimentalPeakInDaltons);
                    FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);
                    previousTheAdd = 1 + spectrum.MassSpectrum.YArray[i] / spectrum.TotalIonCurrent;
                    previousMinRange = CommonParameters.ProductMassTolerance.GetMinimumValue(experimentalPeakInDaltons);
                }
                previousExperimentalPeakInDaltons = experimentalPeakInDaltons;
            }
            previousMaxRange = CommonParameters.ProductMassTolerance.GetMaximumValue(previousExperimentalPeakInDaltons);
            FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);

            //generate experimental complementary ions if specified
            if (addCompIons)
            {
                //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the mass assuming that they're z=1
                int numCompIons = spectrum.MassSpectrum.Size;
                (double mass, double intensity)[] complementaryIons = new(double mass, double intensity)[numCompIons];

                foreach (DissociationType dissociationType in dissociationTypes)
                {
                    if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
                    {
                        double massShiftForComplementaryConversion = thePrecursorMass + protonMassShift; //mass shift needed to reobtain the original product ion for calculating tolerance
                        for (int i = numCompIons - 1; i >= 0; i--)
                            complementaryIons[numCompIons-i-1] = (massShiftForComplementaryConversion - spectrum.MassSpectrum.XArray[i], spectrum.MassSpectrum.YArray[i]);

                        //propogation of error from precursor mass and complementary product mass
                        //IMPLEMENT AbsoluteTolerance expandedFragmentTolerance = new AbsoluteTolerance(Math.Sqrt(Math.Pow(CommonParameters.ProductMassTolerance.Value, 2) + Math.Pow(thePrecursorMass / 1000000 * precursorTolerance.Value, 2)));
                        previousTheAdd = 1 + complementaryIons[0].intensity / spectrum.TotalIonCurrent;
                        //we already subtracted that proton, so don't add it again (unit test should break if you do!)
                        previousExperimentalPeakInDaltons = complementaryIons[0].mass;
                        //need to use original tolerance since it's mass based.
                        double previousOriginalMassInDaltons = massShiftForComplementaryConversion - previousExperimentalPeakInDaltons;
                        previousMinRange = previousExperimentalPeakInDaltons - previousOriginalMassInDaltons + CommonParameters.ProductMassTolerance.GetMinimumValue(previousOriginalMassInDaltons);
                        for (int i = 1; i < complementaryIons.Length; i++)
                        {
                            //we already subtracted that proton when making comp ions, so don't add it again (unit test should break if you do!)
                            double experimentalPeakInDaltons = complementaryIons[i].mass;
                            double originalMassInDaltons = massShiftForComplementaryConversion - experimentalPeakInDaltons;
                            if (CommonParameters.ProductMassTolerance.Within(previousOriginalMassInDaltons, originalMassInDaltons))
                            {
                                previousTheAdd += complementaryIons[i].intensity / spectrum.TotalIonCurrent; //open to debate, currently sum intensities of all peaks within tolerance like it was low res. Classic search takes first intensity.
                            }
                            else
                            {
                                previousMaxRange = previousExperimentalPeakInDaltons - previousOriginalMassInDaltons + CommonParameters.ProductMassTolerance.GetMaximumValue(previousOriginalMassInDaltons);
                                FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);
                                previousTheAdd = 1 + complementaryIons[i].intensity / spectrum.TotalIonCurrent;
                                previousMinRange = experimentalPeakInDaltons - originalMassInDaltons + CommonParameters.ProductMassTolerance.GetMinimumValue(originalMassInDaltons);
                            }
                            previousExperimentalPeakInDaltons = experimentalPeakInDaltons;
                            previousOriginalMassInDaltons = massShiftForComplementaryConversion - previousExperimentalPeakInDaltons;
                        }
                        previousMaxRange = previousExperimentalPeakInDaltons - previousOriginalMassInDaltons + CommonParameters.ProductMassTolerance.GetMaximumValue(previousOriginalMassInDaltons);
                        FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
            }
        }

        protected void FindPeakMatches(double theAdd, double min, double max, double[] peptideScores)
        {
            float closestPeak;
            int ipos = Array.BinarySearch(keys, (float)min);
            if (ipos < 0)
                ipos = ~ipos;

            while (ipos < keys.Length)
            {
                closestPeak = keys[ipos];
                if (closestPeak < max)
                {
                    foreach (int heh in fragmentIndex[ipos])
                        peptideScores[heh] += theAdd;
                }
                else
                    break;
                ipos++;
            }

        }

        #endregion Protected Methods
    }
}