using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace IndexSearchAndAnalyze
{
    internal class TrainingPointsExtractor
    {
        private const int numFragmentsNeeded = 10;

        public static List<LabeledDataPoint> GetDataPoints(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, List<NewPsmWithFDR> identifications, SoftwareLockMassParams p, HashSet<int> matchesToExclude)
        {
            p.po.output("Extracting data points:");
            // The final training point list
            List<LabeledDataPoint> trainingPointsToReturn = new List<LabeledDataPoint>();

            // Set of peaks, identified by m/z and retention time. If a peak is in here, it means it has been a part of an accepted identification, and should be rejected
            HashSet<Tuple<double, double>> peaksAddedFromMS1HashSet = new HashSet<Tuple<double, double>>();

            int numIdentifications = identifications.Count;
            // Loop over identifications
            for (int matchIndex = 0; matchIndex < numIdentifications; matchIndex++)
            {
                var identification = identifications[matchIndex];
                if (identification.QValue > 0.01)
                    break;

                // Progress
                if (numIdentifications < 100 || matchIndex % (numIdentifications / 100) == 0)
                    p.po.ReportProgress(new ProgressEventArgs(100 * matchIndex / numIdentifications, "Looking at identifications..."));

                // Skip decoys, they are for sure not there!
                if (identification.isDecoy)
                    continue;

                // Skip exclusions! These are for testing
                if (matchesToExclude.Contains(matchIndex))
                    continue;

                // Each identification has an MS2 spectrum attached to it.
                int ms2spectrumIndex = identification.thisPSM.newPsm.scanNumber;

                // Get the peptide, don't forget to add the modifications!!!!
                var SequenceWithChemicalFormulas = identification.thisPSM.SequenceWithChemicalFormulas;
                int peptideCharge = identification.thisPSM.newPsm.scanPrecursorCharge;

                #region watch

                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    p.po.output("ms2spectrumIndex: " + ms2spectrumIndex);
                    p.po.output(" peptide: " + SequenceWithChemicalFormulas);
                }

                #endregion watch

                int numFragmentsIdentified = -1;
                List<LabeledDataPoint> candidateTrainingPointsForPeptide = new List<LabeledDataPoint>();
                Proteomics.Peptide coolPeptide = null;
                try
                {
                    coolPeptide = new Peptide(SequenceWithChemicalFormulas);
                }
                catch
                {
                    continue;
                }

                candidateTrainingPointsForPeptide = SearchMS2Spectrum(myMsDataFile.GetOneBasedScan(ms2spectrumIndex), coolPeptide, peptideCharge, p, out numFragmentsIdentified);

                //SoftwareLockMassRunner.WriteDataToFiles(candidateTrainingPointsForPeptide, ms2spectrumIndex.ToString());

                //p.po.RTBoutput(new OutputHandlerEventArgs(numFragmentsIdentified.ToString()));

                // If MS2 has low evidence for peptide, skip and go to next one
                if (numFragmentsIdentified < numFragmentsNeeded)
                    continue;

                // Calculate isotopic distribution of the full peptide

                IsotopicDistribution dist = new IsotopicDistribution(coolPeptide.GetChemicalFormula(), p.fineResolution, 0.001);

                double[] masses = new double[dist.Masses.Count];
                double[] intensities = new double[dist.Intensities.Count];
                for (int i = 0; i < dist.Masses.Count; i++)
                {
                    masses[i] = dist.Masses[i];
                    intensities[i] = dist.Intensities[i];
                }
                Array.Sort(intensities, masses, Comparer<double>.Create((x, y) => y.CompareTo(x)));

                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    p.po.output(" Isotopologue distribution: ");
                    p.po.output(" masses = " + string.Join(", ", masses) + "...");
                    p.po.output(" intensities = " + string.Join(", ", intensities) + "...");
                }

                int lowestMS1ind = SearchMS1Spectra(myMsDataFile, masses, intensities, candidateTrainingPointsForPeptide, ms2spectrumIndex, -1, peaksAddedFromMS1HashSet, p, peptideCharge);

                int highestMS1ind = SearchMS1Spectra(myMsDataFile, masses, intensities, candidateTrainingPointsForPeptide, ms2spectrumIndex, 1, peaksAddedFromMS1HashSet, p, peptideCharge);

                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    p.po.output(" ms1range: " + lowestMS1ind + " to " + highestMS1ind);
                }
                trainingPointsToReturn.AddRange(candidateTrainingPointsForPeptide);
            }
            p.po.output("Number of training points: " + trainingPointsToReturn.Count());
            return trainingPointsToReturn;
        }

        private static int GetOneBasedSpectrumIndexFromID(string spectrumID, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile)
        {
            var lastInt = Convert.ToInt32(Regex.Match(spectrumID, @"\d+$").Value);
            if (lastInt > 0 && lastInt <= myMsDataFile.Count() && myMsDataFile.GetOneBasedScan(lastInt).id.Equals(spectrumID))
                return lastInt;
            if (lastInt - 1 > 0 && lastInt - 1 <= myMsDataFile.Count() && myMsDataFile.GetOneBasedScan(lastInt - 1).id.Equals(spectrumID))
                return lastInt - 1;
            if (lastInt + 1 > 0 && lastInt + 1 <= myMsDataFile.Count() && myMsDataFile.GetOneBasedScan(lastInt + 1).id.Equals(spectrumID))
                return lastInt + 1;

            int oneBasedScanNumber = 0;
            do
            {
                oneBasedScanNumber++;
            } while (!spectrumID.Equals(myMsDataFile.GetOneBasedScan(oneBasedScanNumber).id));
            return oneBasedScanNumber;
        }

        private static List<LabeledDataPoint> SearchMS2Spectrum(IMsDataScan<IMzSpectrum<MzPeak>> ms2DataScan, Proteomics.Peptide peptide, int peptideCharge, SoftwareLockMassParams p, out int candidateFragmentsIdentified)
        {
            List<LabeledDataPoint> myCandidatePoints = new List<LabeledDataPoint>();

            // Key: mz value, Value: error
            Dictionary<double, double> addedPeaks = new Dictionary<double, double>();

            int SelectedIonGuessChargeStateGuess;
            ms2DataScan.TryGetSelectedIonGuessChargeStateGuess(out SelectedIonGuessChargeStateGuess);
            double IsolationMZ;
            ms2DataScan.TryGetIsolationMZ(out IsolationMZ);

            int ms2spectrumIndex = ms2DataScan.OneBasedScanNumber;

            var countForThisMS2 = 0;
            var countForThisMS2a = 0;
            var numFragmentsIdentified = 0;

            var scanWindowRange = ms2DataScan.ScanWindowRange;

            Fragment[] fragmentList = peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true).ToArray();

            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
            {
                p.po.output(" Considering individual fragments:");
            }

            foreach (IHasChemicalFormula fragment in fragmentList)
            {
                bool fragmentIdentified = false;
                bool computedIsotopologues = false;
                double[] masses = new double[0];
                double[] intensities = new double[0];
                // First look for monoisotopic masses, do not compute distribution spectrum!
                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    p.po.output("  Considering fragment " + (fragment as Fragment).Sequence + " with formula " + fragment.ThisChemicalFormula.Formula);
                    //if ((fragment as Fragment).Modifications.Count() > 0)
                    //p.po.RTBoutput("  Modifications: " + string.Join(", ", (fragment as Fragment).Modifications));
                }

                #region loop to determine if need to compute isotopologue distribution

                for (int chargeToLookAt = 1; chargeToLookAt <= peptideCharge; chargeToLookAt++)
                {
                    var monoisotopicMZ = fragment.MonoisotopicMass.ToMassToChargeRatio(chargeToLookAt);
                    if (monoisotopicMZ > scanWindowRange.Maximum)
                        continue;
                    if (monoisotopicMZ < scanWindowRange.Minimum)
                        break;
                    var closestPeakMZ = ms2DataScan.MassSpectrum.GetClosestPeakXvalue(monoisotopicMZ);
                    if (Math.Abs(closestPeakMZ - monoisotopicMZ) < p.toleranceInMZforMS2Search)
                    {
                        if (!computedIsotopologues)
                        {
                            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                p.po.output("    Computing isotopologues because absolute error " + Math.Abs(closestPeakMZ - monoisotopicMZ) + " is smaller than tolerance " + p.toleranceInMZforMS2Search);
                                p.po.output("    Charge was = " + chargeToLookAt + "  closestPeakMZ = " + closestPeakMZ + " while monoisotopicMZ = " + monoisotopicMZ);
                            }

                            IsotopicDistribution dist = new IsotopicDistribution(fragment.ThisChemicalFormula, p.fineResolution, 0.001);

                            masses = new double[dist.Masses.Count];
                            intensities = new double[dist.Intensities.Count];
                            for (int i = 0; i < dist.Masses.Count; i++)
                            {
                                masses[i] = dist.Masses[i];
                                intensities[i] = dist.Intensities[i];
                            }
                            Array.Sort(intensities, masses, Comparer<double>.Create((x, y) => y.CompareTo(x)));
                            computedIsotopologues = true;
                            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                p.po.output("    Isotopologue distribution: ");
                                p.po.output("    masses = " + string.Join(", ", masses) + "...");
                                p.po.output("    intensities = " + string.Join(", ", intensities) + "...");
                            }

                            break;
                        }
                    }
                }

                #endregion loop to determine if need to compute isotopologue distribution

                if (computedIsotopologues)
                {
                    #region actually add training points

                    if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                    {
                        p.po.output("   Considering individual charges, to get training points:");
                    }
                    bool startingToAdd = false;
                    for (int chargeToLookAt = 1; chargeToLookAt <= peptideCharge; chargeToLookAt++)
                    {
                        if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                        {
                            p.po.output("    Considering charge " + chargeToLookAt);
                        }
                        if (masses.First().ToMassToChargeRatio(chargeToLookAt) > scanWindowRange.Maximum)
                        {
                            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                p.po.output("    Out of range: too high");
                            }
                            continue;
                        }
                        if (masses.Last().ToMassToChargeRatio(chargeToLookAt) < scanWindowRange.Minimum)
                        {
                            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                p.po.output("    Out of range: too low");
                            }
                            break;
                        }
                        List<TrainingPoint> trainingPointsToAverage = new List<TrainingPoint>();
                        foreach (double a in masses)
                        {
                            double theMZ = a.ToMassToChargeRatio(chargeToLookAt);
                            var npwr = ms2DataScan.MassSpectrum.NumPeaksWithinRange(theMZ - p.toleranceInMZforMS2Search, theMZ + p.toleranceInMZforMS2Search);
                            if (npwr == 0)
                            {
                                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                                    p.po.output("     Breaking because extracted.Count = " + npwr);
                                break;
                            }
                            if (npwr > 1)
                            {
                                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                                    p.po.output("     Not looking for " + theMZ + " because extracted.Count = " + npwr);
                                continue;
                            }
                            var closestPeak = ms2DataScan.MassSpectrum.GetClosestPeak(theMZ);
                            var closestPeakMZ = closestPeak.MZ;
                            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                p.po.output("     Found       " + closestPeakMZ + "   Error is    " + (closestPeakMZ - theMZ));
                            }
                            if (!addedPeaks.ContainsKey(closestPeakMZ))
                            {
                                addedPeaks.Add(closestPeakMZ, Math.Abs(closestPeakMZ - theMZ));
                                trainingPointsToAverage.Add(new TrainingPoint(new DataPoint(closestPeakMZ, double.NaN, 0, closestPeak.Intensity, double.NaN, double.NaN), closestPeakMZ - theMZ));
                            }
                            else
                            {
                                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                                {
                                    p.po.output("     Not using because already added peak");
                                }
                            }
                        }
                        // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                        if (trainingPointsToAverage.Count == 0 && startingToAdd == true)
                            break;
                        if (trainingPointsToAverage.Count < Math.Min(p.minMS2, intensities.Count()))
                        {
                            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                p.po.output("    Not adding, since not enought isotopes were seen");
                            }
                        }
                        else
                        {
                            startingToAdd = true;
                            if (!fragmentIdentified)
                            {
                                fragmentIdentified = true;
                                numFragmentsIdentified += 1;
                            }

                            countForThisMS2 += trainingPointsToAverage.Count;
                            countForThisMS2a += 1;

                            double addedMZ = trainingPointsToAverage.Select(b => b.dp.mz).Average();
                            double relativeMZ = (addedMZ - ms2DataScan.ScanWindowRange.Minimum) / (ms2DataScan.ScanWindowRange.Maximum - ms2DataScan.ScanWindowRange.Minimum);
                            double[] inputs = new double[9] { 2, addedMZ, ms2DataScan.RetentionTime, trainingPointsToAverage.Select(b => b.dp.intensity).Average(), ms2DataScan.TotalIonCurrent, ms2DataScan.InjectionTime, SelectedIonGuessChargeStateGuess, IsolationMZ, relativeMZ };
                            var a = new LabeledDataPoint(inputs, trainingPointsToAverage.Select(b => b.l).Median());

                            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                p.po.output("    Adding aggregate of " + trainingPointsToAverage.Count + " points FROM MS2 SPECTRUM");
                                p.po.output("    a.dp.mz " + a.inputs[1]);
                                p.po.output("    a.dp.rt " + a.inputs[2]);
                                p.po.output("    a.l     " + a.output);
                            }
                            myCandidatePoints.Add(a);
                        }
                    }

                    #endregion actually add training points
                }
            }

            //p.po.RTBoutput(new OutputHandlerEventArgs("ind = " + ms2spectrumIndex + " count = " + countForThisMS2 + " count2 = " + countForThisMS2a + " fragments = " + numFragmentsIdentified));
            if (p.MS2spectraToWatch.Contains(ms2spectrumIndex))
            {
                p.po.output(" countForThisMS2 = " + countForThisMS2);
                p.po.output(" countForThisMS2a = " + countForThisMS2a);
                p.po.output(" numFragmentsIdentified = " + numFragmentsIdentified);
            }
            candidateFragmentsIdentified = numFragmentsIdentified;
            return myCandidatePoints;
        }

        private static int SearchMS1Spectra(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, double[] originalMasses, double[] originalIntensities, List<LabeledDataPoint> myCandidatePoints, int ms2spectrumIndex, int direction, HashSet<Tuple<double, double>> peaksAddedHashSet, SoftwareLockMassParams p, int peptideCharge)
        {
            int goodIndex = -1;
            List<int> scores = new List<int>();
            var theIndex = -1;
            if (direction == 1)
                theIndex = ms2spectrumIndex;
            else
                theIndex = ms2spectrumIndex - 1;

            bool addedAscan = true;

            int highestKnownChargeForThisPeptide = peptideCharge;
            while (theIndex >= 1 && theIndex <= myMsDataFile.NumSpectra && addedAscan == true)
            {
                int countForThisScan = 0;
                if (myMsDataFile.GetOneBasedScan(theIndex).MsnOrder > 1)
                {
                    theIndex += direction;
                    continue;
                }
                addedAscan = false;
                if (p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex))
                {
                    p.po.output(" Looking in MS1 spectrum " + theIndex + " because of MS2 spectrum " + ms2spectrumIndex);
                }
                List<LabeledDataPoint> myCandidatePointsForThisMS1scan = new List<LabeledDataPoint>();
                var fullMS1scan = myMsDataFile.GetOneBasedScan(theIndex);
                double ms1RetentionTime = fullMS1scan.RetentionTime;
                var scanWindowRange = fullMS1scan.ScanWindowRange;
                var fullMS1spectrum = fullMS1scan.MassSpectrum;
                if (fullMS1spectrum.Count == 0)
                    break;

                bool startingToAddCharges = false;
                int chargeToLookAt = 1;
                do
                {
                    if (p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex))
                    {
                        p.po.output("  Looking at charge " + chargeToLookAt);
                    }
                    if (originalMasses[0].ToMassToChargeRatio(chargeToLookAt) > scanWindowRange.Maximum)
                    {
                        chargeToLookAt++;
                        continue;
                    }
                    if (originalMasses[0].ToMassToChargeRatio(chargeToLookAt) < scanWindowRange.Minimum)
                        break;
                    List<TrainingPoint> trainingPointsToAverage = new List<TrainingPoint>();
                    foreach (double a in originalMasses)
                    {
                        double theMZ = a.ToMassToChargeRatio(chargeToLookAt);

                        var npwr = fullMS1spectrum.NumPeaksWithinRange(theMZ - p.toleranceInMZforMS1Search, theMZ + p.toleranceInMZforMS1Search);
                        if (npwr == 0)
                        {
                            if (p.MS1spectraToWatch.Contains(theIndex))
                                p.po.output("      Breaking because extracted.Count = " + npwr);
                            break;
                        }
                        if (npwr > 1)
                        {
                            if (p.MS1spectraToWatch.Contains(theIndex))
                                p.po.output("      Not looking for " + theMZ + " because extracted.Count = " + npwr);
                            continue;
                        }

                        var closestPeak = fullMS1spectrum.GetClosestPeak(theMZ);
                        var closestPeakMZ = closestPeak.MZ;

                        //if (closestPeak.Intensity / fullMS1scan.TotalIonCurrent < 2e-4)
                        //{
                        //    if (p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex))
                        //        p.po.RTBoutput("      Breaking because intensity fraction is " + closestPeak.Intensity / fullMS1scan.TotalIonCurrent);
                        //    break;
                        //}
                        if ((p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex)) && p.mzRange.Contains(theMZ))
                        {
                            p.po.output("      Looking for " + theMZ + " found " + closestPeakMZ + " error is " + (closestPeakMZ - theMZ));
                        }

                        var theTuple = Tuple.Create(closestPeakMZ, ms1RetentionTime);
                        if (!peaksAddedHashSet.Contains(theTuple))
                        {
                            peaksAddedHashSet.Add(theTuple);
                            highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
                            //if ((p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex)) && p.mzRange.Contains(theMZ))
                            if ((p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex)))
                            {
                                p.po.output("      Found " + closestPeakMZ + ", was looking for " + theMZ + ", e=" + (closestPeakMZ - theMZ) + ", if=" + closestPeak.Intensity / fullMS1scan.TotalIonCurrent);
                            }
                            trainingPointsToAverage.Add(new TrainingPoint(new DataPoint(closestPeakMZ, double.NaN, 1, closestPeak.Intensity, double.NaN, double.NaN), closestPeakMZ - theMZ));
                        }
                        else
                            break;
                    }
                    // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage.Count == 0 && startingToAddCharges == true)
                    {
                        if (p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex))
                        {
                            p.po.output("    Started adding and suddnely stopped, no need to look at higher charges");
                        }
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 0 || (trainingPointsToAverage.Count == 1 && originalIntensities[0] < 0.65)) && (peptideCharge <= chargeToLookAt))
                    {
                        if (p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex))
                        {
                            p.po.output("    Did not find (or found without isotopes) charge " + chargeToLookAt + ", no need to look at higher charges");
                        }
                        break;
                    }
                    if (trainingPointsToAverage.Count == 1 && originalIntensities[0] < 0.65)
                    {
                        if ((p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex)))
                        {
                            p.po.output("    Not adding, since originalIntensities[0] is " + originalIntensities[0] + " which is too low");
                        }
                    }
                    else if (trainingPointsToAverage.Count < Math.Min(p.minMS1, originalIntensities.Count()))
                    {
                        if ((p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex)))
                        {
                            p.po.output("    Not adding, since count " + trainingPointsToAverage.Count + " is too low");
                        }
                    }
                    else
                    {
                        addedAscan = true;
                        startingToAddCharges = true;
                        countForThisScan += 1;
                        double[] inputs = new double[6] { 1, trainingPointsToAverage.Select(b => b.dp.mz).Average(), fullMS1scan.RetentionTime, trainingPointsToAverage.Select(b => b.dp.intensity).Average(), fullMS1scan.TotalIonCurrent, fullMS1scan.InjectionTime };
                        var a = new LabeledDataPoint(inputs, trainingPointsToAverage.Select(b => b.l).Median());
                        //if (a.output > 0)
                        //    p.po.RTBoutput(theIndex + "," + ms2spectrumIndex);
                        if (p.MS2spectraToWatch.Contains(ms2spectrumIndex) || p.MS1spectraToWatch.Contains(theIndex))
                        {
                            p.po.output("    Adding aggregate of " + trainingPointsToAverage.Count + " points FROM MS1 SPECTRUM");
                            p.po.output("    a.dp.mz " + a.inputs[1]);
                            p.po.output("    a.dp.rt " + a.inputs[2]);
                            p.po.output("    a.l     " + a.output);
                        }
                        myCandidatePointsForThisMS1scan.Add(a);
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);

                if (myCandidatePointsForThisMS1scan.Count > 0)
                    goodIndex = theIndex;
                //    SoftwareLockMassRunner.WriteDataToFiles(myCandidatePointsForThisMS1scan, theIndex.ToString());
                myCandidatePoints.AddRange(myCandidatePointsForThisMS1scan);

                scores.Add(countForThisScan);
                theIndex += direction;
            }
            return goodIndex;
        }
    }
}