using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using IO.ThermoRawFileReader;
using Proteomics.AminoAcidPolymer;

namespace EngineLayer
{
    public class IsotopeAnalysis
    {
        public Dictionary<string, MsDataScan[]> _MS1Scans { get; private set; }

        public List<PeptideSpectralMatch> AllPsms { get; set; }

        public Dictionary<string, List<(double, double)>> _modifiedSequenceToIsotopicDistribution {get; set; }

        public static int MinimumNumberIsotopesRequired = 2;

        MassDiffAcceptor MassDiffAcceptor; 

        private void CalculateTheoreticalIsotopeDistributions()
        {
            _modifiedSequenceToIsotopicDistribution = new Dictionary<string, List<(double, double)>>();

            // calculate averagine (used for isotopic distributions for unknown modifications)
            double averageC = 4.9384;
            double averageH = 7.7583;
            double averageO = 1.4773;
            double averageN = 1.3577;
            double averageS = 0.0417;

            double averagineMass =
                PeriodicTable.GetElement("C").AverageMass * averageC +
                PeriodicTable.GetElement("H").AverageMass * averageH +
                PeriodicTable.GetElement("O").AverageMass * averageO +
                PeriodicTable.GetElement("N").AverageMass * averageN +
                PeriodicTable.GetElement("S").AverageMass * averageS;

            // calculate monoisotopic masses and isotopic envelope for the base sequences
            foreach (PeptideSpectralMatch psm in AllPsms)
            {
                if (_modifiedSequenceToIsotopicDistribution.ContainsKey(psm.FullSequence))
                {
                    continue;
                }

                if (psm.PeptideMonisotopicMass == null) continue; //This isn't the best way of handling a null value TODO: Fix this

                ChemicalFormula formula = psm.ModsChemicalFormula;

                var isotopicMassesAndNormalizedAbundances = new List<(double massShift, double abundance)>();

                if (formula == null)
                {
                    double massDiff = (double)psm.PeptideMonisotopicMass;

                    if (!String.IsNullOrEmpty(psm.BaseSequence))
                    {
                        Peptide baseSequence = new Peptide(psm.BaseSequence);
                        formula = baseSequence.GetChemicalFormula();
                        // add averagine for any unknown mass difference (i.e., a modification)
                        massDiff -= baseSequence.MonoisotopicMass;

                        if (Math.Abs(massDiff) > 20)
                        {
                            double averagines = massDiff / averagineMass;

                            formula.Add("C", (int)Math.Round(averagines * averageC, 0));
                            formula.Add("H", (int)Math.Round(averagines * averageH, 0));
                            formula.Add("O", (int)Math.Round(averagines * averageO, 0));
                            formula.Add("N", (int)Math.Round(averagines * averageN, 0));
                            formula.Add("S", (int)Math.Round(averagines * averageS, 0));
                        }
                    }
                    else
                    {
                        double averagines = massDiff / averagineMass;
                        string averagineFormulaString =
                            "C" + (int)Math.Round(averagines * averageC, 0) +
                            "H" + (int)Math.Round(averagines * averageH, 0) +
                            "O" + (int)Math.Round(averagines * averageO, 0) +
                            "N" + (int)Math.Round(averagines * averageN, 0) +
                            "S" + (int)Math.Round(averagines * averageS, 0);
                        formula = ChemicalFormula.ParseFormula(averagineFormulaString);
                    }
                }

                var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);

                double[] masses = isotopicDistribution.Masses.ToArray();
                double[] abundances = isotopicDistribution.Intensities.ToArray();

                for (int i = 0; i < masses.Length; i++)
                {
                    masses[i] += ((double)psm.PeptideMonisotopicMass - formula.MonoisotopicMass);
                }

                double highestAbundance = abundances.Max();
                int highestAbundanceIndex = Array.IndexOf(abundances, highestAbundance);

                for (int i = 0; i < masses.Length; i++)
                {
                    // expected isotopic mass shifts for this peptide
                    masses[i] -= (double)psm.PeptideMonisotopicMass;

                    // normalized abundance of each isotope
                    abundances[i] /= highestAbundance;

                    // look for these isotopes
                    if (isotopicMassesAndNormalizedAbundances.Count < MinimumNumberIsotopesRequired || abundances[i] > 0.1)
                    {
                        isotopicMassesAndNormalizedAbundances.Add((masses[i], abundances[i]));
                    }
                }

                _modifiedSequenceToIsotopicDistribution.Add(psm.FullSequence, isotopicMassesAndNormalizedAbundances);
            }

            var peptideModifiedSequences = AllPsms.GroupBy(p => p.FullSequence);
            foreach (var fullSequenceGroup in peptideModifiedSequences)
            {
                // isotope where normalized abundance is 1
                double mostAbundantIsotopeShift = _modifiedSequenceToIsotopicDistribution[fullSequenceGroup.First().FullSequence].
                    First(p => p.Item2 == 1.0).Item1;

                foreach (PeptideSpectralMatch psm in fullSequenceGroup)
                {
                    psm.PeakFindingMass = (double)psm.PeptideMonisotopicMass + mostAbundantIsotopeShift;
                }
            }
        }

        public void FindAndScoreEnvelopes()
        {
            var psmsByFile = AllPsms.OrderBy(p => p.ScanNumber).GroupBy(p => p.FullFilePath);

            foreach (var psmGroup in psmsByFile)
            {
                IEnumerable<int> precursorScanIndices = psmGroup.Where(p => p.PrecursorScanNumber != null).
                    Select(p => (int)p.PrecursorScanNumber).Distinct();
                Queue<int> precursorScanQueue = new Queue<int>(precursorScanIndices);
                Queue<MsDataScan> ms1Scans = null;

                for (int i = 0; i < _MS1Scans[psmGroup.Key].Length; i++ )
                {
                    if (_MS1Scans[psmGroup.Key][i] != null && _MS1Scans[psmGroup.Key][i].OneBasedScanNumber == precursorScanQueue.Peek())
                    {
                        ms1Scans.Enqueue(_MS1Scans[psmGroup.Key][i]);
                        precursorScanQueue.Dequeue();
                    }
                }

                foreach (PeptideSpectralMatch psm in psmGroup)
                {
                    IsotopicEnvelope envelope = FindEnvelope(psm, ms1Scans.Dequeue());
                    psm.MS1Envelope = envelope;
                }
                
            }
        }

        public void ScoreEnvelope(IsotopicEnvelope envelope, PeptideSpectralMatch psm)
        {
            //PpmTolerance isotopeTolerance = new PpmTolerance(IsotopePpmTolerance);
            var isotopeMassShifts = _modifiedSequenceToIsotopicDistribution[psm.FullSequence];

            double[] experimentalIsotopeIntensities = new double[isotopeMassShifts.Count];
            double[] theoreticalIsotopeMassShifts = isotopeMassShifts.Select(p => p.Item1).ToArray();
            double[] theoreticalIsotopeAbundances = isotopeMassShifts.Select(p => p.Item2).ToArray();

            List<int> directions = new List<int> { -1, 1 };

            var massShiftToIsotopePeaks = new Dictionary<int, List<(double expIntensity, double theorIntensity, double theorMass)>>
            {
                { -1, new List<(double, double, double)>() },
                { 0, new List<(double, double, double)>() },
                { 1, new List<(double, double, double)>() },
            };

            foreach (var shift in massShiftToIsotopePeaks)
            {
                // look for each isotope peak in the data
                foreach (int direction in directions)
                {
                    int start = (direction == -1) ? peakfindingMassIndex - 1 : peakfindingMassIndex;

                    for (int i = start; i < theoreticalIsotopeAbundances.Length && i >= 0; i += direction)
                    {
                        double isotopeMass = identification.MonoisotopicMass + observedMassError + theoreticalIsotopeMassShifts[i] + shift.Key * Constants.C13MinusC12;
                        double theoreticalIsotopeIntensity = theoreticalIsotopeAbundances[i] * peak.Intensity;

                        IndexedMassSpectralPeak isotopePeak = _peakIndexingEngine.GetIndexedPeak(isotopeMass,
                            peak.ZeroBasedMs1ScanIndex, isotopeTolerance, chargeState);

                        if (isotopePeak == null
                            || isotopePeak.Intensity < theoreticalIsotopeIntensity / 4.0 || isotopePeak.Intensity > theoreticalIsotopeIntensity * 4.0)
                        {
                            break;
                        }

                        shift.Value.Add((isotopePeak.Intensity, theoreticalIsotopeIntensity, isotopeMass));
                        if (shift.Key == 0)
                        {
                            experimentalIsotopeIntensities[i] = isotopePeak.Intensity;
                        }
                    }
                }
            }
        }

        public List<IsotopicEnvelope> GetIsotopicEnvelopes(List<IndexedMassSpectralPeak> xic, Identification identification, int chargeState)
        {
            var isotopicEnvelopes = new List<IsotopicEnvelope>();
            var isotopeMassShifts = _modifiedSequenceToIsotopicDistribution[identification.ModifiedSequence];

            if (isotopeMassShifts.Count < NumIsotopesRequired)
            {
                return isotopicEnvelopes;
            }

            PpmTolerance isotopeTolerance = new PpmTolerance(IsotopePpmTolerance);

            double[] experimentalIsotopeIntensities = new double[isotopeMassShifts.Count];
            double[] theoreticalIsotopeMassShifts = isotopeMassShifts.Select(p => p.Item1).ToArray();
            double[] theoreticalIsotopeAbundances = isotopeMassShifts.Select(p => p.Item2).ToArray();
            int peakfindingMassIndex = (int)Math.Round(identification.PeakfindingMass - identification.MonoisotopicMass, 0);

            List<int> directions = new List<int> { -1, 1 };

            var massShiftToIsotopePeaks = new Dictionary<int, List<(double expIntensity, double theorIntensity, double theorMass)>>
            {
                { -1, new List<(double, double, double)>() },
                { 0, new List<(double, double, double)>() },
                { 1, new List<(double, double, double)>() },
            };

            foreach (IndexedMassSpectralPeak peak in xic)
            {
                Array.Clear(experimentalIsotopeIntensities, 0, experimentalIsotopeIntensities.Length);
                foreach (var kvp in massShiftToIsotopePeaks)
                {
                    kvp.Value.Clear();
                }

                // isotope masses are calculated relative to the observed peak
                double observedMass = peak.Mz.ToMass(chargeState);
                double observedMassError = observedMass - identification.PeakfindingMass;

                foreach (var shift in massShiftToIsotopePeaks)
                {
                    // look for each isotope peak in the data
                    foreach (int direction in directions)
                    {
                        int start = direction == -1 ? peakfindingMassIndex - 1 : peakfindingMassIndex;

                        for (int i = start; i < theoreticalIsotopeAbundances.Length && i >= 0; i += direction)
                        {
                            double isotopeMass = identification.MonoisotopicMass + observedMassError + theoreticalIsotopeMassShifts[i] + shift.Key * Constants.C13MinusC12;
                            double theoreticalIsotopeIntensity = theoreticalIsotopeAbundances[i] * peak.Intensity;

                            IndexedMassSpectralPeak isotopePeak = _peakIndexingEngine.GetIndexedPeak(isotopeMass,
                                peak.ZeroBasedMs1ScanIndex, isotopeTolerance, chargeState);

                            if (isotopePeak == null
                                || isotopePeak.Intensity < theoreticalIsotopeIntensity / 4.0 || isotopePeak.Intensity > theoreticalIsotopeIntensity * 4.0)
                            {
                                break;
                            }

                            shift.Value.Add((isotopePeak.Intensity, theoreticalIsotopeIntensity, isotopeMass));
                            if (shift.Key == 0)
                            {
                                experimentalIsotopeIntensities[i] = isotopePeak.Intensity;
                            }
                        }
                    }
                }

                // check number of isotope peaks observed
                if (massShiftToIsotopePeaks[0].Count < NumIsotopesRequired)
                {
                    continue;
                }

                double corr = Correlation.Pearson(massShiftToIsotopePeaks[0].Select(p => p.expIntensity), massShiftToIsotopePeaks[0].Select(p => p.theorIntensity));

                // check correlation of experimental isotope intensities to the theoretical abundances
                foreach (var shift in massShiftToIsotopePeaks)
                {
                    if (!shift.Value.Any())
                    {
                        continue;
                    }

                    double unexpectedMass = shift.Value.Min(p => p.theorMass) - Constants.C13MinusC12;

                    IndexedMassSpectralPeak unexpectedPeak = _peakIndexingEngine.GetIndexedPeak(unexpectedMass,
                                peak.ZeroBasedMs1ScanIndex, isotopeTolerance, chargeState);

                    if (unexpectedPeak == null)
                    {
                        shift.Value.Add((0, 0, unexpectedMass));
                    }
                    else
                    {
                        shift.Value.Add((unexpectedPeak.Intensity, 0, unexpectedMass));
                    }
                }

                double corrWithPadding = Correlation.Pearson(massShiftToIsotopePeaks[0].Select(p => p.expIntensity), massShiftToIsotopePeaks[0].Select(p => p.theorIntensity));
                double corrShiftedLeft = Correlation.Pearson(massShiftToIsotopePeaks[-1].Select(p => p.expIntensity), massShiftToIsotopePeaks[-1].Select(p => p.theorIntensity));
                double corrShiftedRight = Correlation.Pearson(massShiftToIsotopePeaks[1].Select(p => p.expIntensity), massShiftToIsotopePeaks[1].Select(p => p.theorIntensity));

                if (double.IsNaN(corrShiftedLeft))
                {
                    corrShiftedLeft = -1;
                }
                if (double.IsNaN(corrShiftedRight))
                {
                    corrShiftedRight = -1;
                }

                if (corr > 0.7 && (corrShiftedLeft - corrWithPadding < 0.1 && corrShiftedRight - corrWithPadding < 0.1))
                {
                    // impute unobserved isotope peak intensities
                    for (int i = 0; i < experimentalIsotopeIntensities.Length; i++)
                    {
                        if (experimentalIsotopeIntensities[i] == 0)
                        {
                            experimentalIsotopeIntensities[i] = theoreticalIsotopeAbundances[i] * experimentalIsotopeIntensities[peakfindingMassIndex];
                        }
                    }

                    isotopicEnvelopes.Add(new IsotopicEnvelope(peak, chargeState, experimentalIsotopeIntensities.Sum()));
                }
            }

            return isotopicEnvelopes;
        }

        public IsotopicEnvelope FindEnvelope(PeptideSpectralMatch psm, MsDataScan scan)
        {
            DoubleRange monoMass = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass((double)psm.PeptideMonisotopicMass).
                Select(i => i.AllowedInterval).First();
            DoubleRange mostAbundantMass = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(psm.PeakFindingMass).
                Select(i => i.AllowedInterval).First();
            double experimentalMonoMass = 0;
            double experimentalMostAbundantMass = 0;
            MzSpectrum spectrum = scan.MassSpectrum;
            Range? spectrumRange = null;
            int envelopeSize = _modifiedSequenceToIsotopicDistribution[psm.FullSequence].Count;
            for (int i = 0; i < spectrum.XArray.Length; i++)
            {
                if (monoMass.Contains(spectrum.XArray[i]))
                {
                    // For the envelope, we want to look at the peak immediately before the monoisotopic mass,
                    // so as to avoid situations where the "monoisotopic peak" is actually an isotope peak for something else 
                    int rangeStart = i > 0 ? i - 1 : i;
                    // Look at peak after the end of the envelope as well
                    spectrumRange = rangeStart + 1 + envelopeSize < spectrum.XArray.Length ? 
                        new Range(rangeStart, rangeStart + envelopeSize + 1) :
                        new Range(rangeStart, spectrum.XArray.Length);
                    experimentalMonoMass = spectrum.XArray[i];
                }
                if (mostAbundantMass.Contains(spectrum.XArray[i]))
                {
                    experimentalMostAbundantMass = spectrum.XArray[i];
                    if (spectrumRange == null)
                    {
                        int massShiftMonoToAbundant = (int)Math.Round(mostAbundantMass.Mean - monoMass.Mean, 0);
                        // want to see what's before the most abundant, in case mono was slightly outside of tolerance or w/e
                        int rangeStart = i > (massShiftMonoToAbundant + 1) ? i - massShiftMonoToAbundant - 1 : 0;
                        spectrumRange = rangeStart + 1 + envelopeSize < spectrum.XArray.Length ?
                            new Range(rangeStart, rangeStart + envelopeSize + 1) :
                            new Range(rangeStart, spectrum.XArray.Length);
                    }
                    break;
                }
            }

            if (spectrumRange == null) return null;
            IsotopicEnvelope envelope = new(
                MakeTupleList(spectrum.XArray.Take((Range)spectrumRange).ToArray(), spectrum.YArray.Take((Range)spectrumRange).ToArray()),
                experimentalMonoMass, psm.ScanPrecursorCharge, 0, 0, 0);
            if (Math.Abs(envelope.MostAbundantObservedIsotopicMass - experimentalMostAbundantMass) > 0.01) throw new Exception("idk"); //Should change this at some point, but I'm interested to see if it breaks in testing
            return envelope;
        }

        public List< (double, double) > MakeTupleList(double[] mzArray, double[] intensityArray)
        {
            if (mzArray.Length != intensityArray.Length) return null;
            List<(double, double)> tupleList = new();
            for (int i = 0; i < mzArray.Length; i++)
            {
                tupleList.Add((mzArray[i], intensityArray[i]));
            }
            return tupleList;
        }

        public bool ReadMS1Scans(string filePath, bool silent)
        {
            if (!silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            MsDataScan[] msDataScans = null;

            // read spectra file
            string ext = Path.GetExtension(filePath).ToUpperInvariant();
            if (ext.Equals(".MZML"))
            {
                try
                {
                    msDataScans = Mzml.LoadAllStaticData(filePath).GetAllScansList()
                        .OrderBy(p => p.OneBasedScanNumber).ToArray();
                }
                catch (FileNotFoundException)
                {
                    if (!silent)
                    {
                        Console.WriteLine("\nCan't find .mzML file" + filePath + "\n");
                    }

                    return false;
                }
                catch (Exception e)
                {
                    if (!silent)
                    {
                        Console.WriteLine("Problem opening .mzML file " + filePath + "; " +
                                          e.Message);
                    }

                    return false;
                }

                for (int i = 0; i < msDataScans.Length; i++)
                {
                    if (msDataScans[i].MsnOrder > 1)
                    {
                        msDataScans[i] = null;
                    }
                }
            }
            else if (ext.Equals(".RAW"))
            {
                var tempList = new List<MsDataScan>();
                ThermoDynamicData dynamicConnection = null;

                try
                {
                    dynamicConnection = new ThermoDynamicData(filePath);

                    // use thermo dynamic connection to get the ms1 scans and then dispose of the connection
                    for (int i = 0; i < dynamicConnection.MsOrdersByScan.Length; i++)
                    {
                        if (dynamicConnection.MsOrdersByScan[i] == 1)
                        {
                            tempList.Add(dynamicConnection.GetOneBasedScanFromDynamicConnection(i + 1));
                        }
                        else
                        {
                            tempList.Add(null);
                        }
                    }

                    dynamicConnection.CloseDynamicConnection();
                }
                catch (FileNotFoundException)
                {
                    if (dynamicConnection != null)
                    {
                        dynamicConnection.CloseDynamicConnection();
                    }

                    if (!silent)
                    {
                        Console.WriteLine("\nCan't find .raw file" + filePath + "\n");
                    }

                    return false;
                }
                catch (Exception e)
                {
                    if (dynamicConnection != null)
                    {
                        dynamicConnection.CloseDynamicConnection();
                    }

                    if (!silent)
                    {
                        throw new MzLibException("FlashLFQ Error: Problem opening .raw file " + filePath + "; " + e.Message);
                    }
                }

                msDataScans = tempList.ToArray();
            }
            else
            {
                if (!silent)
                {
                    Console.WriteLine("Unsupported file type " + ext);
                    return false;
                }
            }

            _MS1Scans.Add(filePath, msDataScans);
            return true;
        }
    }
}
