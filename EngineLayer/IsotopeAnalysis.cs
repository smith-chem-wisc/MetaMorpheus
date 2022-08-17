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
using MathNet.Numerics.Statistics;

namespace EngineLayer
{
    public class IsotopeAnalysis
    {
        public Dictionary<string, MsDataScan[]> _MS1Scans { get; private set; }
        public List<PeptideSpectralMatch> AllPsms { get; set; }
        public Dictionary<string, List<(double, double)>> _modifiedSequenceToIsotopicDistribution {get; set; }
        public MassDiffAcceptor MassDiffAcceptor { get; set; }

        public readonly int MinimumNumberIsotopesRequired = 2;
        
        public IsotopeAnalysis(List<string> fullFilePaths, List<PeptideSpectralMatch> allPsms, MassDiffAcceptor massDiffAcceptor)
        {
            this._MS1Scans = new();
            this.MassDiffAcceptor = massDiffAcceptor;
            this.AllPsms = allPsms.Where(p => p != null).ToList();
            CalculateTheoreticalIsotopeDistributions();
            foreach (string fullPathToFile in fullFilePaths)
            {
                if (ReadMS1Scans(fullPathToFile, true))
                {
                    FindAndScoreEnvelopes(fullPathToFile);
                }
            }
        }
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
                if (psm.FullSequence == null || _modifiedSequenceToIsotopicDistribution.ContainsKey(psm.FullSequence))
                {
                    continue;
                }

                if (psm.PeptideMonisotopicMass == null) continue; //This isn't the best way of handling a null value TODO: Fix this

                ChemicalFormula formula = null;

                var isotopicMassesAndNormalizedAbundances = new List<(double massShift, double abundance)>();

                if (formula == null)
                {
                    double massDiff = (double)psm.ScanPrecursorMass;

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

            //var peptideModifiedSequences = AllPsms.GroupBy(p => p.FullSequence);
            //foreach (var fullSequenceGroup in peptideModifiedSequences)
            //{
            //    // isotope where normalized abundance is 1
            //    double mostAbundantIsotopeShift = _modifiedSequenceToIsotopicDistribution[fullSequenceGroup.First().FullSequence].
            //        First(p => p.Item2 == 1.0).Item1;

            //    foreach (PeptideSpectralMatch psm in fullSequenceGroup)
            //    {
            //        psm.PeakFindingMass = (double)psm.PeptideMonisotopicMass + mostAbundantIsotopeShift;
            //    }
            //}
        }

        public void FindAndScoreEnvelopes(string fullFilePath)
        {
            var psmsByFile = AllPsms.Where(p => Path.GetFileNameWithoutExtension(p.FullFilePath).Equals(Path.GetFileNameWithoutExtension(fullFilePath))).
                OrderBy(p => p.ScanNumber);

            IEnumerable<int> precursorScanIndices = psmsByFile.Where(p => p.PrecursorScanNumber != null).
                Select(p => (int)p.PrecursorScanNumber).Distinct();
            Queue<int> precursorScanQueue = new Queue<int>(precursorScanIndices);
            Queue<MsDataScan> ms1Scans = new();

            for (int i = 0; i < _MS1Scans[fullFilePath].Length; i++ )
            {
                if (_MS1Scans[fullFilePath][i] != null && _MS1Scans[fullFilePath][i].OneBasedScanNumber == precursorScanQueue.Peek())
                {
                    ms1Scans.Enqueue(_MS1Scans[fullFilePath][i]);
                    precursorScanQueue.Dequeue();
                }
            }

            foreach (PeptideSpectralMatch psm in psmsByFile)
            {
                if (psm.FullSequence == null) continue; // Can't/Won't score ambiguous peptides
                bool scanMatch = false;
                while(!scanMatch)
                {
                    if (!ms1Scans.TryPeek(out MsDataScan scan)) return;
                    int scanNumber = scan.OneBasedScanNumber;
                    int? psmScanNumber = psm.PrecursorScanNumber;
                    if (psmScanNumber != scanNumber) ms1Scans.Dequeue();
                    else scanMatch = true;
                }
                
                IsotopicEnvelope envelope = FindEnvelope(psm, ms1Scans.Peek());
                ScoreEnvelope(envelope, psm);
                //psm.MS1Envelope = envelope;
            }
        }

        public void ScoreEnvelope(IsotopicEnvelope envelope, PeptideSpectralMatch psm)
        {
            //PpmTolerance isotopeTolerance = new PpmTolerance(IsotopePpmTolerance);
            var isotopeMassShifts = _modifiedSequenceToIsotopicDistribution[psm.FullSequence];

            double[] experimentalIsotopeIntensities = new double[isotopeMassShifts.Count];
            double[] theoreticalIsotopeMassShifts = isotopeMassShifts.Select(p => p.Item1).ToArray();
            double[] theoreticalIsotopeAbundances = isotopeMassShifts.Select(p => p.Item2).ToArray();

            var massShiftToIsotopePeaks = new Dictionary<int, List<(double expIntensity, double theorIntensity, double theorMass)>>
            {
                { -1, new List<(double, double, double)>() },
                { 0, new List<(double, double, double)>() },
                { 1, new List<(double, double, double)>() },
            };

            // isotope masses are calculated relative to the observed peak
            double observedMass = psm.ScanPrecursorMass;
            double observedMassError = observedMass - (double)psm.PeptideMonisotopicMass;

            foreach (var shift in massShiftToIsotopePeaks)
            {
                for (int i = 0; i < theoreticalIsotopeAbundances.Length; i++)
                {
                    double isotopeMass = (double)psm.PeptideMonisotopicMass + observedMassError + theoreticalIsotopeMassShifts[i] + shift.Key * Constants.C13MinusC12;
                    double theoreticalIsotopeIntensity = theoreticalIsotopeAbundances[i] *
                        envelope.Peaks.OrderByDescending(p => p.intensity).First().intensity;


                    //if (envelope.Peaks[i + 1 + shift.Key].intensity < theoreticalIsotopeIntensity / 4.0 ||
                    //    envelope.Peaks[i + 1 + shift.Key].intensity > theoreticalIsotopeIntensity * 4.0)
                    //{
                    //    break;
                    //}

                    shift.Value.Add((envelope.Peaks[i + 1 + shift.Key].intensity, theoreticalIsotopeIntensity, isotopeMass));
                    if (shift.Key == 0)
                    {
                        experimentalIsotopeIntensities[i] = envelope.Peaks[i].intensity;
                    }
                }
            }

            double corr = Correlation.Pearson(massShiftToIsotopePeaks[0].Select(p => p.expIntensity), massShiftToIsotopePeaks[0].Select(p => p.theorIntensity));
            double corrShiftedLeft = Correlation.Pearson(massShiftToIsotopePeaks[-1].Select(p => p.expIntensity), massShiftToIsotopePeaks[-1].Select(p => p.theorIntensity));
            double corrShiftedRight = Correlation.Pearson(massShiftToIsotopePeaks[1].Select(p => p.expIntensity), massShiftToIsotopePeaks[1].Select(p => p.theorIntensity));

            //if (corrShiftedLeft - corr > 0.1) psm.AddIsotopeCorrelation(-1 * corrShiftedLeft);
            //else if (corrShiftedRight - corr > 0.1) psm.AddIsotopeCorrelation(-1 * corrShiftedRight);
            //else psm.AddIsotopeCorrelation(corr);
        }

        public IsotopicEnvelope FindEnvelope(PeptideSpectralMatch psm, MsDataScan scan)
        {
            // envelopeSize + 2 because we want to look on either side of the isotopic envelope
            int envelopeSize = _modifiedSequenceToIsotopicDistribution[psm.FullSequence].Count + 2;
            double precursorCharge =(double)psm.ScanPrecursorCharge;
            DoubleRange[] ranges = new DoubleRange[envelopeSize];
            double[] intensityValues = new double[envelopeSize];
            double[] mzValues = new double[envelopeSize];

            
            for (int i = 0; i < envelopeSize; i++)
            {
                DoubleRange interval = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(
                        (double)psm.PeptideMonisotopicMass + (i - 1) * Constants.C13MinusC12).First().AllowedInterval; 
                ranges[i] = interval;
                mzValues[i] = interval.Mean;
            }

            MzSpectrum spectrum = scan.MassSpectrum;

            for (int i = 0; i < spectrum.XArray.Length; i++)
            {
                if ((spectrum.XArray[i] - Constants.ProtonMass) * precursorCharge >= ranges[0].Minimum) 
                {
                    for (int j = 0; j < envelopeSize; j++)
                    {
                        while ((spectrum.XArray[i] - Constants.ProtonMass) * precursorCharge <= ranges[j].Maximum)
                        {
                            if (ranges[j].Contains((spectrum.XArray[i] - Constants.ProtonMass) * precursorCharge) )
                            {
                                mzValues[j] = (spectrum.XArray[i] - Constants.ProtonMass) * precursorCharge;
                                intensityValues[j] = spectrum.YArray[i];
                            }
                            i++;
                        }
                    }
                    break;
                }
            }

            IsotopicEnvelope envelope = new(
                MakeTupleList(mzValues, intensityValues),
                0, psm.ScanPrecursorCharge,
                0, 0, 0);

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
