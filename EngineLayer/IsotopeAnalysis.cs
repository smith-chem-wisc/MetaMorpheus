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

        public void ScoreEnvelopes()
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
                    
                }
                
            }
        }

        public void FindEnvelope(PeptideSpectralMatch psm, MsDataScan scan)
        {
            DoubleRange monoMass = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass((double)psm.PeptideMonisotopicMass).
                Select(i => i.AllowedInterval).First();
            DoubleRange mostAbundantMass = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(psm.PeakFindingMass).
                Select(i => i.AllowedInterval).First();
            MzSpectrum spectrum = scan.MassSpectrum;
            for (int i = 0; i < 200; i++)
            {
                if (monoMass.Contains(spectrum.XArray[i])) 
                {

                }
            }
            bool y = x.Contains(2.0);
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
