using Chemistry;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    public class NonSpecificEnzymeSequencesToActualPeptides : SequencesToActualProteinPeptidesEngine
    {
        #region Private Fields

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private readonly List<MassDiffAcceptor> massDiffAcceptors;

        #endregion Private Fields

        #region Public Constructors

        public NonSpecificEnzymeSequencesToActualPeptides(List<Psm>[] allPsms, List<Protein> proteinList, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> variableModifications, TerminusType terminusType, CommonParameters CommonParameters, List<MassDiffAcceptor> massDiffAcceptors, List<string> nestedIds) : base(allPsms, proteinList, fixedModifications, variableModifications, terminusType, CommonParameters, nestedIds)
        {
            this.massDiffAcceptors = massDiffAcceptors;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double precursorTolerance = 0;
            if (massDiffAcceptors.Count() > 1)
            {
                if (massDiffAcceptors[0].ToString().Contains("ppmAroundZero"))
                {
                    string name = massDiffAcceptors[0].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorTolerance = Convert.ToDouble(name.Substring(0, index));
                }
                else if (massDiffAcceptors[1].ToString().Contains("ppmAroundZero"))
                {
                    string name = massDiffAcceptors[1].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorTolerance = Convert.ToDouble(name.Substring(0, index));
                }
            }
            else { throw new MetaMorpheusException("Multiple searches required to run NonSpecific Search"); }
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            Dictionary<CompactPeptideBase, List<double>> compactPeptideToMassMatching = new Dictionary<CompactPeptideBase, List<double>>();

            //Looking at the search results, generate a dictionary of keys for each unique CompactPeptide with empty values
            foreach (var psmListForASpecificSearchMode in allPsms) //should only be one
                if (psmListForASpecificSearchMode != null)
                    foreach (var psm in psmListForASpecificSearchMode)
                        if (psm != null)
                        {
                            foreach (var cp in psm.CompactPeptides)
                            {
                                if (compactPeptideToMassMatching.TryGetValue(cp.Key, out List<double> ld))
                                {
                                    ld.Add(psm.ScanPrecursorMass);
                                }
                                else
                                {
                                    compactPeptideToProteinPeptideMatching.Add(cp.Key as CompactPeptideBase, new HashSet<PeptideWithSetModifications>()); //populate dictionary with all keys
                                    compactPeptideToMassMatching.Add(cp.Key, new List<double> { psm.ScanPrecursorMass });
                                }
                            }
                        }

            //CP==CompactPeptide
            //CPWM==CompactPeptideWithMass (Patched to respresent a double)
            //PWSM==PeptideWithSetModification
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> CPWMtoPWSM = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            int totalProteins = proteinList.Count;
            int proteinsSeen = 0;
            int old_progress = 0;
            var obj = new object();
            //Status("Adding possible sources to peptide dictionary...", new List<string> { taskId });
            //Populate the dictionary with possible sources for those ions
            //particularly tricky for single proteases, since each is more scan specific.
            if (terminusType == TerminusType.N)
            {
                Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
                {
                    //Digest protein into large peptide fragments and store in local1
                    Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> localCPtoPWSM = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key as CompactPeptideBase, b => new HashSet<PeptideWithSetModifications>());
                    for (int i = fff.Item1; i < fff.Item2; i++)
                    {
                        foreach (var peptideWithPossibleModifications in proteinList[i].Digest(CommonParameters.Protease, (int)CommonParameters.MaxMissedCleavages, CommonParameters.MinPeptideLength, CommonParameters.MaxPeptideLength, CommonParameters.InitiatorMethionineBehavior, fixedModifications))
                        {
                            foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, (int)CommonParameters.MaxModificationIsoforms, (int)CommonParameters.Max_mods_for_peptide))
                            {
                                if (localCPtoPWSM.TryGetValue(new CompactPeptide(peptideWithSetModifications, terminusType), out HashSet<PeptideWithSetModifications> v))
                                    v.Add(peptideWithSetModifications);
                            }
                        }
                    }
                    //Foreach large peptide in localCPtoPWSM, find the precursor masses it's associated with and attempt to find other terminus. Store new compact peptide in local2
                    //CP==CompactPeptide
                    //CPWM==CompactPeptideWithMass (Patched to respresent a double)
                    //PWSM==PeptideWithSetModificationDictionary
                    Dictionary<CompactPeptideWithModifiedMass, HashSet<PeptideWithSetModifications>> localCPWMtoPWSM = new Dictionary<CompactPeptideWithModifiedMass, HashSet<PeptideWithSetModifications>>();
                    foreach (KeyValuePair<CompactPeptideBase, HashSet<PeptideWithSetModifications>> kvp in localCPtoPWSM) //foreach theoretical kvp
                    {
                        if (compactPeptideToMassMatching.TryGetValue(kvp.Key, out List<double> listScanPrecursorMasses)) //get list of theoretical precursor masses that have been found and are associated with compactPeptide
                        {
                            foreach (PeptideWithSetModifications pwsm in kvp.Value)
                            {
                                //Determine if the precursor mass can be obtained within the acceptable margin of error.
                                double initialMass = 0;
                                if (pwsm.allModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                                    foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                                        initialMass = pep_n_term_variable_mod.monoisotopicMass - nl;
                                else
                                    initialMass = 0;
                                double[] finalMass = new double[1];
                                foreach (double precursorMass in listScanPrecursorMasses) //foreach precursor
                                {
                                    finalMass[0] = initialMass + waterMonoisotopicMass; //This is the starting mass of the final mass
                                    int index = ComputePeptideIndexes(pwsm, finalMass, 1, 1, precursorMass, precursorTolerance);
                                    if (index >= 0 && (!CommonParameters.MinPeptideLength.HasValue | index >= CommonParameters.MinPeptideLength))
                                    {
                                        //generate correct sequence
                                        PeptideWithSetModifications tempPWSM = new PeptideWithSetModifications(pwsm, pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedStartResidueInProtein + index - 1);
                                        double modifiedMass = finalMass[0];
                                        CompactPeptideWithModifiedMass tempCPWM = new CompactPeptideWithModifiedMass(kvp.Key, modifiedMass);
                                        tempCPWM.SwapMonoisotopicMassWithModifiedMass();
                                        if (localCPWMtoPWSM.TryGetValue(tempCPWM, out HashSet<PeptideWithSetModifications> tempPWSMHashSet))
                                        {
                                            tempPWSMHashSet.Add(tempPWSM);
                                        }
                                        else
                                        {
                                            localCPWMtoPWSM.Add(tempCPWM, new HashSet<PeptideWithSetModifications> { tempPWSM });
                                        }
                                    }
                                }
                            }
                        }
                    }
                    lock (obj)
                    {
                        //PopulateCPWMtoPWSM
                        foreach (KeyValuePair<CompactPeptideWithModifiedMass, HashSet<PeptideWithSetModifications>> kvp in localCPWMtoPWSM)
                        {
                            if (CPWMtoPWSM.TryGetValue(kvp.Key, out HashSet<PeptideWithSetModifications> tempPWSMHashSet))
                            {
                                foreach (PeptideWithSetModifications PWSM in kvp.Value)
                                {
                                    if (!tempPWSMHashSet.Contains(PWSM))
                                    {
                                        tempPWSMHashSet.Add(PWSM);
                                    }
                                }
                            }
                            else
                            {
                                tempPWSMHashSet = new HashSet<PeptideWithSetModifications>();
                                foreach (PeptideWithSetModifications PWSM in kvp.Value)
                                {
                                    if (!tempPWSMHashSet.Contains(PWSM))
                                    {
                                        tempPWSMHashSet.Add(PWSM);
                                    }
                                }
                                CPWMtoPWSM.Add(kvp.Key, tempPWSMHashSet);
                            }
                        }

                        //Everything has been populated!
                        //now need to
                        //record progress
                        proteinsSeen += fff.Item2 - fff.Item1;
                        var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                        if (new_progress > old_progress)
                        {
                            //ReportProgress(new ProgressEventArgs(new_progress, "In adding possible sources to peptide dictionary loop", nestedIds));
                            old_progress = new_progress;
                        }
                    }
                });
            }
            else //if (terminusType==TerminusType.C)
            {
                Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
                {
                    //Digest protein into large peptide fragments and store in local1
                    Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> localCPtoPWSM = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                    for (int i = fff.Item1; i < fff.Item2; i++)
                    {
                        foreach (var peptideWithPossibleModifications in proteinList[i].Digest(CommonParameters.Protease, (int)CommonParameters.MaxMissedCleavages, CommonParameters.MinPeptideLength, CommonParameters.MaxPeptideLength, CommonParameters.InitiatorMethionineBehavior, fixedModifications))
                        {
                            foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, (int)CommonParameters.MaxModificationIsoforms, (int)CommonParameters.Max_mods_for_peptide))
                            {
                                if (localCPtoPWSM.TryGetValue(new CompactPeptide(peptideWithSetModifications, terminusType), out HashSet<PeptideWithSetModifications> v))
                                    v.Add(peptideWithSetModifications);
                            }
                        }
                    }
                    //Foreach large peptide in localCPtoPWSM, find the precursor masses it's associated with and attempt to find other terminus. Store new compact peptide in local2
                    //CP==CompactPeptide
                    //CPWM==CompactPeptideWithMass (Patched to respresent a double)
                    //PWSM==PeptideWithSetModificationDictionary
                    Dictionary<CompactPeptideWithModifiedMass, HashSet<PeptideWithSetModifications>> localCPWMtoPWSM = new Dictionary<CompactPeptideWithModifiedMass, HashSet<PeptideWithSetModifications>>();
                    foreach (KeyValuePair<CompactPeptideBase, HashSet<PeptideWithSetModifications>> kvp in localCPtoPWSM) //foreach theoretical kvp
                    {
                        if (compactPeptideToMassMatching.TryGetValue(kvp.Key, out List<double> listScanPrecursorMasses)) //do peaks match? Then lets modify double[] into compactpeptide
                        {
                            foreach (PeptideWithSetModifications pwsm in kvp.Value)
                            {
                                //Determine if the precursor mass can be obtained within the acceptable margin of error.
                                double initialMass = 0;
                                if (pwsm.allModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                                    foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                                        initialMass = pep_n_term_variable_mod.monoisotopicMass - nl;
                                else
                                    initialMass = 0;
                                double[] finalMass = new double[1];

                                foreach (double precursorMass in listScanPrecursorMasses)
                                {
                                    finalMass[0] = initialMass + waterMonoisotopicMass;
                                    int index = ComputePeptideIndexes(pwsm, finalMass, pwsm.Length, -1, precursorMass, precursorTolerance);
                                    if (index >= 0 && (!CommonParameters.MinPeptideLength.HasValue | (pwsm.OneBasedEndResidueInProtein - (pwsm.OneBasedStartResidueInProtein + index - 2)) >= CommonParameters.MinPeptideLength))
                                    {
                                        //generate correct sequence
                                        PeptideWithSetModifications tempPWSM = new PeptideWithSetModifications(pwsm, pwsm.OneBasedStartResidueInProtein + index - 1, pwsm.OneBasedEndResidueInProtein);
                                        double modifiedMass = finalMass[0];
                                        CompactPeptideWithModifiedMass tempCPWM = new CompactPeptideWithModifiedMass(kvp.Key, modifiedMass);
                                        tempCPWM.SwapMonoisotopicMassWithModifiedMass();
                                        if (localCPWMtoPWSM.TryGetValue(tempCPWM, out HashSet<PeptideWithSetModifications> tempPWSMHashSet))
                                        {
                                            tempPWSMHashSet.Add(tempPWSM);
                                        }
                                        else
                                        {
                                            localCPWMtoPWSM.Add(tempCPWM, new HashSet<PeptideWithSetModifications> { tempPWSM });
                                        }
                                    }
                                }
                            }
                        }
                    }
                    lock (obj)
                    {
                        //PopulateCPWMtoPWSM
                        int i = 0;
                        foreach (KeyValuePair<CompactPeptideWithModifiedMass, HashSet<PeptideWithSetModifications>> kvp in localCPWMtoPWSM)
                        {
                            i++;
                            if (CPWMtoPWSM.TryGetValue(kvp.Key, out HashSet<PeptideWithSetModifications> tempPWSMHashSet))
                            {
                                foreach (PeptideWithSetModifications PWSM in kvp.Value)
                                {
                                    if (!tempPWSMHashSet.Contains(PWSM))
                                    {
                                        tempPWSMHashSet.Add(PWSM);
                                    }
                                }
                            }
                            else
                            {
                                tempPWSMHashSet = new HashSet<PeptideWithSetModifications>();
                                foreach (PeptideWithSetModifications PWSM in kvp.Value)
                                {
                                    if (!tempPWSMHashSet.Contains(PWSM))
                                    {
                                        tempPWSMHashSet.Add(PWSM);
                                    }
                                }
                                CPWMtoPWSM.Add(kvp.Key, tempPWSMHashSet);
                            }
                        }

                        //Everything has been populated!
                        //now need to
                        //record progress
                        proteinsSeen += fff.Item2 - fff.Item1;
                        var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                        if (new_progress > old_progress)
                        {
                            //ReportProgress(new ProgressEventArgs(new_progress, "In adding possible sources to peptide dictionary loop", nestedIds));
                            old_progress = new_progress;
                        }
                    }
                });
            }
            //with filled CPtoCPWM and CPWMtoPWSM, convert psm objects to corrected CP mass
            foreach (var psmListForASpecificSearchMode in allPsms) //should only be one
                if (psmListForASpecificSearchMode != null)
                    foreach (var psm in psmListForASpecificSearchMode)
                        if (psm != null)
                        {
                            foreach (KeyValuePair<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>> kvp in psm.CompactPeptides)
                            {
                                (kvp.Key as CompactPeptideWithModifiedMass).SwapMonoisotopicMassWithModifiedMass();
                                //Change CPWM to reflect actual CP
                                if (CPWMtoPWSM.TryGetValue(kvp.Key, out HashSet<PeptideWithSetModifications> misplacedPWSMs))
                                {
                                    (kvp.Key as CompactPeptideWithModifiedMass).CropTerminalMasses(terminusType);
                                    if (CPWMtoPWSM.TryGetValue(kvp.Key, out HashSet<PeptideWithSetModifications> wellPlacedPWSMs))
                                    {
                                        foreach (PeptideWithSetModifications PWSM in misplacedPWSMs)
                                            wellPlacedPWSMs.Add(PWSM);
                                    }
                                    else
                                    {
                                        CPWMtoPWSM.Add(kvp.Key, misplacedPWSMs);
                                    }
                                }
                            }
                            psm.CompactCompactPeptides();
                        }
            return new SequencesToActualProteinPeptidesEngineResults(this, CPWMtoPWSM);
        }

        #endregion Protected Methods

        #region Private Methods

        private int ComputePeptideIndexes(PeptideWithSetModifications yyy, double[] prevMass, int oneBasedIndexToLookAt, int direction, double precursorMass, double ppmTolerance)
        {
            ModificationWithMass residue_variable_mod = null;
            do
            {
                prevMass[0] += Residue.ResidueMonoisotopicMass[yyy[oneBasedIndexToLookAt - 1]];

                yyy.allModsOneIsNterminus.TryGetValue(oneBasedIndexToLookAt + 1, out residue_variable_mod);
                if (residue_variable_mod == null)
                {
                    if (Math.Abs((precursorMass - prevMass[0]) / (prevMass[0]) * 1e6) < ppmTolerance)
                    {
                        return oneBasedIndexToLookAt;
                    }
                }
                else if (residue_variable_mod.neutralLosses.Count == 1)
                {
                    prevMass[0] += residue_variable_mod.monoisotopicMass - residue_variable_mod.neutralLosses.First();
                    if (Math.Abs((precursorMass - prevMass[0]) / (prevMass[0]) * 1e6) < ppmTolerance)
                    {
                        return oneBasedIndexToLookAt;
                    }
                }
                else
                {
                    foreach (double nl in residue_variable_mod.neutralLosses)
                    {
                        prevMass[0] = prevMass[0] + residue_variable_mod.monoisotopicMass - nl;
                        if (Math.Abs((precursorMass - prevMass[0]) / (prevMass[0]) * 1e6) < ppmTolerance)
                        {
                            return oneBasedIndexToLookAt;
                        }
                        if ((direction == 1 && oneBasedIndexToLookAt + direction < yyy.Length) ||
                            (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                            return ComputePeptideIndexes(yyy, prevMass, oneBasedIndexToLookAt + direction, direction, precursorMass, ppmTolerance);
                    }
                    break;
                }
                oneBasedIndexToLookAt += direction;
            } while ((oneBasedIndexToLookAt >= 1 && direction == -1) || (oneBasedIndexToLookAt <= yyy.Length && direction == 1));
            return -1;
        }

        #endregion Private Methods
    }
}