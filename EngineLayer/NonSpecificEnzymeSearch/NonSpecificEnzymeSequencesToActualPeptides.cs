using System;
using Proteomics;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Chemistry;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    class NonSpecificEnzymeSequencesToActualPeptides : SequencesToActualProteinPeptidesEngine
    {
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        public NonSpecificEnzymeSequencesToActualPeptides(List<Psm>[] allPsms, List<Protein> proteinList, List<MassDiffAcceptor> massDiffAcceptors, Protease protease, int maxMissedCleavages, int? minPeptideLength, int? maxPeptideLength, InitiatorMethionineBehavior initiatorMethionineBehavior, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> variableModifications, int maxModificationIsoforms, List<string> nestedIds, bool addCompIons, List<ProductType> lp) : base(allPsms, proteinList, massDiffAcceptors, protease, maxMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications, variableModifications, maxModificationIsoforms, nestedIds, addCompIons, lp) { }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            bool classicAntigens = false;
            double precursorTolerance = 0;
            if (massDiffAcceptors.Count() > 1)
            {
                if (massDiffAcceptors[0].ToString().Contains("ppmAroundZero"))
                {
                    string name = massDiffAcceptors[0].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorTolerance = Convert.ToDouble(name.Substring(0, index));
                    classicAntigens = true;
                }
                else if (massDiffAcceptors[1].ToString().Contains("ppmAroundZero"))
                {
                    string name = massDiffAcceptors[1].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorTolerance = Convert.ToDouble(name.Substring(0, index));
                    classicAntigens = true;
                }
            }
            else { throw new MetaMorpheusException("Multiple searches required to run fast nonspecific enzyme searches"); }
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
            Dictionary<CompactPeptide, List<double>> compactPeptideToMassMatching = new Dictionary<CompactPeptide, List<double>>();

            //myAnalysisResults.AddText("Starting compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            //Status("Adding observed peptides to dictionary...", new List<string> { taskId });


            //Looking at the search results, generate a dictionary of keys for each unique CompactPeptide with empty values
            foreach (var psmListForASpecificSearchMode in allPsms) //should only be one
                if (psmListForASpecificSearchMode != null)
                    foreach (var psm in psmListForASpecificSearchMode)
                        if (psm != null)
                        {
                            foreach (var cp in psm.compactPeptides)
                            {
                                List<double> ld = new List<double>();
                                if (compactPeptideToMassMatching.TryGetValue(cp.Key, out ld))
                                {
                                    ld.Add(psm.ScanPrecursorMass);
                                }
                                else
                                {
                                    compactPeptideToProteinPeptideMatching.Add(cp.Key, new HashSet<PeptideWithSetModifications>()); //populate dictionary with all keys
                                    compactPeptideToMassMatching.Add(cp.Key, new List<double> { psm.ScanPrecursorMass });
                                }
                            }
                        }
            //myAnalysisResults.AddText("Ending compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            Dictionary<CompactPeptide, HashSet<CompactPeptideWithMass>> CPtoCPWM = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<CompactPeptideWithMass>()); //cannot populate values yet, as we want the theoretical masses of the peptides
            Dictionary<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>> CPWMtoPWSM = new Dictionary<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>>();
            int totalProteins = proteinList.Count;
            int proteinsSeen = 0;
            int old_progress = 0;
            var obj = new object();
            //Status("Adding possible sources to peptide dictionary...", new List<string> { taskId });
            //Populate the dictionary with possible sources for those ions
            //particularly tricky for single proteases, since each is more scan specific.
            if (protease.Name.Equals("singleN"))
            {
                Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
                {
                        //Digest protein into large peptide fragments and store in local1
                        Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> localCPtoPWSM = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                    for (int i = fff.Item1; i < fff.Item2; i++)
                    {
                        foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maxMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications, addCompIons))
                        {
                            foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModificationIsoforms, max_mods_for_peptide))
                            {
                                HashSet<PeptideWithSetModifications> v;
                                if (localCPtoPWSM.TryGetValue(new CompactPeptide(peptideWithSetModifications, addCompIons), out v))
                                    v.Add(peptideWithSetModifications);
                            }
                        }
                    }
                        //Foreach large peptide in localCPtoPWSM, find the precursor masses it's associated with and attempt to find other terminus. Store new compact peptide in local2
                        Dictionary<CompactPeptide, HashSet<CompactPeptideWithMass>> localCPtoCPWM = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<CompactPeptideWithMass>());
                    Dictionary<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>> localCPWMtoPWSM = new Dictionary<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>>();
                    foreach (KeyValuePair<CompactPeptide, HashSet<PeptideWithSetModifications>> kvp in localCPtoPWSM)
                    {
                        List<double> listScanPrecursorMasses;
                        if (compactPeptideToMassMatching.TryGetValue(kvp.Key, out listScanPrecursorMasses)) //do peaks match? Then lets modify double[] into compactpeptide
                            {
                            foreach (PeptideWithSetModifications pwsm in kvp.Value)
                            {
                                    //Determine if the precursor mass can be obtained within the acceptable margin of error.
                                    ModificationWithMass pep_n_term_variable_mod;
                                double initialMass = 0;
                                if (pwsm.allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
                                    foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                                        initialMass = pep_n_term_variable_mod.monoisotopicMass - nl;
                                else
                                    initialMass = 0;
                                double[] finalMass = new double[1];
                                foreach (double precursorMass in listScanPrecursorMasses)
                                {
                                    finalMass[0] = initialMass + waterMonoisotopicMass;
                                    int index = ComputePeptideIndexes(pwsm, finalMass, 1, 1, precursorMass, precursorTolerance);
                                    if (index >= 0 && (!minPeptideLength.HasValue | index >= minPeptideLength))
                                    {
                                            //generate correct sequence
                                            PeptideWithPossibleModifications tempPWPM = new PeptideWithPossibleModifications(pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedStartResidueInProtein + index - 1, pwsm.Protein, 0, pwsm.PeptideDescription, pwsm.modPep.thisDictionaryOfFixedMods, addCompIons);
                                        PeptideWithSetModifications tempPWSM = new PeptideWithSetModifications(tempPWPM, pwsm.allModsOneIsNterminus, pwsm.numFixedMods, addCompIons);
                                        HashSet<CompactPeptideWithMass> tempCPWMHashSet;
                                        CompactPeptideWithMass tempCPWM = new CompactPeptideWithMass(kvp.Key, finalMass[0]);
                                        if (localCPtoCPWM.TryGetValue(kvp.Key, out tempCPWMHashSet))
                                        {
                                            tempCPWMHashSet.Add(tempCPWM); //populate dictionary with all keys
                                            }
                                        HashSet<PeptideWithSetModifications> tempPWSMHashSet;
                                        if (localCPWMtoPWSM.TryGetValue(tempCPWM, out tempPWSMHashSet))
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
                            //we have our locals, now we need to merge into the master
                            //localCPtoPWSM is junk, because the PWSM is the large, uninformative peptide
                            //localCPWMtoPWSM has the good stuff in it

                            //Populate CPtoCPWM
                            foreach (KeyValuePair<CompactPeptide, HashSet<CompactPeptideWithMass>> kvp in localCPtoCPWM)
                        {
                            HashSet<CompactPeptideWithMass> tempCPWMHashSet;
                            if (CPtoCPWM.TryGetValue(kvp.Key, out tempCPWMHashSet))
                            {
                                foreach (CompactPeptideWithMass CPWM in kvp.Value)
                                {
                                    if (!tempCPWMHashSet.Contains(CPWM))
                                    {
                                        tempCPWMHashSet.Add(CPWM);
                                    }
                                }
                            }
                            else
                            {
                                tempCPWMHashSet = new HashSet<CompactPeptideWithMass>();
                                foreach (CompactPeptideWithMass CPWM in kvp.Value)
                                {
                                    if (!tempCPWMHashSet.Contains(CPWM))
                                    {
                                        tempCPWMHashSet.Add(CPWM);
                                    }
                                }
                                CPtoCPWM.Add(kvp.Key, tempCPWMHashSet);
                            }
                        }
                            //PopulateCPWMtoPWSM
                            foreach (KeyValuePair<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>> kvp in localCPWMtoPWSM)
                        {
                            HashSet<PeptideWithSetModifications> tempPWSMHashSet;
                            if (CPWMtoPWSM.TryGetValue(kvp.Key, out tempPWSMHashSet))
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
            else //if protease.Name.Equals("singleC")
            {
                Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
                {
                        //Digest protein into large peptide fragments and store in local1
                        Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> localCPtoPWSM = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                    for (int i = fff.Item1; i < fff.Item2; i++)
                    {
                        foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maxMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications, addCompIons))
                        {
                            foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModificationIsoforms, max_mods_for_peptide))
                            {
                                HashSet<PeptideWithSetModifications> v;
                                if (localCPtoPWSM.TryGetValue(new CompactPeptide(peptideWithSetModifications, addCompIons), out v))
                                    v.Add(peptideWithSetModifications);
                            }
                        }
                    }
                        //Foreach large peptide in localCPtoPWSM, find the precursor masses it's associated with and attempt to find other terminus. Store new compact peptide in local2
                        Dictionary<CompactPeptide, HashSet<CompactPeptideWithMass>> localCPtoCPWM = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<CompactPeptideWithMass>());
                    Dictionary<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>> localCPWMtoPWSM = new Dictionary<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>>();
                    foreach (KeyValuePair<CompactPeptide, HashSet<PeptideWithSetModifications>> kvp in localCPtoPWSM)
                    {
                        List<double> listScanPrecursorMasses;
                        if (compactPeptideToMassMatching.TryGetValue(kvp.Key, out listScanPrecursorMasses)) //do peaks match? Then lets modify double[] into compactpeptide
                            {
                            foreach (PeptideWithSetModifications pwsm in kvp.Value)
                            {
                                    //Determine if the precursor mass can be obtained within the acceptable margin of error.
                                    ModificationWithMass pep_n_term_variable_mod;
                                double initialMass = 0;
                                if (pwsm.allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
                                    foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                                        initialMass = pep_n_term_variable_mod.monoisotopicMass - nl;
                                else
                                    initialMass = 0;
                                double[] finalMass = new double[1];

                                foreach (double precursorMass in listScanPrecursorMasses)
                                {
                                    if (precursorMass > 2563.2148 && precursorMass < 2563.2150)
                                    { }
                                    finalMass[0] = initialMass + waterMonoisotopicMass;
                                    int index = ComputePeptideIndexes(pwsm, finalMass, pwsm.Length, -1, precursorMass, precursorTolerance);
                                    if (index >= 0 && (!minPeptideLength.HasValue | (pwsm.OneBasedEndResidueInProtein - (pwsm.OneBasedStartResidueInProtein + index - 2)) >= minPeptideLength))
                                    {
                                            //generate correct sequence
                                            PeptideWithPossibleModifications tempPWPM = new PeptideWithPossibleModifications(pwsm.OneBasedStartResidueInProtein + index - 1, pwsm.OneBasedEndResidueInProtein, pwsm.Protein, 0, pwsm.PeptideDescription, pwsm.modPep.thisDictionaryOfFixedMods, addCompIons);
                                        PeptideWithSetModifications tempPWSM = new PeptideWithSetModifications(tempPWPM, pwsm.allModsOneIsNterminus, pwsm.numFixedMods, addCompIons);
                                        HashSet<CompactPeptideWithMass> tempCPWMHashSet;
                                        CompactPeptideWithMass tempCPWM = new CompactPeptideWithMass(kvp.Key, finalMass[0]);
                                        if (localCPtoCPWM.TryGetValue(kvp.Key, out tempCPWMHashSet))
                                        {
                                            tempCPWMHashSet.Add(tempCPWM); //populate dictionary with all keys
                                            }
                                        HashSet<PeptideWithSetModifications> tempPWSMHashSet;
                                        if (localCPWMtoPWSM.TryGetValue(tempCPWM, out tempPWSMHashSet))
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
                            //we have our locals, now we need to merge into the master
                            //localCPtoPWSM is junk, because the PWSM is the large, uninformative peptide
                            //localCPWMtoPWSM has the good stuff in it

                            //Populate CPtoCPWM
                            foreach (KeyValuePair<CompactPeptide, HashSet<CompactPeptideWithMass>> kvp in localCPtoCPWM)
                        {
                            HashSet<CompactPeptideWithMass> tempCPWMHashSet;
                            if (CPtoCPWM.TryGetValue(kvp.Key, out tempCPWMHashSet))
                            {
                                foreach (CompactPeptideWithMass CPWM in kvp.Value)
                                {
                                    if (!tempCPWMHashSet.Contains(CPWM))
                                    {
                                        tempCPWMHashSet.Add(CPWM);
                                    }
                                    else { }
                                }
                            }
                            else
                            {
                                tempCPWMHashSet = new HashSet<CompactPeptideWithMass>();
                                foreach (CompactPeptideWithMass CPWM in kvp.Value)
                                {
                                    if (!tempCPWMHashSet.Contains(CPWM))
                                    {
                                        tempCPWMHashSet.Add(CPWM);
                                    }
                                    else { }
                                }
                                CPtoCPWM.Add(kvp.Key, tempCPWMHashSet);
                            }
                        }
                            //PopulateCPWMtoPWSM
                            int i = 0;
                        foreach (KeyValuePair<CompactPeptideWithMass, HashSet<PeptideWithSetModifications>> kvp in localCPWMtoPWSM)
                        {
                            i++;
                            HashSet<PeptideWithSetModifications> tempPWSMHashSet;
                            if (CPWMtoPWSM.TryGetValue(kvp.Key, out tempPWSMHashSet))
                            {
                                foreach (PeptideWithSetModifications PWSM in kvp.Value)
                                {
                                    if (!tempPWSMHashSet.Contains(PWSM))
                                    {
                                        tempPWSMHashSet.Add(PWSM);
                                    }
                                    else { }
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
                                    else { }
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
            //Add these corrected CPs to a new Dictionary for output
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> MHCOutput = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
            foreach (var psmListForAspecificSerchMode in allPsms) //should only be one
                if (psmListForAspecificSerchMode != null)
                    foreach (var psm in psmListForAspecificSerchMode)
                        if (psm != null)
                        {
                            List<CompactPeptide> keysToUpdate = new List<CompactPeptide>();
                            foreach (var cp in psm.compactPeptides)
                            {
                                keysToUpdate.Add(cp.Key);
                            }

                            foreach (CompactPeptide keyToUpdate in keysToUpdate)
                            {
                                Tuple<int, HashSet<PeptideWithSetModifications>> tempTuple;
                                psm.compactPeptides.TryGetValue(keyToUpdate, out tempTuple);
                                psm.compactPeptides.Remove(keyToUpdate);
                                //get new CompactPeptide
                                HashSet<CompactPeptideWithMass> tempCWPMHashSet;
                                if (CPtoCPWM.TryGetValue(keyToUpdate, out tempCWPMHashSet))
                                {
                                    foreach (CompactPeptideWithMass cpwm in tempCWPMHashSet)
                                    {
                                        if (Math.Abs((psm.ScanPrecursorMass - cpwm.precursorMass) / (cpwm.precursorMass) * 1e6) < precursorTolerance)
                                        {
                                            CompactPeptide newKey = new CompactPeptide(keyToUpdate.CTerminalMasses, keyToUpdate.NTerminalMasses, cpwm.precursorMass, keyToUpdate.addCompIons);
                                            if (!psm.compactPeptides.ContainsKey(newKey))
                                            {
                                                psm.compactPeptides.Add(newKey, tempTuple);
                                            }
                                            HashSet<PeptideWithSetModifications> localTempPWSMHashSet;
                                            if (CPWMtoPWSM.TryGetValue(cpwm, out localTempPWSMHashSet))
                                            {
                                                HashSet<PeptideWithSetModifications> masterTempPWSMHashSet;
                                                if (MHCOutput.TryGetValue(newKey, out masterTempPWSMHashSet))
                                                {
                                                    foreach (PeptideWithSetModifications pwsm in localTempPWSMHashSet)
                                                        masterTempPWSMHashSet.Add(pwsm);
                                                }
                                                else
                                                {
                                                    MHCOutput.Add(newKey, localTempPWSMHashSet);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if (psm.compactPeptides.Count() == 0)
                            {
                                throw new MissingMemberException("There was an incorrect mass calculation resulting in the loss of scan " + psm.ScanNumber);
                            }
                        }
            return new SequencesToActualProteinPeptidesEngineResults(this, MHCOutput);
        }

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
    }
}
