using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Chemistry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Analysis
{
    public class AnalysisEngine : MyEngine
    {

        #region Private Fields

        private const int max_mods_for_peptide = 3;
        private readonly double binTol;
        private readonly int maximumMissedCleavages;
        private readonly int maxModIsoforms;
        private readonly PsmParent[][] newPsms;
        private readonly List<Protein> proteinList;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly List<ModificationWithMass> fixedModifications;
        private readonly List<ModificationWithMass> localizeableModifications;
        private readonly Protease protease;
        private readonly List<SearchMode> searchModes;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private readonly Tolerance fragmentTolerance;
        private readonly Action<BinTreeStructure, string> writeHistogramPeaksAction;
        private readonly Action<List<NewPsmWithFdr>, string> writePsmsAction;
        private readonly Action<List<ProteinGroup>, string> writeProteinGroupsAction;
        private readonly bool doParsimony;
        private readonly bool noOneHitWonders;
        private readonly bool doHistogramAnalysis;
        private readonly bool quantify;
        private readonly double quantifyRtTol;
        private readonly double quantifyPpmTol;
        private readonly List<ProductType> lp;
        private readonly InitiatorMethionineBehavior initiatorMethionineBehavior;
        private readonly List<string> nestedIds;
        private Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;

        #endregion Private Fields

        #region Public Constructors

        public AnalysisEngine(PsmParent[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFdr>, string> action2, Action<List<ProteinGroup>, string> action3, bool doParsimony, bool noOneHitWonders, int maximumMissedCleavages, int maxModIsoforms, bool doHistogramAnalysis, List<ProductType> lp, double binTol, InitiatorMethionineBehavior initiatorMethionineBehavior, List<string> nestedIds, bool Quantify, double QuantifyRtTol, double QuantifyPpmTol)
        {
            this.doParsimony = doParsimony;
            this.noOneHitWonders = noOneHitWonders;
            this.doHistogramAnalysis = doHistogramAnalysis;
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.protease = protease;
            this.searchModes = searchModes;
            this.myMsDataFile = myMSDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.writeHistogramPeaksAction = action1;
            this.writePsmsAction = action2;
            this.writeProteinGroupsAction = action3;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.maxModIsoforms = maxModIsoforms;
            this.lp = lp;
            this.binTol = binTol;
            this.initiatorMethionineBehavior = initiatorMethionineBehavior;
            this.quantify = Quantify;
            this.quantifyRtTol = QuantifyRtTol;
            this.quantifyPpmTol = QuantifyPpmTol;
            this.nestedIds = nestedIds;
        }

        #endregion Public Constructors

        #region Public Methods

        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> ApplyProteinParsimony(out List<ProteinGroup> proteinGroups)
        {
            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var parsimonyDict = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var proteinsWithUniquePeptides = new Dictionary<Protein, HashSet<CompactPeptide>>();

            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                // if a peptide is associated with a decoy protein, remove all target protein associations with the peptide
                var decoyPeptidesHere = new HashSet<PeptideWithSetModifications>(kvp.Value.Where(p => p.Protein.IsDecoy));
                if (decoyPeptidesHere.Any())
                {
                    var pepWithSetModsToRemove = new HashSet<PeptideWithSetModifications>(kvp.Value.Except(decoyPeptidesHere));

                    foreach (var pepWithSetMods in pepWithSetModsToRemove)
                        kvp.Value.Remove(pepWithSetMods);
                }

                // if a peptide is associated with a contaminant protein, remove all target protein associations with the peptide
                var contaminantPeptidesHere = new HashSet<PeptideWithSetModifications>(kvp.Value.Where(p => p.Protein.IsContaminant));
                if (contaminantPeptidesHere.Any())
                {
                    var pepWithSetModsToRemove = new HashSet<PeptideWithSetModifications>(kvp.Value.Except(contaminantPeptidesHere));

                    foreach (var pepWithSetMods in pepWithSetModsToRemove)
                        kvp.Value.Remove(pepWithSetMods);
                }

                // finds unique peptides (peptides that can belong to only one protein)
                var proteinListHere = new HashSet<Protein>(kvp.Value.Select(p => p.Protein));

                if (proteinListHere.Count == 1)
                {
                    if (!proteinsWithUniquePeptides.ContainsKey(proteinListHere.First()))
                    {
                        var peptides = new HashSet<CompactPeptide>();
                        peptides.Add(kvp.Key);
                        proteinsWithUniquePeptides.Add(proteinListHere.First(), peptides);
                    }
                    else
                    {
                        HashSet<CompactPeptide> peptides;
                        proteinsWithUniquePeptides.TryGetValue(proteinListHere.First(), out peptides);
                        peptides.Add(kvp.Key);
                    }
                }
            }

            // makes dictionary with proteins as keys and list of associated peptides as the value (makes parsimony algo easier)
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinToPeptidesMatching.ContainsKey(peptide.Protein))
                    {
                        var peptides = new HashSet<CompactPeptide>();
                        peptides.Add(kvp.Key);
                        proteinToPeptidesMatching.Add(peptide.Protein, peptides);
                    }
                    else
                    {
                        HashSet<CompactPeptide> peptides;
                        proteinToPeptidesMatching.TryGetValue(peptide.Protein, out peptides);
                        peptides.Add(kvp.Key);
                    }
                }
            }

            // build protein list for each peptide before parsimony has been applied
            var peptideBaseSeqProteinListMatch = new Dictionary<string, HashSet<Protein>>();
            foreach (var kvp in proteinToPeptidesMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    var proteinListHere = new HashSet<Protein>();
                    string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (!peptideBaseSeqProteinListMatch.ContainsKey(peptideBaseSequence))
                    {
                        proteinListHere.Add(kvp.Key);
                        peptideBaseSeqProteinListMatch.Add(peptideBaseSequence, proteinListHere);
                    }
                    else
                    {
                        peptideBaseSeqProteinListMatch.TryGetValue(peptideBaseSequence, out proteinListHere);
                        proteinListHere.Add(kvp.Key);
                    }
                }
            }

            // dictionary associates proteins w/ unused base seqs
            var algDictionary = new Dictionary<Protein, HashSet<string>>();
            foreach (var kvp in proteinToPeptidesMatching)
            {
                var newPeptideBaseSeqs = new HashSet<string>(kvp.Value.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)));
                algDictionary.Add(kvp.Key, newPeptideBaseSeqs);
            }

            int numNewSeqs;
            while (true)
            {
                numNewSeqs = algDictionary.Max(p => p.Value.Count);

                if (numNewSeqs == 0)
                    break;

                // gets list of proteins with the most unaccounted-for peptide base sequences
                var possibleBestProteinList = algDictionary.Where(p => p.Value.Count == numNewSeqs).ToList();
                Protein bestProtein = possibleBestProteinList.First().Key;

                HashSet<string> newSeqs;
                algDictionary.TryGetValue(bestProtein, out newSeqs);
                newSeqs = new HashSet<string>(newSeqs);

                // may need to select different protein
                if (possibleBestProteinList.Count > 1)
                {
                    var proteinsWithTheseBaseSeqs = new HashSet<Protein>();

                    foreach (var kvp in possibleBestProteinList)
                    {
                        if (newSeqs.IsSubsetOf(kvp.Value))
                            proteinsWithTheseBaseSeqs.Add(kvp.Key);
                    }

                    if (proteinsWithTheseBaseSeqs.Count > 1)
                    {
                        var dict = new Dictionary<Protein, HashSet<string>>();
                        foreach (var protein in proteinsWithTheseBaseSeqs)
                        {
                            HashSet<CompactPeptide> t;
                            proteinToPeptidesMatching.TryGetValue(protein, out t);
                            var seqs = new HashSet<string>(t.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)));
                            dict.Add(protein, seqs);
                        }

                        bestProtein = dict.OrderByDescending(kvp => kvp.Value.Count).First().Key;
                    }
                }
                
                HashSet<CompactPeptide> l;
                proteinToPeptidesMatching.TryGetValue(bestProtein, out l);
                parsimonyDict.Add(bestProtein, l);
                
                // remove used peptides from their proteins
                foreach (var newBaseSeq in newSeqs)
                {
                    HashSet<Protein> proteinsWithThisPeptide;
                    peptideBaseSeqProteinListMatch.TryGetValue(newBaseSeq, out proteinsWithThisPeptide);

                    foreach (var protein in proteinsWithThisPeptide)
                    {
                        HashSet<string> thisProteinsBaseSeqs;
                        algDictionary.TryGetValue(protein, out thisProteinsBaseSeqs);
                        thisProteinsBaseSeqs.Remove(newBaseSeq);
                    }
                }
            }

            // build protein group after parsimony (each group only has 1 protein at this point)
            proteinGroups = new List<ProteinGroup>();
            foreach (var kvp in parsimonyDict)
            {
                HashSet<CompactPeptide> uniquePeptidesHere;
                proteinsWithUniquePeptides.TryGetValue(kvp.Key, out uniquePeptidesHere);
                if (uniquePeptidesHere == null)
                    uniquePeptidesHere = new HashSet<CompactPeptide>();
                proteinGroups.Add(new ProteinGroup(new HashSet<Protein>() { kvp.Key }, kvp.Value, uniquePeptidesHere));
            }

            // grab indistinguishable proteins ("if" conditions are to narrow search space)
            foreach (var proteinGroup in proteinGroups)
            {
                if (!proteinGroup.TotalUniquePeptideList.Any())
                {
                    foreach (var kvp in proteinToPeptidesMatching)
                    {
                        if (!proteinsWithUniquePeptides.ContainsKey(kvp.Key))
                        {
                            // prevents looking at itself
                            if (!parsimonyDict.ContainsKey(kvp.Key))
                            {
                                if (kvp.Value.Count == proteinGroup.TotalPeptideList.Count)
                                {
                                    if (kvp.Value.SetEquals(proteinGroup.TotalPeptideList))
                                    {
                                        proteinGroup.Proteins.Add(kvp.Key);
                                        parsimonyDict.Add(kvp.Key, kvp.Value);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            // build protein list for each peptide after parsimony has been applied (helps make return dictionary)
            var peptideProteinListMatch = new Dictionary<CompactPeptide, HashSet<Protein>>();
            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    var proteinListHere = new HashSet<Protein>();

                    if (!peptideProteinListMatch.ContainsKey(peptide))
                    {
                        proteinListHere.Add(kvp.Key);
                        peptideProteinListMatch.Add(peptide, proteinListHere);
                    }
                    else
                        peptideProteinListMatch[peptide].Add(kvp.Key);
                }
            }

            // constructs return dictionary (only use parsimony proteins for the new PeptideWithSetModifications list)
            var answer = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!answer.ContainsKey(peptide))
                    {
                        // find CompactPeptide's original (unparsimonious) peptide matches
                        var oldPepsWithSetMods = compactPeptideToProteinPeptideMatching[peptide];

                        // get the CompactPeptide's protein list after parsimony
                        var proteinListHere = peptideProteinListMatch[peptide];

                        // get the peptides that belong to the post-parsimony protein(s) only
                        var newPeptides = new HashSet<PeptideWithSetModifications>(oldPepsWithSetMods.Where(p => proteinListHere.Contains(p.Protein)));

                        answer.Add(peptide, newPeptides);
                    }
                }
            }
            
            Status("Finished Parsimony", nestedIds);


            compactPeptideToProteinPeptideMatching = answer;

            // returns for test class (TODO**: remove)
            return answer;
        }

        public void ScoreProteinGroups(List<ProteinGroup> proteinGroups, List<NewPsmWithFdr> psmList)
        {
            Status("Scoring protein groups...", nestedIds);

            Dictionary<CompactPeptide, HashSet<ProteinGroup>> peptideToProteinGroupMatching = new Dictionary<CompactPeptide, HashSet<ProteinGroup>>();
            HashSet<CompactPeptide> allRazorPeptides = new HashSet<CompactPeptide>();
            HashSet<ProteinGroup> proteinGroupsToRemove = new HashSet<ProteinGroup>();

            // add each protein group's psm's
            Dictionary<CompactPeptide, NewPsmWithFdr> psmToCompactPeptideMatching = new Dictionary<CompactPeptide, NewPsmWithFdr>();
            foreach (var psm in psmList)
            {
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications);
                if (!psmToCompactPeptideMatching.ContainsKey(peptide))
                    psmToCompactPeptideMatching.Add(peptide, psm);
            }

            foreach (var proteinGroup in proteinGroups)
            {
                foreach (var peptide in proteinGroup.TotalPeptideList)
                {
                    // build PSM list for scoring
                    NewPsmWithFdr psm;
                    psmToCompactPeptideMatching.TryGetValue(peptide, out psm);
                    if (psm != null)
                        proteinGroup.TotalPsmList.Add(psm);

                    // build PeptideWithSetMod list to calc sequence coverage
                    HashSet<PeptideWithSetModifications> peptidesWithSetMods;
                    compactPeptideToProteinPeptideMatching.TryGetValue(peptide, out peptidesWithSetMods);
                    foreach (var pep in peptidesWithSetMods)
                    {
                        proteinGroup.TotalPeptideWithSetModsList.Add(pep);
                    }
                }
            }

            // score the group and get number of protein groups per peptide
            foreach (var proteinGroup in proteinGroups)
            {
                // score the group (scoring algorithm defined in the ProteinGroup class)
                proteinGroup.ScoreThisProteinGroup(variableModifications, localizeableModifications, fixedModifications);

                // for finding razor peptides later
                foreach (var peptide in proteinGroup.StrictPeptideList)
                {
                    HashSet<ProteinGroup> proteinGroupsHere = new HashSet<ProteinGroup>();
                    if (peptideToProteinGroupMatching.ContainsKey(peptide))
                    {
                        peptideToProteinGroupMatching.TryGetValue(peptide, out proteinGroupsHere);
                        proteinGroupsHere.Add(proteinGroup);
                    }
                    else
                    {
                        proteinGroupsHere.Add(proteinGroup);
                        peptideToProteinGroupMatching.Add(peptide, proteinGroupsHere);
                    }
                }
            }

            // merge protein groups that are indistinguishable after scoring
            var pg = proteinGroups.OrderByDescending(p => p.proteinGroupScore).ToList();
            for (int i = 0; i < (pg.Count - 1); i++)
            {
                if (pg[i].proteinGroupScore == pg[i + 1].proteinGroupScore && pg[i].proteinGroupScore != 0)
                {
                    // get all protein groups with the exact same score
                    var pgsWithThisScore = pg.Where(p => p.proteinGroupScore == pg[i].proteinGroupScore).ToList();

                    // check to make sure they have the same peptides, then merge them
                    foreach (var p in pgsWithThisScore)
                    {
                        if (p != pg[i] && p.StrictPeptideList.SetEquals(pg[i].StrictPeptideList))
                        {
                            pg[i].MergeProteinGroupWith(p);
                        }
                    }
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                if (proteinGroup.proteinGroupScore == 0)
                    proteinGroupsToRemove.Add(proteinGroup);
            }

            // remove empty protein groups (peptides were too poor quality and group doesn't exist anymore)
            foreach (var proteinGroup in proteinGroupsToRemove)
            {
                proteinGroups.Remove(proteinGroup);
            }

            // build razor peptide list (peptides that have >1 protein groups in the final, scored protein group list)
            foreach (var kvp in peptideToProteinGroupMatching)
            {
                if (kvp.Value.Count > 1)
                {
                    allRazorPeptides.Add(kvp.Key);
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                foreach (var peptide in proteinGroup.TotalPeptideList)
                {
                    // build razor peptide list for each protein group
                    if (allRazorPeptides.Contains(peptide))
                    {
                        // TODO**
                        // if the razor pep is associated with >1 protein group, it's a razor only for the group with the most ID'd peptides
                        // var sortedProteinGroups =
                        proteinGroup.StrictRazorPeptideList.Add(peptide);
                    }
                }

                // calculate sequence coverage for each protein in the group
                proteinGroup.CalculateSequenceCoverage();
            }
        }

        public List<ProteinGroup> DoProteinFdr(List<ProteinGroup> proteinGroups)
        {
            Status("Calculating protein FDR...", nestedIds);

            // order protein groups by score
            var sortedProteinGroups = proteinGroups.OrderByDescending(b => b.proteinGroupScore).ToList();
            List<ProteinGroup> proteinGroupsToRemove = new List<ProteinGroup>();

            // do fdr
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;
            foreach (var proteinGroup in sortedProteinGroups)
            {
                if (proteinGroup.isDecoy)
                {
                    cumulativeDecoy++;
                }
                else
                {
                    cumulativeTarget++;
                }

                proteinGroup.cumulativeTarget = cumulativeTarget;
                proteinGroup.cumulativeDecoy = cumulativeDecoy;
                proteinGroup.QValue = ((double)cumulativeDecoy / (cumulativeTarget + cumulativeDecoy));

                // remove one hit wonders if option enabled
                if (noOneHitWonders && !proteinGroup.isDecoy && !proteinGroup.isContaminant && proteinGroup.StrictPeptideList.Count == 1)
                    proteinGroupsToRemove.Add(proteinGroup);
            }

            foreach (var pg in proteinGroupsToRemove)
            {
                sortedProteinGroups.Remove(pg);
            }

            return sortedProteinGroups;
        }

        public void RunQuantification(List<NewPsmWithFdr> psms, double rtTolerance, double ppmTolerance)
        {
            // key is rough m/z (m/z rounded to 2nd decimal), KVP is value; key of this is the peak, value is RT
            var roughMzToPeakMatching = new Dictionary<double, List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>>();
            //var peakToScanMatching = new Dictionary<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>();
            var fullSeqToPsmMatching = psms.GroupBy(p => p.thisPSM.FullSequence).ToList();
            var thisPeptidesNormalizedIsotopicAbundances = new Dictionary<string, List<KeyValuePair<double, double>>>();

            var minChargeState = psms.Select(p => p.thisPSM.newPsm.scanPrecursorCharge).Min();
            var maxChargeState = psms.Select(p => p.thisPSM.newPsm.scanPrecursorCharge).Max();
            var chargeStates = Enumerable.Range(minChargeState, maxChargeState - 1);

            // build theoretical m/z bins
            foreach (var pepGrouping in fullSeqToPsmMatching)
            {
                var pepWithFormula = new Proteomics.Peptide(pepGrouping.First().thisPSM.BaseSequence);
                var pepIsotopicDistribution = IsotopicDistribution.GetDistribution(pepWithFormula.GetChemicalFormula(), 0.0001, 0.01);
                var modmass = pepGrouping.First().thisPSM.PeptideMonoisotopicMass - pepWithFormula.MonoisotopicMass;
                var masses = pepIsotopicDistribution.Masses.ToList();
                var abundances = pepIsotopicDistribution.Intensities.ToList();
                var maxAbundance = abundances.Max();
                var temp = new List<KeyValuePair<double, double>>();

                for(int i = 0; i < masses.Count; i++)
                {
                    // normalize abundances
                    abundances[i] = abundances[i] / maxAbundance;

                    masses[i] += modmass;

                    // add to isotopic list to check if abundance >20% normalized abundance
                    if(abundances[i] > 0.2)
                        temp.Add(new KeyValuePair<double, double> (masses[i], abundances[i]));
                }

                thisPeptidesNormalizedIsotopicAbundances.Add(pepGrouping.First().thisPSM.FullSequence, temp);

                var thisPeptidesMass = temp.Where(x => x.Value == 1).First().Key;

                //var thisPeptidesMass = pepGrouping.First().thisPSM.PeptideMonoisotopicMass;
                //thisPeptidesMass += ((int)Math.Round(pepWithFormula.ElementCountWithIsotopes("C") / 92.764, 0)) * 1.003355;
                //thisPeptidesMass += ((int)Math.Round(pepWithFormula.ElementCountWithIsotopes("H") / 8695.65, 0)) * 1.006277;
                //thisPeptidesMass += ((int)Math.Round(pepWithFormula.ElementCountWithIsotopes("N") / 271.739, 0)) * 0.997035;
                //thisPeptidesMass += ((int)Math.Round(pepWithFormula.ElementCountWithIsotopes("S") / 131.579, 0)) * 0.999386;
                //thisPeptidesMass += ((int)Math.Round(pepWithFormula.ElementCountWithIsotopes("O") / 2631.58, 0)) * 1.004217;
                //thisPeptidesMass += ((int)Math.Round(pepWithFormula.ElementCountWithIsotopes("S") / 23.31, 0)) * 1.995796;
                //thisPeptidesMass += ((int)Math.Round(pepWithFormula.ElementCountWithIsotopes("O") / 487.805, 0)) * 2.004245;

                foreach (var pep in pepGrouping)
                    pep.thisPSM.newPsm.mostAbundantMass = thisPeptidesMass;

                foreach (var chargeState in chargeStates)
                {
                    var t = Chemistry.ClassExtensions.ToMz(thisPeptidesMass, chargeState);
                    var m = Math.Round(t, 3);
                    if (!roughMzToPeakMatching.ContainsKey(m))
                        roughMzToPeakMatching.Add(m, new List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>());
                }
            }

            // populate m/z bins with observed m/z's
            var allMs1Scans = myMsDataFile.Where(s => s.MsnOrder == 1).ToList();
            foreach (var scan in allMs1Scans)
            {
                foreach (var peak in scan.MassSpectrum)
                {
                    if (peak.Intensity > 5000)
                    {
                        List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>> mzBin;
                        double floorMz = Math.Floor(peak.Mz * 1000) / 1000;
                        double ceilingMz = Math.Ceiling(peak.Mz * 1000) / 1000;

                        if (roughMzToPeakMatching.TryGetValue(floorMz, out mzBin))
                            mzBin.Add(new KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>(peak, scan));
                        if (roughMzToPeakMatching.TryGetValue(ceilingMz, out mzBin))
                            mzBin.Add(new KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>(peak, scan));
                    }
                }
            }

            // remove theoretical m/z bins not observed
            roughMzToPeakMatching = roughMzToPeakMatching.Where(x => x.Value.Count != 0).ToDictionary(x => x.Key, x => x.Value);

            // order bins by intensity, so the first match for each charge state will be the best
            roughMzToPeakMatching = roughMzToPeakMatching.ToDictionary(x => x.Key, x => x.Value.OrderByDescending(p => p.Key.Intensity).ToList());

            // find apex intensity of each peptide
            foreach (var pepGrouping in fullSeqToPsmMatching)
            {
                var verfiedPeaks = new List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>();
                
                // find peaks within specified tolerances in the m/z bins
                foreach (var chargeState in chargeStates)
                {
                    bool doneWithThisChargeState = false;
                    double theorMzHere = Chemistry.ClassExtensions.ToMz(pepGrouping.First().thisPSM.newPsm.mostAbundantMass, chargeState);
                    double mzTolHere = ((ppmTolerance / 1e6) * pepGrouping.First().thisPSM.newPsm.mostAbundantMass) / chargeState;

                    double floorMz = Math.Floor(theorMzHere * 1000) / 1000;
                    double ceilingMz = Math.Ceiling(theorMzHere * 1000) / 1000;
                    
                    IEnumerable<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>> binPeaks = new List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>();
                    List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>> t;
                    if (roughMzToPeakMatching.TryGetValue(floorMz, out t))
                        binPeaks = binPeaks.Concat(t);
                    if (roughMzToPeakMatching.TryGetValue(ceilingMz, out t))
                        binPeaks = binPeaks.Concat(t);

                    foreach (var peak in binPeaks)
                    {
                        // check ppm tolerance
                        if (Math.Abs(peak.Key.Mz - theorMzHere) < mzTolHere)
                        {
                            // check rt window tolerance
                            foreach (var pep in pepGrouping)
                            {
                                if (Math.Abs(pep.thisPSM.newPsm.scanRetentionTime - peak.Value.RetentionTime) < rtTolerance)
                                {
                                    // TODO** check isotopic distribution
                                    // TODO** check SNR
                                    // peak has been verified - add it
                                    if (!verfiedPeaks.Contains(peak))
                                    {
                                        verfiedPeaks.Add(peak);
                                        doneWithThisChargeState = true;
                                        break;
                                    }
                                }
                            }
                            if (doneWithThisChargeState)
                                break;
                        }
                    }
                }
            
                double apexIntensity = 0;
                double apexRT = 0;
                double apexMZ = 0;

                if (verfiedPeaks.Any())
                {
                    apexIntensity = verfiedPeaks.Select(x => x.Key.Intensity).Max();
                    KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>> p = verfiedPeaks.Where(x => x.Key.Intensity == apexIntensity).First();
                    apexRT = p.Value.RetentionTime;
                    apexMZ = p.Key.Mz;
                }

                foreach (var pep in pepGrouping)
                {
                    pep.thisPSM.newPsm.apexIntensity = apexIntensity;
                    pep.thisPSM.newPsm.apexRT = apexRT;
                    //pep.thisPSM.newPsm.apexMZ = apexMZ;
                }
            }

            // TODO** error checking for peptides that use the same apex peak
            // assign peptide with closest (ms2RT-apexRT) to apex, try to find other peaks
            /*
            var v = new Dictionary<double, IGrouping<string, NewPsmWithFdr>>();
            IGrouping<string, NewPsmWithFdr> v3;
            var badMatchList = new List<IGrouping<string, NewPsmWithFdr>>();
            foreach(var pepGrouping in fullSeqToPsmMatching)
            {
                if (v.TryGetValue(pepGrouping.First().thisPSM.newPsm.apexRT, out v3))
                {
                    badMatchList.Add(v3);
                    badMatchList.Add(pepGrouping);
                }
                else
                    v.Add(pepGrouping.First().thisPSM.newPsm.apexRT, pepGrouping);
            }

            foreach(var badMatch in badMatchList)
            {
                foreach(var pep in badMatch)
                {
                    pep.thisPSM.newPsm.apexIntensity = 0;
                    pep.thisPSM.newPsm.apexRT = 0;
                }
            }
            */
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            AnalysisResults myAnalysisResults = new AnalysisResults(this);
            Status("Running analysis engine!", nestedIds);
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy

            #region Match Seqeunces to PeptideWithSetModifications

            myAnalysisResults.AddText("Starting compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            Status("Adding observed peptides to dictionary...", nestedIds);
            foreach (var psmListForAspecificSerchMode in newPsms)
                if (psmListForAspecificSerchMode != null)
                    foreach (var psm in psmListForAspecificSerchMode)
                        if (psm != null)
                        {
                            var cp = psm.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications);
                            if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                                compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                        }
            myAnalysisResults.AddText("Ending compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            int totalProteins = proteinList.Count;
            int proteinsSeen = 0;
            int old_progress = 0;
            var obj = new object();
            Status("Adding possible sources to peptide dictionary...", nestedIds);
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                for (int i = fff.Item1; i < fff.Item2; i++)
                    foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maximumMissedCleavages, initiatorMethionineBehavior, fixedModifications))
                    {
                        if (peptideWithPossibleModifications.Length <= 1)
                            continue;
                        foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModIsoforms, max_mods_for_peptide))
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications, variableModifications, localizeableModifications, fixedModifications), out v))
                                v.Add(peptideWithSetModifications);
                        }
                    }
                lock (obj)
                {
                    foreach (var ye in local)
                    {
                        HashSet<PeptideWithSetModifications> v;
                        if (compactPeptideToProteinPeptideMatching.TryGetValue(ye.Key, out v))
                            foreach (var huh in ye.Value)
                                v.Add(huh);
                    }
                    proteinsSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In adding possible sources to peptide dictionary loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });

            #endregion Match Seqeunces to PeptideWithSetModifications

            List<ProteinGroup>[] proteinGroups = null;
            if (doParsimony)
            {
                Status("Applying protein parsimony...", nestedIds);
                proteinGroups = new List<ProteinGroup>[searchModes.Count];
                // TODO**: make this faster (only apply parsimony once but make multiple instances of the same ProteinGroups
                for (int i = 0; i < searchModes.Count; i++)
                    ApplyProteinParsimony(out proteinGroups[i]);
            }

            List<NewPsmWithFdr>[] allResultingIdentifications = new List<NewPsmWithFdr>[searchModes.Count];

            for (int j = 0; j < searchModes.Count; j++)
            {
                if (newPsms[j] != null)
                {
                    PsmWithMultiplePossiblePeptides[] psmsWithProteinHashSet = new PsmWithMultiplePossiblePeptides[newPsms[0].Length];
                    for (int i = 0; i < newPsms[0].Length; i++)
                    {
                        var huh = newPsms[j][i];
                        if (huh != null && huh.score >= 1)
                            psmsWithProteinHashSet[i] = new PsmWithMultiplePossiblePeptides(huh, compactPeptideToProteinPeptideMatching[huh.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications)], fragmentTolerance, myMsDataFile, lp);
                    }

                    var orderedPsmsWithPeptides = psmsWithProteinHashSet.Where(b => b != null).OrderByDescending(b => b.Score);

                    Status("Running FDR analysis...", nestedIds);
                    var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides, searchModes[j]);

                    if (quantify)
                    {
                        Status("Quantifying peptides...", nestedIds);
                        RunQuantification(orderedPsmsWithFDR, quantifyRtTol, quantifyPpmTol);
                    }

                    writePsmsAction(orderedPsmsWithFDR, searchModes[j].FileNameAddition);

                    if (doHistogramAnalysis)
                    {
                        var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.qValue <= 0.01)).ToList();
                        if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                        {
                            Status("Running histogram analysis...", nestedIds);
                            var myTreeStructure = new BinTreeStructure();
                            myTreeStructure.GenerateBins(limitedpsms_with_fdr, binTol);
                            writeHistogramPeaksAction(myTreeStructure, searchModes[j].FileNameAddition);
                        }
                    }

                    else
                    {
                        Status("Running FDR analysis on unique peptides...", nestedIds);
                        writePsmsAction(DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()), searchModes[j]), "uniquePeptides" + searchModes[j].FileNameAddition);
                    }

                    if (doParsimony)
                    {
                        ScoreProteinGroups(proteinGroups[j], orderedPsmsWithFDR);
                        proteinGroups[j] = DoProteinFdr(proteinGroups[j]);
                        writeProteinGroupsAction(proteinGroups[j], searchModes[j].FileNameAddition);
                    }

                    allResultingIdentifications[j] = orderedPsmsWithFDR;
                }
            }

            myAnalysisResults.AllResultingIdentifications = allResultingIdentifications;
            myAnalysisResults.ProteinGroups = proteinGroups;
            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<NewPsmWithFdr> DoFalseDiscoveryRateAnalysis(IEnumerable<PsmWithMultiplePossiblePeptides> items, SearchMode sm)
        {
            var ids = new List<NewPsmWithFdr>();
            foreach (PsmWithMultiplePossiblePeptides item in items)
                ids.Add(new NewPsmWithFdr(item));

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches];

            for (int i = 0; i < ids.Count; i++)
            {
                var item = ids[i];
                var isDecoy = item.IsDecoy;
                int notch = item.thisPSM.newPsm.notch;
                if (isDecoy)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                if (isDecoy)
                    cumulative_decoy_per_notch[notch]++;
                else
                    cumulative_target_per_notch[notch]++;

                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                double temp_q_value_for_notch = (double)cumulative_decoy_per_notch[notch] / (cumulative_target_per_notch[notch] + cumulative_decoy_per_notch[notch]);
                item.SetValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch);
            }

            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[sm.NumNotches];
            for (int i = 0; i < sm.NumNotches; i++)
                min_q_value_notch[i] = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                NewPsmWithFdr id = ids[i];
                if (id.qValue > min_q_value)
                    id.qValue = min_q_value;
                else if (id.qValue < min_q_value)
                    min_q_value = id.qValue;

                int notch = id.thisPSM.newPsm.notch;
                if (id.qValueNotch > min_q_value_notch[notch])
                    id.qValueNotch = min_q_value_notch[notch];
                else if (id.qValueNotch < min_q_value_notch[notch])
                    min_q_value_notch[notch] = id.qValueNotch;
            }

            return ids;
        }

        #endregion Private Methods

    }
}