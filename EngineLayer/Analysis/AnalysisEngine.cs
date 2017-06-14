using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Analysis
{
    public class AnalysisEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private const int max_mods_for_peptide = 3;
        private readonly double binTol;
        private readonly int maximumMissedCleavages;
        private readonly int? minPeptideLength;
        private readonly int? maxPeptideLength;
        private readonly int maxModIsoforms;
        private readonly PsmParent[][] newPsms;
        private readonly List<Protein> proteinList;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly List<ModificationWithMass> fixedModifications;
        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;
        private readonly Protease protease;
        private readonly List<SearchMode> searchModes;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private readonly Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans;
        private readonly Tolerance fragmentTolerance;
        private readonly Action<BinTreeStructure, string> writeHistogramPeaksAction;
        private readonly Action<List<NewPsmWithFdr>, string, List<string>> writePsmsAction;
        private readonly Action<List<ProteinGroup>, string, List<string>> writeProteinGroupsAction;
        private readonly Action<List<NewPsmWithFdr>, List<ProteinGroup>, SearchMode, string, List<string>> writeMzIdentmlAction;
        private readonly bool doParsimony;
        private readonly bool noOneHitWonders;
        private readonly bool doHistogramAnalysis;
        private readonly bool treatModPeptidesAsDifferentPeptides;
        private readonly bool quantify;
        private readonly double quantifyRtTol;
        private readonly double quantifyPpmTol;
        private readonly List<ProductType> lp;
        private readonly InitiatorMethionineBehavior initiatorMethionineBehavior;
        private readonly List<string> nestedIds;
        private Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;

        #endregion Private Fields

        #region Public Constructors

        public AnalysisEngine(PsmParent[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, Protease protease, List<SearchMode> searchModes, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFdr>, string, List<string>> action2, Action<List<ProteinGroup>, string, List<string>> action3, Action<List<NewPsmWithFdr>, List<ProteinGroup>, SearchMode, string, List<string>> action4,  bool doParsimony, bool noOneHitWonders, bool modPeptidesAreUnique, int maximumMissedCleavages, int? minPeptideLength, int? maxPeptideLength, int maxModIsoforms, bool doHistogramAnalysis, List<ProductType> lp, double binTol, InitiatorMethionineBehavior initiatorMethionineBehavior, List<string> nestedIds, bool Quantify, double QuantifyRtTol, double QuantifyPpmTol, Dictionary<ModificationWithMass, ushort> modsDictionary, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile)
        {
            this.doParsimony = doParsimony;
            this.noOneHitWonders = noOneHitWonders;
            this.doHistogramAnalysis = doHistogramAnalysis;
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.protease = protease;
            this.searchModes = searchModes;
            this.myMsDataFile = myMSDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.writeHistogramPeaksAction = action1;
            this.writePsmsAction = action2;
            this.writeProteinGroupsAction = action3;
            this.writeMzIdentmlAction = action4;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.minPeptideLength = minPeptideLength;
            this.maxPeptideLength = maxPeptideLength;
            this.maxModIsoforms = maxModIsoforms;
            this.lp = lp;
            this.binTol = binTol;
            this.initiatorMethionineBehavior = initiatorMethionineBehavior;
            this.quantify = Quantify;
            this.quantifyRtTol = QuantifyRtTol;
            this.quantifyPpmTol = QuantifyPpmTol;
            this.nestedIds = nestedIds;
            this.treatModPeptidesAsDifferentPeptides = modPeptidesAreUnique;
            this.modsDictionary = modsDictionary;
            this.arrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
        }

        #endregion Public Constructors

        #region Public Methods

        public void ApplyProteinParsimony(out List<ProteinGroup> proteinGroups)
        {
            // digesting an XML database results in a non-mod-agnostic digestion; need to fix this if mod-agnostic parsimony enabled
            if (!treatModPeptidesAsDifferentPeptides)
            {
                var groups = compactPeptideToProteinPeptideMatching.GroupBy(p => p.Value.First().BaseSequence);
                
                foreach(var baseseq in groups)
                {
                    var proteinsForAllCompactPeptidesWithThisBaseSeq = baseseq.SelectMany(p => p.Value.Select(v => v.Protein)).Distinct().ToList();

                    foreach (var compactPeptide in baseseq)
                    {
                        var proteinsHere = compactPeptide.Value.Select(p => p.Protein).Distinct();

                        if (proteinsHere.Count() != proteinsForAllCompactPeptidesWithThisBaseSeq.Count)
                        {
                            var missingProteins = proteinsForAllCompactPeptidesWithThisBaseSeq.Except(proteinsHere);

                            foreach (var missingProtein in missingProteins)
                            {
                                var templatePepWithSetMod = baseseq.SelectMany(p => p.Value).Where(v => v.Protein.Equals(missingProtein)).First();
                                compactPeptideToProteinPeptideMatching[compactPeptide.Key].Add(new PeptideWithSetModifications(compactPeptide.Value.First(), templatePepWithSetMod));
                            }
                        }
                    }
                }
            }

            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var parsimonyProteinList = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var proteinsWithUniquePeptides = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            // peptide matched to fullseq (used depending on user preference)
            var compactPeptideToFullSeqMatch = compactPeptideToProteinPeptideMatching.ToDictionary(x => x.Key, x => x.Value.First().Sequence);

            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                // finds unique peptides (peptides that can belong to only one protein)
                HashSet<Protein> proteinsAssociatedWithThisPeptide = new HashSet<Protein>(kvp.Value.Select(p => p.Protein));
                if (proteinsAssociatedWithThisPeptide.Count() == 1)
                {
                    var peptides = new HashSet<PeptideWithSetModifications>();
                    if (!proteinsWithUniquePeptides.TryGetValue(kvp.Value.First().Protein, out peptides))
                        proteinsWithUniquePeptides.Add(kvp.Value.First().Protein, new HashSet<PeptideWithSetModifications>(kvp.Value));
                    else
                        peptides.UnionWith(kvp.Value);
                }

                // if a peptide is associated with a decoy protein, remove all target protein associations with the peptide
                if (kvp.Value.Where(p => p.Protein.IsDecoy).Any())
                    kvp.Value.RemoveWhere(p => !p.Protein.IsDecoy);

                // if a peptide is associated with a contaminant protein, remove all target protein associations with the peptide
                if (kvp.Value.Where(p => p.Protein.IsContaminant).Any())
                    kvp.Value.RemoveWhere(p => !p.Protein.IsContaminant);
            }

            // makes dictionary with proteins as keys and list of associated peptides as the value (makes parsimony algo easier)
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    var peptides = new HashSet<CompactPeptide>();
                    if (!proteinToPeptidesMatching.TryGetValue(peptide.Protein, out peptides))
                        proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<CompactPeptide>() { kvp.Key });
                    else
                        peptides.Add(kvp.Key);
                }
            }

            // build protein list for each peptide before parsimony has been applied
            var peptideSeqProteinListMatch = new Dictionary<string, HashSet<Protein>>();
            foreach (var kvp in proteinToPeptidesMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    var proteinListHere = new HashSet<Protein>();
                    string pepSequence;
                    if (!treatModPeptidesAsDifferentPeptides)
                        pepSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));
                    else
                        pepSequence = compactPeptideToFullSeqMatch[peptide];

                    if (!peptideSeqProteinListMatch.TryGetValue(pepSequence, out proteinListHere))
                        peptideSeqProteinListMatch.Add(pepSequence, new HashSet<Protein>() { kvp.Key });
                    else
                        proteinListHere.Add(kvp.Key);
                }
            }

            // dictionary associates proteins w/ unused base seqs (list will shrink over time)
            var algDictionary = new Dictionary<Protein, HashSet<string>>();
            foreach (var kvp in peptideSeqProteinListMatch)
            {
                foreach (var protein in kvp.Value)
                {
                    HashSet<string> newPeptideBaseSeqs;
                    if (algDictionary.TryGetValue(protein, out newPeptideBaseSeqs))
                        newPeptideBaseSeqs.Add(kvp.Key);
                    else
                        algDictionary.Add(protein, new HashSet<string> { kvp.Key });
                }
            }

            // dictionary associates proteins w/ unused base seqs (list will NOT shrink over time)
            var proteinToPepSeqMatch = algDictionary.ToDictionary(x => x.Key, x => x.Value);

            // *** main parsimony loop
            bool uniquePeptidesLeft = false;
            if (proteinsWithUniquePeptides.Any())
                uniquePeptidesLeft = true;
            int numNewSeqs = algDictionary.Max(p => p.Value.Count);

            while (numNewSeqs != 0)
            {
                var possibleBestProteinList = new List<KeyValuePair<Protein, HashSet<string>>>();

                if (uniquePeptidesLeft)
                {
                    var proteinsWithUniquePeptidesLeft = algDictionary.Where(p => proteinsWithUniquePeptides.ContainsKey(p.Key));
                    if (proteinsWithUniquePeptidesLeft.Any())
                        possibleBestProteinList.Add(proteinsWithUniquePeptidesLeft.First());
                    else
                        uniquePeptidesLeft = false;
                }

                // gets list of proteins with the most unaccounted-for peptide base sequences
                if (!uniquePeptidesLeft)
                    possibleBestProteinList = algDictionary.Where(p => p.Value.Count == numNewSeqs).ToList();

                Protein bestProtein = possibleBestProteinList.First().Key;
                HashSet<string> newSeqs = new HashSet<string>(algDictionary[bestProtein]);

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
                        var proteinsOrderedByTotalPeptideCount = new Dictionary<Protein, HashSet<string>>();
                        foreach (var protein in proteinsWithTheseBaseSeqs)
                            proteinsOrderedByTotalPeptideCount.Add(protein, proteinToPepSeqMatch[protein]);

                        bestProtein = proteinsOrderedByTotalPeptideCount.OrderByDescending(kvp => kvp.Value.Count).First().Key;
                    }
                }

                parsimonyProteinList.Add(bestProtein, proteinToPeptidesMatching[bestProtein]);

                // remove used peptides from their proteins
                foreach (var newBaseSeq in newSeqs)
                {
                    HashSet<Protein> proteinsWithThisPeptide = peptideSeqProteinListMatch[newBaseSeq];

                    foreach (var protein in proteinsWithThisPeptide)
                        algDictionary[protein].Remove(newBaseSeq);
                }

                algDictionary.Remove(bestProtein);
                if (algDictionary.Any())
                    numNewSeqs = algDictionary.Max(p => p.Value.Count);
                else
                    numNewSeqs = 0;
            }

            // *** done with parsimony

            // add indistinguishable proteins
            var proteinsGroupedByNumPeptides = proteinToPeptidesMatching.GroupBy(p => p.Value.Count);
            var parsimonyProteinsGroupedByNumPeptides = parsimonyProteinList.GroupBy(p => p.Value.Count);
            var indistinguishableProteins = new Dictionary<Protein, HashSet<CompactPeptide>>();

            foreach (var group in proteinsGroupedByNumPeptides)
            {
                var parsimonyProteinsWithSameNumPeptides = parsimonyProteinsGroupedByNumPeptides.Where(p => p.Key == group.Key).FirstOrDefault();
                if (parsimonyProteinsWithSameNumPeptides != null)
                {
                    foreach (var protein in group)
                    {
                        foreach (var parsimonyProteinWithThisNumPeptides in parsimonyProteinsWithSameNumPeptides)
                        {
                            if (parsimonyProteinWithThisNumPeptides.Key != protein.Key)
                                if (proteinToPeptidesMatching[parsimonyProteinWithThisNumPeptides.Key].SetEquals(proteinToPeptidesMatching[protein.Key]))
                                    indistinguishableProteins.Add(protein.Key, proteinToPeptidesMatching[protein.Key]);
                        }
                    }
                }
            }
            foreach (var protein in indistinguishableProteins)
                parsimonyProteinList.Add(protein.Key, protein.Value);

            foreach (var kvp in compactPeptideToProteinPeptideMatching)
                kvp.Value.RemoveWhere(p => !parsimonyProteinList.ContainsKey(p.Protein));

            proteinGroups = ConstructProteinGroups(new HashSet<PeptideWithSetModifications>(proteinsWithUniquePeptides.Values.SelectMany(p => p)), new HashSet<PeptideWithSetModifications>(compactPeptideToProteinPeptideMatching.Values.SelectMany(p => p)));

            Status("Finished Parsimony", nestedIds);
        }

        public List<ProteinGroup> ConstructProteinGroups(HashSet<PeptideWithSetModifications> uniquePeptides, HashSet<PeptideWithSetModifications> allPeptides)
        {
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            foreach (var peptide in allPeptides)
            {
                HashSet<PeptideWithSetModifications> peptidesHere;
                if (proteinToPeptidesMatching.TryGetValue(peptide.Protein, out peptidesHere))
                    peptidesHere.Add(peptide);
                else
                    proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<PeptideWithSetModifications> { peptide });
            }

            // build protein group after parsimony (each group only has 1 protein at this point, indistinguishables will be merged during scoring)
            foreach (var kvp in proteinToPeptidesMatching)
            {
                var allPeptidesHere = proteinToPeptidesMatching[kvp.Key];
                var uniquePeptidesHere = new HashSet<PeptideWithSetModifications>(allPeptidesHere.Where(p => uniquePeptides.Contains(p)));

                proteinGroups.Add(new ProteinGroup(new HashSet<Protein>() { kvp.Key }, allPeptidesHere, uniquePeptidesHere));
            }

            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.AllPeptides.RemoveWhere(p => !proteinGroup.Proteins.Contains(p.Protein));
                proteinGroup.DisplayModsOnPeptides = treatModPeptidesAsDifferentPeptides;
            }

            return proteinGroups;
        }

        public void ScoreProteinGroups(List<ProteinGroup> proteinGroups, List<NewPsmWithFdr> psmList)
        {
            Status("Scoring protein groups...", nestedIds);

            // add each protein groups PSMs
            var peptideToPsmMatching = new Dictionary<PeptideWithSetModifications, HashSet<NewPsmWithFdr>>();
            foreach (var psm in psmList)
            {
                if (psm.QValue <= 0.01)
                {
                    foreach (var pepWithSetMods in psm.thisPSM.Pli.PeptidesWithSetModifications)
                    {
                        var psmsForThisPeptide = new HashSet<NewPsmWithFdr>();

                        if (!peptideToPsmMatching.TryGetValue(pepWithSetMods, out psmsForThisPeptide))
                            peptideToPsmMatching.Add(pepWithSetMods, new HashSet<NewPsmWithFdr> { psm });
                        else
                            psmsForThisPeptide.Add(psm);
                    }
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                List<PeptideWithSetModifications> pepsToRemove = new List<PeptideWithSetModifications>();
                foreach (var peptide in proteinGroup.AllPeptides)
                {
                    // build PSM list for scoring
                    HashSet<NewPsmWithFdr> psms;
                    if (peptideToPsmMatching.TryGetValue(peptide, out psms))
                        proteinGroup.AllPsmsBelowOnePercentFDR.UnionWith(psms);
                    else
                        pepsToRemove.Add(peptide);
                }

                proteinGroup.AllPeptides.ExceptWith(pepsToRemove);
                proteinGroup.UniquePeptides.ExceptWith(pepsToRemove);
            }

            // score the group
            foreach (var proteinGroup in proteinGroups)
                proteinGroup.Score();

            // merge protein groups that are indistinguishable after scoring
            var pg = proteinGroups.OrderByDescending(p => p.ProteinGroupScore).ToList();
            for (int i = 0; i < (pg.Count - 1); i++)
            {
                if (pg[i].ProteinGroupScore == pg[i + 1].ProteinGroupScore && pg[i].ProteinGroupScore != 0)
                {
                    var pgsWithThisScore = pg.Where(p => p.ProteinGroupScore == pg[i].ProteinGroupScore).ToList();

                    // check to make sure they have the same peptides, then merge them
                    foreach (var p in pgsWithThisScore)
                    {
                        var seqs1 = new HashSet<string>(p.AllPeptides.Select(x => x.Sequence));
                        var seqs2 = new HashSet<string>(pg[i].AllPeptides.Select(x => x.Sequence));

                        if (p != pg[i] && seqs1.SetEquals(seqs2))
                        {
                            pg[i].MergeProteinGroupWith(p);
                        }
                    }
                }
            }

            // remove empty protein groups (peptides were too poor quality or group was merged)
            proteinGroups.RemoveAll(p => p.ProteinGroupScore == 0);

            // calculate sequence coverage
            foreach (var proteinGroup in proteinGroups)
                proteinGroup.CalculateSequenceCoverage();

            // distribute razor peptides
            var sharedPepWithProteinGroups = new Dictionary<string, HashSet<ProteinGroup>>();
            var seqToPeptide = new Dictionary<string, HashSet<PeptideWithSetModifications>>();
            foreach (var proteinGroup in proteinGroups)
            {
                var sharedPeps = proteinGroup.AllPeptides.Except(proteinGroup.UniquePeptides);
                foreach (var sharedPep in sharedPeps)
                {
                    HashSet<ProteinGroup> v;
                    if (sharedPepWithProteinGroups.TryGetValue(sharedPep.Sequence, out v))
                        v.Add(proteinGroup);
                    else
                        sharedPepWithProteinGroups.Add(sharedPep.Sequence, new HashSet<ProteinGroup> { proteinGroup });

                    HashSet<PeptideWithSetModifications> v1;
                    if (seqToPeptide.TryGetValue(sharedPep.Sequence, out v1))
                        v1.Add(sharedPep);
                    else
                        seqToPeptide.Add(sharedPep.Sequence, new HashSet<PeptideWithSetModifications> { sharedPep });
                }
            }

            foreach (var kvp in sharedPepWithProteinGroups)
            {
                int i = kvp.Value.Select(p => p.AllPeptides.Select(x => x.BaseSequence).Count()).Max();
                HashSet<ProteinGroup> t = new HashSet<ProteinGroup>(kvp.Value.Where(p => p.AllPeptides.Select(x => x.BaseSequence).Count() == i));
                foreach (var proteinGroup in t)
                {
                    var peptides = seqToPeptide[kvp.Key];
                    proteinGroup.RazorPeptides.UnionWith(peptides.Where(p => proteinGroup.Proteins.Contains(p.Protein)));
                }
            }

            foreach (var proteinGroup in proteinGroups)
                proteinGroup.Quantify();
        }

        public List<ProteinGroup> DoProteinFdr(List<ProteinGroup> proteinGroups)
        {
            Status("Calculating protein FDR...", nestedIds);

            if (noOneHitWonders)
            {
                if (treatModPeptidesAsDifferentPeptides)
                    proteinGroups = proteinGroups.Where(p => p.isDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.Sequence)).Count() > 1).ToList();
                else
                    proteinGroups = proteinGroups.Where(p => p.isDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.BaseSequence)).Count() > 1).ToList();
            }

            // order protein groups by score
            var sortedProteinGroups = proteinGroups.OrderByDescending(b => b.ProteinGroupScore).ToList();

            // do fdr
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;
            foreach (var proteinGroup in sortedProteinGroups)
            {
                if (proteinGroup.isDecoy)
                    cumulativeDecoy++;
                else
                    cumulativeTarget++;

                proteinGroup.CumulativeTarget = cumulativeTarget;
                proteinGroup.CumulativeDecoy = cumulativeDecoy;
                proteinGroup.QValue = ((double)cumulativeDecoy / (cumulativeTarget + cumulativeDecoy));
            }

            return sortedProteinGroups;
        }

        public void RunQuantification(List<NewPsmWithFdr> psms, double rtTolerance, double ppmTolerance)
        {
            // key is rough m/z (m/z rounded to 2nd decimal), value contains peak with its scan
            var mzBins = new Dictionary<double, List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>>();

            var peptideGroups = psms.GroupBy(p => p.thisPSM.Pli.FullSequence).ToList();
            var lengthToIsotopicDistribution = new Dictionary<int, List<KeyValuePair<double, double>>>();
            HashSet<int> peptideLengths = new HashSet<int>(peptideGroups.Select(p => p.First().thisPSM.Pli.BaseSequence.Length));

            foreach (var length in peptideLengths)
            {
                // averagine has 4.9384 carbon atoms
                var formula = "C" + (int)(4.9384 * length);
                var isotopicDistribution = IsotopicDistribution.GetDistribution(ChemicalFormula.ParseFormula(formula), 0.0001, 0.01);

                var masses = isotopicDistribution.Masses.ToArray();
                var abundances = isotopicDistribution.Intensities.ToArray();
                var massesWithAbundances = new List<KeyValuePair<double, double>>();

                var monoIsotopicMass = masses.Min();
                var highestAbundance = abundances.Max();

                for (int i = 0; i < masses.Length; i++)
                {
                    // expected isotopic mass shifts for peptide of this length
                    masses[i] -= monoIsotopicMass;

                    // normalized abundance of each mass shift
                    abundances[i] /= highestAbundance;

                    if (abundances[i] > 0.2)
                        massesWithAbundances.Add(new KeyValuePair<double, double>(masses[i], abundances[i]));
                }

                lengthToIsotopicDistribution.Add(length, massesWithAbundances);
            }

            var minChargeState = psms.Select(p => p.thisPSM.ScanPrecursorCharge).Min();
            var maxChargeState = psms.Select(p => p.thisPSM.ScanPrecursorCharge).Max();
            var chargeStates = Enumerable.Range(minChargeState, maxChargeState - 1);

            // build theoretical m/z bins
            foreach (var pepGrouping in peptideGroups)
            {
                var mostCommonIsotopeShift = lengthToIsotopicDistribution[pepGrouping.First().thisPSM.Pli.BaseSequence.Length].Where(p => p.Value == 1).First().Key;
                var thisPeptidesMass = pepGrouping.First().thisPSM.Pli.PeptideMonoisotopicMass + mostCommonIsotopeShift;

                foreach (var pep in pepGrouping)
                    pep.thisPSM.MostAbundantMass = thisPeptidesMass;

                foreach (var chargeState in chargeStates)
                {
                    var t = Chemistry.ClassExtensions.ToMz(thisPeptidesMass, chargeState);
                    var m = Math.Round(t, 2);
                    if (!mzBins.ContainsKey(m))
                        mzBins.Add(m, new List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>());
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
                        double floorMz = Math.Floor(peak.Mz * 100) / 100;
                        double ceilingMz = Math.Ceiling(peak.Mz * 100) / 100;

                        if (mzBins.TryGetValue(floorMz, out mzBin))
                            mzBin.Add(new KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>(peak, scan));
                        if (mzBins.TryGetValue(ceilingMz, out mzBin))
                            mzBin.Add(new KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>(peak, scan));
                    }
                }
            }

            // remove theoretical m/z bins not observed
            mzBins = mzBins.Where(x => x.Value.Count != 0).ToDictionary(x => x.Key, x => x.Value);

            // find apex intensity of each peptide
            foreach (var pepGrouping in peptideGroups)
            {
                var verfiedPeaks = new List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>();

                // find peaks within specified tolerances in the m/z bins
                foreach (var chargeState in chargeStates)
                {
                    double theorMzHere = Chemistry.ClassExtensions.ToMz(pepGrouping.First().thisPSM.MostAbundantMass, chargeState);
                    double mzTolHere = ((ppmTolerance / 1e6) * pepGrouping.First().thisPSM.MostAbundantMass) / chargeState;

                    double floorMz = Math.Floor(theorMzHere * 100) / 100;
                    double ceilingMz = Math.Ceiling(theorMzHere * 100) / 100;

                    IEnumerable<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>> binPeaks = new List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>>();
                    List<KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>> t;
                    if (mzBins.TryGetValue(floorMz, out t))
                        binPeaks = binPeaks.Concat(t);
                    if (mzBins.TryGetValue(ceilingMz, out t))
                        binPeaks = binPeaks.Concat(t);

                    binPeaks = binPeaks.OrderByDescending(p => p.Key.Intensity);

                    foreach (var peakWithScan in binPeaks)
                    {
                        if (!verfiedPeaks.Contains(peakWithScan))
                        {
                            // check ppm tolerance
                            if (Math.Abs(peakWithScan.Key.Mz - theorMzHere) < mzTolHere)
                            {
                                // check rt tolerance
                                var validRTs = pepGrouping.Select(p => p.thisPSM.ScanRetentionTime).Where(p => Math.Abs(peakWithScan.Value.RetentionTime - p) < rtTolerance);

                                if (validRTs.Any())
                                {
                                    // check isotopic distribution
                                    var temp = lengthToIsotopicDistribution[pepGrouping.First().thisPSM.Pli.BaseSequence.Length];
                                    var isotopes = temp.Select(p => new KeyValuePair<double, double>(p.Key + pepGrouping.First().thisPSM.Pli.PeptideMonoisotopicMass, p.Value)).ToList();

                                    var lowestMassIsotope = isotopes.Select(p => p.Key).Min();
                                    var highestMassIsotope = isotopes.Select(p => p.Key).Max();
                                    IMzPeak[] isotopePeaks = new IMzPeak[isotopes.Count];

                                    lowestMassIsotope = Chemistry.ClassExtensions.ToMz(lowestMassIsotope, chargeState);
                                    lowestMassIsotope -= (ppmTolerance / 1e6) * lowestMassIsotope;

                                    highestMassIsotope = Chemistry.ClassExtensions.ToMz(highestMassIsotope, chargeState);
                                    highestMassIsotope += (ppmTolerance / 1e6) * highestMassIsotope;

                                    var otherPeaksInThisScan = peakWithScan.Value.MassSpectrum.Where(p => p.Mz > lowestMassIsotope && p.Mz < highestMassIsotope).ToList();
                                    bool isotopeDistributionCheck = true;

                                    foreach (var isotope in isotopes)
                                    {
                                        double theorIsotopeMz = Chemistry.ClassExtensions.ToMz(isotope.Key, chargeState);
                                        double isotopeMzTol = ((3.0 / 1e6) * isotope.Key) / chargeState;
                                        bool thisIsotopeSeen = false;

                                        foreach (var otherPeak in otherPeaksInThisScan)
                                        {
                                            if (Math.Abs(otherPeak.Mz - theorIsotopeMz) < isotopeMzTol)
                                            {
                                                thisIsotopeSeen = true;
                                                isotopePeaks[isotopes.IndexOf(isotope)] = otherPeak;
                                                break;
                                            }
                                        }

                                        if (!thisIsotopeSeen)
                                            isotopeDistributionCheck = false;
                                    }

                                    if (isotopeDistributionCheck)
                                    {
                                        for (int i = 0; i < isotopes.Count; i++)
                                        {
                                            if (isotopes[i].Value < 0.8)
                                                isotopePeaks[i] = null;
                                        }
                                        isotopePeaks = isotopePeaks.Where(v => v != null).ToArray();
                                        double maxIsotopeIntensity = isotopePeaks.Select(v => v.Intensity).Max();
                                        IMzPeak maxIsotopicPeak = isotopePeaks.Where(x => x.Intensity == maxIsotopeIntensity).First();
                                        verfiedPeaks.Add(new KeyValuePair<IMzPeak, IMsDataScan<IMzSpectrum<IMzPeak>>>(maxIsotopicPeak, peakWithScan.Value));
                                        // done with this charge state, move on to the next one
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                double apexIntensity = 0;
                double apexRT = 0;
                double apexMZ = 0;

                if (verfiedPeaks.Any())
                {
                    apexIntensity = verfiedPeaks.Select(x => x.Key.Intensity).Max();
                    var p = verfiedPeaks.Where(x => x.Key.Intensity == apexIntensity).First();
                    apexRT = p.Value.RetentionTime;
                    apexMZ = p.Key.Mz;
                }

                foreach (var pep in pepGrouping)
                {
                    if (apexIntensity != 0)
                    {
                        pep.thisPSM.QuantIntensity = new double[] { apexIntensity };
                        pep.thisPSM.QuantRT = apexRT;
                        pep.thisPSM.ApexMz = apexMZ;
                    }
                    else
                    {
                        pep.thisPSM.QuantIntensity = new double[] { pep.thisPSM.ScanPrecursorMonoisotopicPeak.Intensity };
                        pep.thisPSM.QuantRT = pep.thisPSM.ScanRetentionTime;
                        pep.thisPSM.ApexMz = pep.thisPSM.ScanPrecursorMonoisotopicPeak.Mz;
                    }
                }
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
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
                            var cp = psm.GetCompactPeptide(modsDictionary);
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
                    foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications))
                    {
                        if (peptideWithPossibleModifications.Length <= 1)
                            continue;
                        foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModIsoforms, max_mods_for_peptide))
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications, modsDictionary), out v))
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
            Dictionary<string, int>[] allModsSeen = new Dictionary<string, int>[searchModes.Count];
            Dictionary<string, int>[] allModsOnPeptides = new Dictionary<string, int>[searchModes.Count];

            for (int j = 0; j < searchModes.Count; j++)
            {
                if (newPsms[j] != null)
                {
                    Status("Computing info about actual peptides with modifications...", nestedIds);
                    for (int myScanWithMassIndex = 0; myScanWithMassIndex < newPsms[0].Length; myScanWithMassIndex++)
                    {
                        var huh = newPsms[j][myScanWithMassIndex];
                        if (huh != null && huh.Pli == null)
                            huh.ComputeProteinLevelInfo(compactPeptideToProteinPeptideMatching, fragmentTolerance, arrayOfSortedMS2Scans?[huh.ScanIndex], lp, modsDictionary);
                    }

                    Status("Sorting and grouping psms..", nestedIds);
                    var orderedPsmsWithPeptides = newPsms[j].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.Pli.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, string>(b.FileName, b.ScanNumber, b.Pli.FullSequence)).Select(b => b.First());

                    Status("Running FDR analysis...", nestedIds);
                    var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides, searchModes[j]);

                    Status("Running modification analysis...", nestedIds);

                    Dictionary<string, int> modsSeen = new Dictionary<string, int>();
                    Dictionary<string, int> modsOnPeptides = new Dictionary<string, int>();

                    foreach (var highConfidencePSM in orderedPsmsWithFDR.Where(b => b.QValue <= 0.01 && !b.IsDecoy).GroupBy(b => b.thisPSM.Pli.PeptidesWithSetModifications.First().Sequence).Select(b => b.FirstOrDefault()))
                    {
                        var singlePeptide = highConfidencePSM.thisPSM.Pli.PeptidesWithSetModifications.First();
                        var modsIdentified = singlePeptide.allModsOneIsNterminus;
                        foreach (var modSeen in modsIdentified)
                        {
                            if (modsSeen.ContainsKey(modSeen.Value.id))
                                modsSeen[modSeen.Value.id]++;
                            else
                                modsSeen.Add(modSeen.Value.id, 1);
                        }
                        var modsInProtein = singlePeptide.Protein.OneBasedPossibleLocalizedModifications.Where(b => b.Key >= singlePeptide.OneBasedStartResidueInProtein && b.Key <= singlePeptide.OneBasedEndResidueInProtein).SelectMany(b => b.Value);
                        foreach (var modInProtein in modsInProtein)
                        {
                            if (modsOnPeptides.ContainsKey(modInProtein.id))
                                modsOnPeptides[modInProtein.id]++;
                            else
                                modsOnPeptides.Add(modInProtein.id, 1);
                        }
                    }
                    allModsSeen[j] = modsSeen;
                    allModsOnPeptides[j] = modsOnPeptides;

                    if (quantify)
                    {
                        Status("Quantifying peptides...", nestedIds);
                        RunQuantification(orderedPsmsWithFDR, quantifyRtTol, quantifyPpmTol);
                    }

                    writePsmsAction?.Invoke(orderedPsmsWithFDR, "allPSMs_" + searchModes[j].FileNameAddition, nestedIds);
                    if (!doParsimony && writeMzIdentmlAction != null) writeMzIdentmlAction.Invoke(orderedPsmsWithFDR, null, searchModes[j], searchModes[j].FileNameAddition, nestedIds);
                    if (doHistogramAnalysis)
                    {
                        var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.QValue <= 0.01)).ToList();
                        if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                        {
                            Status("Running histogram analysis...", nestedIds);
                            var myTreeStructure = new BinTreeStructure();
                            myTreeStructure.GenerateBins(limitedpsms_with_fdr, binTol);
                            writeHistogramPeaksAction(myTreeStructure, searchModes[j].FileNameAddition);
                        }
                    }
                    Status("Running FDR analysis on unique peptides...", nestedIds);
                    var peptidesWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides.GroupBy(b => b.Pli.FullSequence).Select(b => b.FirstOrDefault()), searchModes[j]);
                    List<NewPsmWithFdr> peptidesWithGoodFDR = new List<NewPsmWithFdr>();
                    Console.WriteLine(peptidesWithFDR[0].QValue);
                    foreach(NewPsmWithFdr PSM in peptidesWithFDR)
                    {
                        if(PSM.QValue <= .01)
                        peptidesWithGoodFDR.Add(PSM);
                    }
                    myAnalysisResults.AddText("Unique PSMS within 1% FDR: " + peptidesWithGoodFDR.Count);
                    writePsmsAction?.Invoke(peptidesWithFDR, "uniquePeptides_" + searchModes[j].FileNameAddition, nestedIds);
                    

                    // individual (for single-file search) or aggregate results
                    if (doParsimony && writeProteinGroupsAction != null)
                    {
                        ScoreProteinGroups(proteinGroups[j], orderedPsmsWithFDR);
                        proteinGroups[j] = DoProteinFdr(proteinGroups[j]);

                        // call multifile protein quantification helper function (need all the filenames to organize results properly)
                        var files = orderedPsmsWithFDR.Select(p => p.thisPSM.FileName).Distinct().ToList();
                        foreach (var pg in proteinGroups[j])
                            pg.AggregateQuantifyHelper(files);
                        writeProteinGroupsAction(proteinGroups[j], searchModes[j].FileNameAddition, nestedIds);
                        writeMzIdentmlAction(orderedPsmsWithFDR, proteinGroups[j], searchModes[j], searchModes[j].FileNameAddition, nestedIds);
                    }


                    // write individual file results based on aggregate results
                    if (arrayOfSortedMS2Scans == null && doParsimony)
                    {
                        var psmsGroupedByFilename = orderedPsmsWithFDR.GroupBy(p => p.thisPSM.FileName);

                        // individual psm files (with global psm fdr, global parsimony)
                        foreach (var group in psmsGroupedByFilename)
                        {
                            var fileName = System.IO.Path.GetFileNameWithoutExtension(group.First().thisPSM.FileName);
                            writePsmsAction(group.ToList(), fileName + "_allPSMS_" + searchModes[j].FileNameAddition, new List<string>(nestedIds.Concat(new List<string> { "Individual Searches", group.First().thisPSM.FileName })));
                        }

                        // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                        var fileNames = psmsGroupedByFilename.Select(p => p.Key).Distinct();
                        foreach (var file in fileNames)
                        {
                            var subsetProteinGroupsForThisFile = new List<ProteinGroup>();
                            foreach (var pg in proteinGroups[j])
                            {
                                var subsetPg = pg.ConstructSubsetProteinGroup(file);
                                subsetPg.Score();

                                if (subsetPg.ProteinGroupScore != 0)
                                {
                                    subsetPg.CalculateSequenceCoverage();
                                    subsetPg.Quantify();
                                    subsetProteinGroupsForThisFile.Add(subsetPg);
                                }
                            }

                            subsetProteinGroupsForThisFile = DoProteinFdr(subsetProteinGroupsForThisFile);
                            writeProteinGroupsAction(subsetProteinGroupsForThisFile, file + "_" + searchModes[j].FileNameAddition, new List<string>(nestedIds.Concat(new List<string> { "Individual Searches", file })));
                            writeMzIdentmlAction(psmsGroupedByFilename.Where(p => p.Key == file).SelectMany(g => g).ToList(), subsetProteinGroupsForThisFile, searchModes[j], file + "_" + searchModes[j].FileNameAddition, new List<string>(nestedIds.Concat(new List<string> { "Individual Searches", file })));
                        }
                    }


                    allResultingIdentifications[j] = orderedPsmsWithFDR;
                }
            }

            myAnalysisResults.AllResultingIdentifications = allResultingIdentifications;
            myAnalysisResults.ProteinGroups = proteinGroups;
            myAnalysisResults.allModsSeen = allModsSeen;
            myAnalysisResults.allModsOnPeptides = allModsOnPeptides;
            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<NewPsmWithFdr> DoFalseDiscoveryRateAnalysis(IEnumerable<PsmParent> items, SearchMode sm)
        {
            var ids = new List<NewPsmWithFdr>();
            foreach (PsmParent item in items)
                ids.Add(new NewPsmWithFdr(item));

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches];

            for (int i = 0; i < ids.Count; i++)
            {
                var item = ids[i];
                var isDecoy = item.IsDecoy;
                int notch = item.thisPSM.Notch;
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
                if (id.QValue > min_q_value)
                    id.QValue = min_q_value;
                else if (id.QValue < min_q_value)
                    min_q_value = id.QValue;

                int notch = id.thisPSM.Notch;
                if (id.QValueNotch > min_q_value_notch[notch])
                    id.QValueNotch = min_q_value_notch[notch];
                else if (id.QValueNotch < min_q_value_notch[notch])
                    min_q_value_notch[notch] = id.QValueNotch;
            }

            return ids;
        }

        #endregion Private Methods

    }
}
