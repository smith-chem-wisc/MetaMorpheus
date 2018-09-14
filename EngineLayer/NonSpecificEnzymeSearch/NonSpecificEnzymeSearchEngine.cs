using Chemistry;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    public class NonSpecificEnzymeSearchEngine : ModernSearchEngine
    {
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private readonly List<int>[] fragmentIndexPrecursor;
        private readonly int MinimumPeptideLength;

        public NonSpecificEnzymeSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex, List<int>[] fragmentIndex, List<int>[] fragmentIndexPrecursor, int currentPartition, CommonParameters CommonParameters, MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds) : base(globalPsms, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, CommonParameters, massDiffAcceptor, maximumMassThatFragmentIonScoreIsDoubled, nestedIds)
        {
            this.fragmentIndexPrecursor = fragmentIndexPrecursor;
            MinimumPeptideLength = commonParameters.DigestionParams.MinPeptideLength;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing nonspecific search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)commonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedMs2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, range =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                HashSet<int> idsOfPeptidesPossiblyObserved = new HashSet<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = ListOfSortedMs2Scans[i];

                    //get bins to add points to
                    List<int> allBinsToSearch = GetBinsToSearch(scan);

                    for (int j = 0; j < allBinsToSearch.Count; j++)
                    {
                        FragmentIndex[allBinsToSearch[j]].ForEach(id => scoringTable[id]++);
                    }

                    //populate ids of possibly observed with those containing allowed precursor masses
                    List<int> binsToSearch = new List<int>();
                    int obsPrecursorFloorMz = (int)Math.Floor(commonParameters.PrecursorMassTolerance.GetMinimumValue(scan.PrecursorMass) * FragmentBinsPerDalton);
                    int obsPrecursorCeilingMz = (int)Math.Ceiling(commonParameters.PrecursorMassTolerance.GetMaximumValue(scan.PrecursorMass) * FragmentBinsPerDalton);
                    for (int fragmentBin = obsPrecursorFloorMz; fragmentBin <= obsPrecursorCeilingMz; fragmentBin++)
                    {
                        binsToSearch.Add(fragmentBin);
                    }

                    foreach (ProductType pt in DissociationTypeCollection.ProductsFromDissociationType[commonParameters.DissociationType].Intersect(TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[commonParameters.FragmentationTerminus]).ToList())
                    {
                        //TODO: check that this is correct (Zach)
                        int binShift = (int)Math.Round((WaterMonoisotopicMass - DissociationTypeCollection.GetMassShiftFromProductType(pt)) * FragmentBinsPerDalton);//if unit test fails try subtracting water

                        for (int j = 0; j < binsToSearch.Count; j++)
                        {
                            int bin = binsToSearch[j] - binShift;
                            if (bin < FragmentIndex.Length && FragmentIndex[bin] != null)
                            {
                                FragmentIndex[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                            }
                        }
                    }

                    for (int j = 0; j < binsToSearch.Count; j++)
                    {
                        int bin = binsToSearch[j];
                        if (bin < fragmentIndexPrecursor.Length && fragmentIndexPrecursor[bin] != null)
                        {
                            fragmentIndexPrecursor[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                        }
                    }

                    // done with initial scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        PeptideSpectralMatch PSM = null;
                        int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]) + 1;
                        while (maxInitialScore > commonParameters.ScoreCutoff)
                        {
                            maxInitialScore--;
                            foreach (var id in idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore))
                            {
                                PeptideWithSetModifications peptide = PeptideIndex[id];
                                List<Product> peptideTheorProducts = new List<Product>();
                                lock (peptide)
                                {
                                    peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, commonParameters.FragmentationTerminus).ToList();
                                }

                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan.MassSpectrum, peptideTheorProducts, commonParameters);

                                if (commonParameters.AddCompIons)
                                {
                                    MzSpectrum complementarySpectrum = GenerateComplementarySpectrum(scan.TheScan.MassSpectrum, scan.PrecursorMass, commonParameters.DissociationType);
                                    matchedIons.AddRange(MatchFragmentIons(complementarySpectrum, peptideTheorProducts, commonParameters));
                                }

                                double thisScore = CalculatePeptideScore(scan.TheScan, matchedIons, MaxMassThatFragmentIonScoreIsDoubled);

                                Tuple<int, double> notchAndPrecursor = Accepts(scan.PrecursorMass, peptide, commonParameters.FragmentationTerminus, MassDiffAcceptor);
                                if (notchAndPrecursor.Item1 >= 0)
                                {
                                    if (PSM == null)
                                    {
                                        PSM = new PeptideSpectralMatch(peptide, notchAndPrecursor.Item1, thisScore, i, scan, commonParameters.DigestionParams, matchedIons);
                                    }
                                    else
                                    {
                                        PSM.AddOrReplace(peptide, thisScore, notchAndPrecursor.Item1, commonParameters.ReportAllAmbiguity, matchedIons);
                                    }
                                }
                            }

                            if (PSM != null) //if we have a match
                            {
                                List<(int notch, PeptideWithSetModifications pwsm)> originalPwsmsWithNotches = PSM.BestMatchingPeptideWithSetMods.ToList();
                                List<(int notch, PeptideWithSetModifications pwsm)> updatedPwsmsWithNotches = new List<(int notch, PeptideWithSetModifications pwsm)>();

                                foreach ((int notch, PeptideWithSetModifications pwsm) originalPwsmWithNotch in originalPwsmsWithNotches)
                                {
                                    PeptideWithSetModifications pwsm = originalPwsmWithNotch.pwsm;
                                    List<double> initialMasses = new List<double>();
                                    if (pwsm.AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod)) //get terminal mods
                                    {
                                        foreach (double nl in pep_n_term_variable_mod.NeutralLosses[commonParameters.DissociationType])
                                        {
                                            double monoisotopicMass = pep_n_term_variable_mod.MonoisotopicMass ?? 0;
                                            initialMasses.Add(monoisotopicMass - nl);
                                        }
                                    }
                                    else
                                    {
                                        initialMasses.Add(0);
                                    }

                                    foreach (double initialMass in initialMasses)
                                    {
                                        double finalMass = initialMass + WaterMonoisotopicMass;

                                        //generate correct sequence
                                        if (commonParameters.FragmentationTerminus == FragmentationTerminus.N)
                                        {
                                            int index = ComputePeptideIndexes(pwsm, ref finalMass, 1, 1, scan.PrecursorMass, MassDiffAcceptor);

                                            if (index >= 0 && index >= MinimumPeptideLength)
                                            {
                                                Dictionary<int, Modification> allModsOneIsNTerminus = pwsm.AllModsOneIsNterminus
                                                .Where(b => b.Key > 1 && b.Key <= (1 + index)).ToDictionary(b => b.Key, b => b.Value);

                                                PeptideWithSetModifications actualPWSM = new PeptideWithSetModifications(pwsm.Protein, commonParameters.DigestionParams, pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedStartResidueInProtein + index - 1, pwsm.PeptideDescription, pwsm.MissedCleavages, allModsOneIsNTerminus, pwsm.NumFixedMods);
                                                updatedPwsmsWithNotches.Add((originalPwsmWithNotch.notch, actualPWSM));
                                                break;
                                            }
                                        }
                                        else //if C terminus
                                        {
                                            int index = ComputePeptideIndexes(pwsm, ref finalMass, pwsm.Length, -1, scan.PrecursorMass, MassDiffAcceptor);

                                            if (index >= 0 && (pwsm.OneBasedEndResidueInProtein - (pwsm.OneBasedStartResidueInProtein + index - 2)) >= MinimumPeptideLength)
                                            {
                                                Dictionary<int, Modification> allModsOneIsNTerminus = pwsm.AllModsOneIsNterminus
                                                .Where(b => b.Key > index && b.Key <= (2 + pwsm.OneBasedEndResidueInProtein - pwsm.OneBasedStartResidueInProtein)).ToDictionary(b => (b.Key + index - 1), b => b.Value);

                                                PeptideWithSetModifications actualPWSM = new PeptideWithSetModifications(pwsm.Protein, commonParameters.DigestionParams, pwsm.OneBasedStartResidueInProtein + index - 1, pwsm.OneBasedEndResidueInProtein, pwsm.PeptideDescription, pwsm.MissedCleavages, allModsOneIsNTerminus, pwsm.NumFixedMods);
                                                updatedPwsmsWithNotches.Add((originalPwsmWithNotch.notch, actualPWSM));
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (updatedPwsmsWithNotches.Count != 0) //need a unit test for this, such as ECD (both y and zdot), where a y mass hits a zdot shift but not a ydot shift
                                {
                                    var firstPwsmWithNotch = updatedPwsmsWithNotches.First();
                                    PeptideSpectralMatch updatedPSM = new PeptideSpectralMatch(firstPwsmWithNotch.pwsm, firstPwsmWithNotch.notch, PSM.Score, i, scan, commonParameters.DigestionParams, PSM.PeptidesToMatchingFragments[originalPwsmsWithNotches.First().pwsm]);
                                    for (int pwsmIndex = 1; pwsmIndex < updatedPwsmsWithNotches.Count; pwsmIndex++)
                                    {
                                        var currentPwsmWithNotch = updatedPwsmsWithNotches[pwsmIndex];

                                        //TODO: This is unnecesary, should be able to back cleave fragments
                                        List<Product> peptideTheorProducts = new List<Product>();
                                        lock (currentPwsmWithNotch.pwsm)
                                        {
                                            peptideTheorProducts = currentPwsmWithNotch.pwsm.Fragment(commonParameters.DissociationType, commonParameters.FragmentationTerminus).ToList();
                                        }
                                        List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan.MassSpectrum, peptideTheorProducts, commonParameters);

                                        if (commonParameters.AddCompIons)
                                        {
                                            MzSpectrum complementarySpectrum = GenerateComplementarySpectrum(scan.TheScan.MassSpectrum, scan.PrecursorMass, commonParameters.DissociationType);
                                            matchedIons.AddRange(MatchFragmentIons(complementarySpectrum, peptideTheorProducts, commonParameters));
                                        }

                                        updatedPSM.AddOrReplace(currentPwsmWithNotch.pwsm, PSM.Score, currentPwsmWithNotch.notch, commonParameters.ReportAllAmbiguity, matchedIons);
                                    }
                                    PeptideSpectralMatches[i] = updatedPSM;
                                    break;
                                }
                            }
                        }
                    }
                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing nonspecific search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        //public static IEnumerable<Product> Fragment( DissociationType dissociationType, FragmentationTerminus fragmentationTerminus, PeptideWithSetModifications pwsm)
        //{
        //    // molecular ion
        //    //yield return new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.None, this.MonoisotopicMass, Length, Length), 0);

        //    var productCollection = TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[fragmentationTerminus].Intersect(DissociationTypeCollection.ProductsFromDissociationType[dissociationType]);

        //    List<(ProductType, int)> skippers = new List<(ProductType, int)>();
        //    foreach (var product in productCollection.Where(f => f != ProductType.zPlusOne))
        //    {
        //        skippers.Add((product, pwsm.BaseSequence.Length));
        //    }

        //    switch (dissociationType)
        //    {
        //        case DissociationType.CID:
        //            skippers.Add((ProductType.b, 1));
        //            break;

        //        case DissociationType.ETD:
        //        case DissociationType.ECD:
        //            break;
        //    }
        //         Dictionary<FragmentationTerminus, CompactPeptide> _compactPeptides = new Dictionary<FragmentationTerminus, CompactPeptide>();

        //    foreach (var productType in productCollection)
        //    {
        //        // we're separating the N and C terminal masses and computing a separate compact peptide for each one
        //        // this speeds calculations up without producing unnecessary terminus fragment info
        //        FragmentationTerminus temporaryFragmentationTerminus = TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[productType];

        //NeutralTerminusFragment[] terminalMasses = CompactPeptide(_compactPeptides, temporaryFragmentationTerminus, pwsm).TerminalMasses;

        //        for (int f = 0; f < terminalMasses.Length; f++)
        //        {
        //            // fragments with neutral loss
        //            if (pwsm.AllModsOneIsNterminus.TryGetValue(terminalMasses[f].AminoAcidPosition + 1, out Modification mod) && mod.NeutralLosses != null
        //                && mod.NeutralLosses.TryGetValue(dissociationType, out List<double> neutralLosses))
        //            {
        //                foreach (double neutralLoss in mod.NeutralLosses[dissociationType])
        //                {
        //                    if (neutralLoss == 0)
        //                    {
        //                        continue;
        //                    }

        //                    for (int n = f; n < terminalMasses.Length; n++)
        //                    {
        //                        if (!skippers.Contains((productType, terminalMasses[n].FragmentNumber)))
        //                        {
        //                            yield return new Product(productType, terminalMasses[n], neutralLoss);
        //                        }
        //                    }
        //                }
        //            }

        //            // "normal" fragment without neutral loss
        //            if (!skippers.Contains((productType, terminalMasses[f].FragmentNumber)))
        //            {
        //                yield return new Product(productType, terminalMasses[f], 0);
        //            }
        //        }
        //    }

        //    if (pwsm.AllModsOneIsNterminus != null)
        //    {
        //        foreach (Modification mod in pwsm.AllModsOneIsNterminus.Values)
        //        {
        //            // molecular ion minus neutral losses
        //            if (mod.NeutralLosses != null && mod.NeutralLosses.TryGetValue(dissociationType, out List<double> losses))
        //            {
        //                foreach (double neutralLoss in losses)
        //                {
        //                    if (neutralLoss != 0)
        //                    {
        //                        yield return new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, pwsm.MonoisotopicMass, 0, 0), neutralLoss);
        //                    }
        //                }
        //            }

        //            // diagnostic ions
        //            if (mod.DiagnosticIons != null && mod.DiagnosticIons.TryGetValue(dissociationType, out List<double> diagnosticIons))
        //            {
        //                foreach (double diagnosticIon in diagnosticIons)
        //                {
        //                    // the diagnostic ion is assumed to be annotated in the mod info as the *neutral mass* of the diagnostic ion, not the ionized species
        //                    yield return new Product(ProductType.D, new NeutralTerminusFragment(FragmentationTerminus.Both, diagnosticIon, 0, 0), 0);
        //                }
        //            }
        //        }
        //    }
        //}


        //public static CompactPeptide CompactPeptide(Dictionary<FragmentationTerminus, CompactPeptide> _compactPeptides, FragmentationTerminus fragmentationTerminus, PeptideWithSetModifications pwsm)
        //{
        //    // need this for deserialization
        //    if (_compactPeptides == null)
        //    {
        //        _compactPeptides = new Dictionary<FragmentationTerminus, CompactPeptide>();
        //    }

        //    if (_compactPeptides.TryGetValue(fragmentationTerminus, out CompactPeptide compactPeptide))
        //    {
        //        return compactPeptide;
        //    }
        //    else
        //    {
        //        CompactPeptide cp = new CompactPeptide(pwsm, fragmentationTerminus);
        //        _compactPeptides.Add(fragmentationTerminus, cp);
        //        return cp;
        //    }
        //}


        private Tuple<int, double> Accepts(double scanPrecursorMass, PeptideWithSetModifications peptide, FragmentationTerminus fragmentationTerminus, MassDiffAcceptor searchMode)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = commonParameters.DigestionParams.MinPeptideLength;
            List<Product> fragments = peptide.Fragment(commonParameters.DissociationType, fragmentationTerminus).ToList();//.Select(m => m.NeutralMass).ToArray();

            for (int i = localminPeptideLength; i < fragments.Count(); i++)
            {
                Product fragment = fragments[i];
                double theoMass = fragment.NeutralMass - DissociationTypeCollection.GetMassShiftFromProductType(fragment.ProductType) + WaterMonoisotopicMass;
                int notch = searchMode.Accepts(scanPrecursorMass, theoMass);
                if (notch >= 0)
                {
                    return new Tuple<int, double>(notch, theoMass);
                }
                else if (theoMass > scanPrecursorMass)
                {
                    break;
                }
            }
            //if the theoretical and experimental have the same mass
            if (fragments.Count > localminPeptideLength)
            {
                double totalMass = peptide.MonoisotopicMass;// + Constants.ProtonMass;
                int notch = searchMode.Accepts(scanPrecursorMass, totalMass);
                if (notch >= 0)
                {
                    return new Tuple<int, double>(notch, totalMass);
                }
            }
            return new Tuple<int, double>(-1, -1);
        }

        private int ComputePeptideIndexes(PeptideWithSetModifications yyy, ref double prevMass, int oneBasedIndexToLookAt, int direction, double precursorMass, MassDiffAcceptor massDiffAcceptor)
        {
            Modification residue_variable_mod = null;
            do
            {
                prevMass += Residue.ResidueMonoisotopicMass[yyy[oneBasedIndexToLookAt - 1]];

                yyy.AllModsOneIsNterminus.TryGetValue(oneBasedIndexToLookAt + 1, out residue_variable_mod);
                if (residue_variable_mod == null)
                {
                    if (massDiffAcceptor.Accepts(precursorMass, prevMass) >= 0)
                    {
                        return oneBasedIndexToLookAt;
                    }
                }
                else if (residue_variable_mod.NeutralLosses == null || residue_variable_mod.NeutralLosses.Count == 0)
                {
                    // TODO: missing case: mod is present and no neutral loss??
                }
                else if (residue_variable_mod.NeutralLosses.Count == 1)
                {
                    double monoisotopic = residue_variable_mod.MonoisotopicMass ?? 0;
                    prevMass += monoisotopic - residue_variable_mod.NeutralLosses[commonParameters.DissociationType].First();
                    if (massDiffAcceptor.Accepts(precursorMass, prevMass) >= 0)
                    {
                        return oneBasedIndexToLookAt;
                    }
                }
                else
                {
                    foreach (double nl in residue_variable_mod.NeutralLosses[commonParameters.DissociationType])
                    {
                        double monoisotopic = residue_variable_mod.MonoisotopicMass ?? 0;
                        prevMass += monoisotopic - nl;
                        if (massDiffAcceptor.Accepts(precursorMass, prevMass) >= 0)
                        {
                            return oneBasedIndexToLookAt;
                        }
                        if ((direction == 1 && oneBasedIndexToLookAt + direction < yyy.Length) ||
                            (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                        {
                            return ComputePeptideIndexes(yyy, ref prevMass, oneBasedIndexToLookAt + direction, direction, precursorMass, massDiffAcceptor);
                        }
                    }
                    break;
                }
                oneBasedIndexToLookAt += direction;
            } while ((oneBasedIndexToLookAt >= 1 && direction == -1) || (oneBasedIndexToLookAt <= yyy.Length && direction == 1));
            return -1;
        }
    }
}