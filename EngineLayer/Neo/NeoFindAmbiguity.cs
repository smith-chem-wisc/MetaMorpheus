using Chemistry;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Neo
{
    public static class NeoFindAmbiguity
    {
        //private BackgroundWorker worker = null;
        public static List<Protein> TheoreticalProteins;

        public static HashSet<string> NotFoundSequences = new HashSet<string>();
        public static Dictionary<string, List<Protein>> FoundSequences = new Dictionary<string, List<Protein>>();
        public static Dictionary<double, string[]> MassDict = new Dictionary<double, string[]>();
        public static double[] Keys;
        public static double ProductMassTolerancePpm = 20;

        //(Ppm)
        public static double PrecursorMassTolerancePpm = 5;

        //(Ppm)
        public static List<ProductType> IonsUsed = new List<ProductType> { ProductType.B, ProductType.Y };

        public static Dictionary<double, char> MassesToResidues = new Dictionary<double, char>();
        public static List<double> SingleAminoAcidMasses = new List<double>();
        public static double MaxDifference;
        public static Dictionary<string, List<string>> NTermDictionary = new Dictionary<string, List<string>>();
        public static Dictionary<string, List<string>> CTermDictionary = new Dictionary<string, List<string>>();
        public static Dictionary<string, List<Protein>> ProtDictionary = new Dictionary<string, List<Protein>>();
        public static char[] AANames = new char[20] { 'G', 'A', 'S', 'P', 'V', 'T', 'L', 'I', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'C', 'Y', 'W' };

        private const int MaxMissingConsecutivePeaks = 2;
        private const int MaxNumPossibleSequences = 2000;
        private const int DecimalDigitsForFragmentMassRounding = 3;
        private static readonly double WaterMonoisotopicMass = Math.Round(PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass, DecimalDigitsForFragmentMassRounding);
        //20 common AA, ordered by mass assuming carbamido

        public static void FindAmbiguity(List<NeoPsm> candidates, List<Protein> theoreticalProteins, Ms2ScanWithSpecificMass[] spectra, string databaseFileName)
        {
            PopulateSequenceLookUpDictionaries(databaseFileName, theoreticalProteins);
            Ms2ScanWithSpecificMass[] indexedSpectra = new Ms2ScanWithSpecificMass[spectra.Max(x => x.OneBasedScanNumber) + 1];
            foreach (Ms2ScanWithSpecificMass scan in spectra)
            {
                indexedSpectra[scan.OneBasedScanNumber] = scan;
            }
            for (int i = 0; i < candidates.Count; i++) //must be mutable while iterating
            {
                NeoPsm psm = candidates[i];
                Ms2ScanWithSpecificMass spectrum = indexedSpectra[psm.ScanNumber];
                psm.FusionType = FusionType.TS; //for some maddening reason, this is not arriving here as trans, but instead translated
                if (IsTooMessy(psm, spectrum)) //having explosion of combinations when greater than 3 consequtive peaks producing tens of thousands of sequences ids, causes hanging
                {
                    candidates.RemoveAt(i);
                    i--;
                }
                else
                {
                    if (!GeneratePossibleSequences(psm, spectrum)) //return true if fewer than specified number of ambiguities
                    {
                        candidates.RemoveAt(i);
                        i--;
                    }
                }
            }
            NTermDictionary.Clear();
            CTermDictionary.Clear();
            ProtDictionary.Clear();
        }

        //use ion hits to know where peaks have been found by morpheus and where there is ambiguity
        public static void FindIons(FusionCandidate fusionCandidate, NeoPsm psm, Ms2ScanWithSpecificMass spectrum)
        {
            fusionCandidate.makeFoundIons();
            string candSeq = fusionCandidate.Seq;
            bool[] foundIons = fusionCandidate.FoundIons;

            //find which aa have peaks
            for (int i = 0; i < foundIons.Length - 1; i++)
            {
                //B IONS//
                if (IonsUsed.Contains(ProductType.B))
                {
                    double bTheoMass = ClassExtensions.ToMz(NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i)) - NeoConstants.WATER_MONOISOTOPIC_MASS, 1);
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, bTheoMass, ProductMassTolerancePpm))
                        {
                            foundIons[i + 1] = true;
                        }
                    }
                }
                //Y IONS//
                if (IonsUsed.Contains(ProductType.Y))
                {
                    double yTheoMass = ClassExtensions.ToMz(NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1)), 1);
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, yTheoMass, ProductMassTolerancePpm))
                        {
                            foundIons[foundIons.Length - 1 - i] = true;
                        }
                    }
                }
                //C IONS//
                if (IonsUsed.Contains(ProductType.C))
                {
                    double cTheoMass = ClassExtensions.ToMz(NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i)) - NeoConstants.WATER_MONOISOTOPIC_MASS + NeoConstants.nitrogenMonoisotopicMass + 3 * NeoConstants.hydrogenMonoisotopicMass, 1);
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, cTheoMass, ProductMassTolerancePpm))
                        {
                            foundIons[i + 1] = true;
                        }
                    }
                }
                //ZDOT IONS//
                if (IonsUsed.Contains(ProductType.Zdot))
                {
                    double zdotTheoMass = ClassExtensions.ToMz(NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1)) - NeoConstants.nitrogenMonoisotopicMass - 2 * NeoConstants.hydrogenMonoisotopicMass, 1);
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, zdotTheoMass, ProductMassTolerancePpm))
                        {
                            foundIons[foundIons.Count() - 1 - i] = true;
                        }
                    }
                }
            }
            foundIons[0] = true;//|A|B|C|D|E|F|K where the whole peptide peak is always placed arbitrarily at the n term
        }

        public static void ReadMassDictionary()
        {
            List<double> tempKeys = new List<double>();

            using (StreamReader masses = new StreamReader(Path.Combine(GlobalVariables.DataDir, @"Neo\\Data\Dictionary" + MaxMissingConsecutivePeaks + ".txt"))) //file located in Morpheus folder
            {
                while (masses.Peek() != -1)
                {
                    string line = masses.ReadLine();
                    string[] fields = line.Split('\t');
                    double key = Convert.ToDouble(fields[0]);
                    string[] sequences = fields[1].Split(';');
                    MassDict.Add(key, sequences);
                    tempKeys.Add(key);
                }
            }
            Keys = new double[tempKeys.Count];
            for (int i = 0; i < tempKeys.Count; i++)
            {
                Keys[i] = tempKeys[i];
            }
        }

        private static bool IsTooMessy(NeoPsm psm, Ms2ScanWithSpecificMass spectrum) //return true if too messy for confident identification
        {
            List<string> baseSequences = new List<string>();
            int currentBestScore = 0;
            for (int index = 0; index < psm.Candidates.Count; index++)
            {
                if (psm.Candidates[index].Seq.Length < 5)
                {
                    psm.Candidates.RemoveAt(index);
                    index--;
                    continue;
                }
                bool badID = false;
                FusionCandidate fc = psm.Candidates[index];
                FindIons(fc, psm, spectrum);
                int consecutiveMissedCounter = 0;
                int totalHitCounter = 0;
                if (fc.FoundIons.Count(x => x) * 2 < fc.FoundIons.Length)
                {
                    psm.Candidates.RemoveAt(index);
                    index--;
                    continue;
                }
                for (int b = 0; b < fc.FoundIons.Length; b++)
                {
                    if (consecutiveMissedCounter > MaxMissingConsecutivePeaks) //if too many permutations possible because of an unmapped region
                        badID = true;
                    else if (!fc.FoundIons[b])
                        consecutiveMissedCounter++;
                    else
                    {
                        totalHitCounter++;
                        consecutiveMissedCounter = 0;
                    } //only care about consecutive
                }

                if (totalHitCounter > currentBestScore && !badID)//the others were worse, so delete them
                {
                    for (int i = 0; i < index; i = 0)
                    {
                        psm.Candidates.Remove(psm.Candidates[0]);
                        index--;
                    }
                    currentBestScore = totalHitCounter;
                    baseSequences = new List<string> { psm.Candidates[index].Seq };
                }
                else if (totalHitCounter < currentBestScore || badID || baseSequences.Contains(psm.Candidates[index].Seq))
                {
                    psm.Candidates.Remove(psm.Candidates[index]);
                    index--;
                }
            }
            //If there's anything left
            if (psm.Candidates.Count > 0) //It wasn't too messy! Yay!
            {
                List<FusionCandidate> nTermCands = new List<FusionCandidate>();
                foreach (FusionCandidate oldCand in psm.Candidates)
                {
                    oldCand.FoundIons[1] = false;
                    oldCand.FoundIons[2] = false;
                    if (!oldCand.FoundIons[3])
                    {
                        FusionCandidate fc = new FusionCandidate(oldCand.Seq);
                        FindIons(fc, psm, spectrum);
                        fc.FoundIons[1] = true;
                        fc.FoundIons[2] = false;
                        if (!oldCand.FoundIons[4])
                        {
                            fc.FoundIons[4] = true;
                        }
                        nTermCands.Add(fc);
                        fc = new FusionCandidate(oldCand.Seq);
                        FindIons(fc, psm, spectrum);
                        fc.FoundIons[1] = false;
                        fc.FoundIons[2] = true;
                        nTermCands.Add(fc);
                        fc = new FusionCandidate(oldCand.Seq);
                        FindIons(fc, psm, spectrum);
                        fc.FoundIons[1] = false;
                        fc.FoundIons[2] = false;
                        fc.FoundIons[3] = true;
                        nTermCands.Add(fc);
                    }
                }
                foreach (FusionCandidate cand in nTermCands)
                    psm.Candidates.Add(cand);
                return false;
            }
            else //this might be a fusion peptide, but we won't get enough valuable information from this spectra, so discard it
                return true;
        }

        private static void ReadAminoAcids()
        {
            for (int i = 0; i < Residue.ResidueMonoisotopicMass.Length; i++)
            {
                if (!double.IsNaN(Residue.ResidueMonoisotopicMass[i]) && !SingleAminoAcidMasses.Contains(Residue.ResidueMonoisotopicMass[i]))
                {
                    MassesToResidues.Add(Residue.ResidueMonoisotopicMass[i], (char)i);
                    SingleAminoAcidMasses.Add(Residue.ResidueMonoisotopicMass[i]);
                }
            }
            SingleAminoAcidMasses.Sort();
            MaxDifference = SingleAminoAcidMasses[SingleAminoAcidMasses.Count - 1];
        }

        private static bool GeneratePossibleSequences(NeoPsm psm, Ms2ScanWithSpecificMass spectrum) //returns false if over the specified number of sequences are generated
        {
            List<string> ambiguousCandidates = new List<string>();
            foreach (FusionCandidate fusionCandidate in psm.Candidates)
            {
                List<string> nMatchedSequences = new List<string>();
                List<string> cMatchedSequences = new List<string>();
                List<string> partialSequences = new List<string> { "" };
                double theoreticalMass = NeoMassCalculator.MonoIsoptopicMass(fusionCandidate.Seq);
                //FindIons(fusionCandidate, psm, spectrum, out string error_message1); //populate the foundIons array
                //error_message += error_message1;

                int globalIndex = 0;
                //C-Terminus
                //Increment through the amino acids to determine all possible sequences
                //once a sequence is no longer findable (is not in database), then stop generating possible sequences for that sequence
                //stop when reaching end of peptide or no more present in database
                while (fusionCandidate.FoundIons.Length > globalIndex)
                {
                    int reverseIndex = fusionCandidate.FoundIons.Length - globalIndex - 1;
                    int mostRecent = fusionCandidate.FoundIons.Length; //most recent Ion found prior to this one (start point)
                    if (fusionCandidate.FoundIons[reverseIndex])
                    {
                        for (int i = reverseIndex + 1; i < fusionCandidate.FoundIons.Length; i++) //identify start point
                        {
                            if (fusionCandidate.FoundIons[i])
                            {
                                mostRecent = i; //save most recent hit, exclusive of the current index
                                break;
                            }
                        }

                        //get combos
                        string ambiguousFrag = fusionCandidate.Seq.Substring(reverseIndex, mostRecent - reverseIndex);
                        double key = NeoMassCalculator.MonoIsoptopicMass(ambiguousFrag);

                        List<string> combinations = GetCombinations(key, theoreticalMass - key, psm.ExpMass);

                        List<string> tempPartialSequences = new List<string>();
                        foreach (string old in partialSequences)
                            foreach (string neu in combinations)
                                tempPartialSequences.Add(neu + old);

                        partialSequences.Clear();

                        if (tempPartialSequences.Count == 0)
                        {
                            break;
                        }

                        for (int i = tempPartialSequences.Count - 1; i >= 0; i--)
                        {
                            string partialSeq = tempPartialSequences[i];
                            if (partialSeq.Length > 3)
                            {
                                //can it be made?
                                if (CTermDictionary.TryGetValue(partialSeq.Substring(partialSeq.Length - 4, 4), out List<string> cEntry))
                                {
                                    int n = 4;
                                    while (cEntry.Count != 0 && n <= partialSeq.Length)
                                    {
                                        string currentSeq = partialSeq.Substring(partialSeq.Length - n, n);
                                        if (!cMatchedSequences.Any(x => x.Length >= n && x.Substring(x.Length - n, n).Equals(currentSeq)))
                                        {
                                            cEntry = cEntry.AsParallel().Where(seq => seq.Length >= n && seq.Substring(seq.Length - n, n).Equals(currentSeq)).ToList();
                                        }
                                        n++;
                                    }
                                    if (cEntry.Count == 0)
                                    {
                                        string strToAdd = partialSeq.Substring(partialSeq.Length - n + 2, n - 2);
                                        if (!cMatchedSequences.Any(x => x.Equals(strToAdd)))
                                        {
                                            for (int j = cMatchedSequences.Count - 1; j >= 0; j--)
                                            {
                                                if (strToAdd.Contains(cMatchedSequences[j]))
                                                {
                                                    cMatchedSequences.RemoveAt(j);
                                                }
                                            }
                                            if (!cMatchedSequences.Any(x => x.Contains(strToAdd)))
                                            {
                                                cMatchedSequences.Add(strToAdd);
                                            }
                                        }
                                    }
                                    else
                                    {
                                        partialSequences.Add(partialSeq); //keep adding stuff to it
                                    }
                                }
                                else
                                {
                                    // do nothing
                                }
                            }
                            else
                            {
                                partialSequences.Add(partialSeq); //keep adding stuff to it
                            }
                        }
                    }
                    globalIndex++;
                }

                foreach (string seq in partialSequences)
                    cMatchedSequences.Add(seq);

                partialSequences.Clear();
                partialSequences.Add("");
                globalIndex = 1;

                //N-terminus
                while (fusionCandidate.FoundIons.Length >= globalIndex)
                {
                    int mostRecent = 0; //most recent Ion found prior to this one (start point)
                    if (fusionCandidate.FoundIons.Length == globalIndex || fusionCandidate.FoundIons[globalIndex])
                    {
                        for (int i = globalIndex - 1; i >= 0; i--) //identify start point
                            if (fusionCandidate.FoundIons[i])
                            {
                                mostRecent = i; //save most recent hit, exclusive of the current index
                                break;
                            }

                        //get combos
                        string ambiguousFrag = fusionCandidate.Seq.Substring(mostRecent, globalIndex - mostRecent);
                        double key = NeoMassCalculator.MonoIsoptopicMass(ambiguousFrag);

                        List<string> combinations = GetCombinations(key, theoreticalMass - key, psm.ExpMass);

                        List<string> tempPartialSequences = new List<string>();
                        foreach (string old in partialSequences)
                            foreach (string neu in combinations)
                                tempPartialSequences.Add(old + neu);

                        partialSequences.Clear();

                        for (int i = tempPartialSequences.Count - 1; i >= 0; i--)
                        {
                            string partialSeq = tempPartialSequences[i];

                            if (partialSeq.Length > 3)
                            {
                                //can it be made?
                                int n = 4;

                                if (NTermDictionary.TryGetValue(partialSeq.Substring(0, 4), out List<string> nEntry))
                                {
                                    while (nEntry.Count != 0 && n <= partialSeq.Length)
                                    {
                                        string currentSeq = partialSeq.Substring(0, n);
                                        if (!nMatchedSequences.Any(x => x.Length >= n && x.Substring(0, n).Equals(currentSeq)))
                                            nEntry = nEntry.AsParallel().Where(seq => seq.Length >= n && seq.Substring(0, n).Equals(currentSeq)).ToList();
                                        n++;
                                    }
                                    if (nEntry.Count == 0)
                                    {
                                        string strToAdd = partialSeq.Substring(0, n - 2);
                                        if (!nMatchedSequences.Any(x => x.Equals(strToAdd)))
                                        {
                                            for (int j = nMatchedSequences.Count - 1; j >= 0; j--)
                                                if (strToAdd.Contains(nMatchedSequences[j]))
                                                    nMatchedSequences.RemoveAt(j);
                                            if (!nMatchedSequences.AsParallel().Any(x => x.Contains(strToAdd)))
                                                nMatchedSequences.Add(strToAdd);
                                        }
                                    }
                                    else
                                        partialSequences.Add(partialSeq);
                                }
                                else { }
                            }
                            else
                                partialSequences.Add(partialSeq);
                        }
                    }
                    globalIndex++;
                }
                foreach (string seq in partialSequences)
                    nMatchedSequences.Add(seq);

                //Splice surviving fragments
                Parallel.ForEach(nMatchedSequences, n =>
                {
                    double nMass = NeoMassCalculator.MonoIsoptopicMass(n);
                    foreach (string c in cMatchedSequences)
                    {
                        lock (ambiguousCandidates)
                        {
                            if (ambiguousCandidates.Count > MaxNumPossibleSequences)
                                break;
                        }
                        double cMass = NeoMassCalculator.MonoIsoptopicMass(c);
                        double totalMass = nMass + cMass - NeoConstants.WATER_MONOISOTOPIC_MASS;
                        if (totalMass + 1 > psm.ExpMass)
                        {
                            if (NeoMassCalculator.IdenticalMasses(psm.ExpMass, totalMass, PrecursorMassTolerancePpm))
                            {
                                lock (ambiguousCandidates)
                                {
                                    if (!ambiguousCandidates.Contains(n + c))
                                        ambiguousCandidates.Add(n + c);
                                }
                            }
                            else
                            {
                                double originalTotalMass = totalMass;
                                totalMass -= NeoMassCalculator.GetMonoisotopicMass(n[n.Length - 1], n);
                                string nFrag = n.Substring(0, n.Length - 1);
                                string cFrag = c;
                                while (true)
                                {
                                    if (NeoMassCalculator.IdenticalMasses(psm.ExpMass, totalMass, PrecursorMassTolerancePpm))
                                    {
                                        lock (ambiguousCandidates)
                                        {
                                            if (!ambiguousCandidates.Contains(nFrag + cFrag))
                                                ambiguousCandidates.Add(nFrag + cFrag);
                                        }
                                        if (cFrag.Length == 0)
                                            break;
                                        nFrag = n;
                                        originalTotalMass -= NeoMassCalculator.GetMonoisotopicMass(cFrag[0], cFrag);
                                        totalMass = originalTotalMass;
                                        cFrag = cFrag.Substring(1, cFrag.Length - 1);
                                    }
                                    else if (totalMass + 10 < psm.ExpMass)
                                    {
                                        if (nFrag.Equals(n))
                                            break;
                                        else
                                        {
                                            nFrag = n;
                                            if (cFrag.Length == 0)
                                                break;
                                            originalTotalMass -= NeoMassCalculator.GetMonoisotopicMass(cFrag[0], cFrag);
                                            totalMass = originalTotalMass;
                                            cFrag = cFrag.Substring(1, cFrag.Length - 1);
                                        }
                                    }
                                    else
                                    {
                                        if (nFrag.Length == 0)
                                            break;
                                        totalMass -= NeoMassCalculator.GetMonoisotopicMass(nFrag[nFrag.Length - 1], n);
                                        nFrag = nFrag.Substring(0, nFrag.Length - 1);
                                    }
                                }
                            }
                        }
                    }
                });
            }

            psm.Candidates.Clear();
            if (ambiguousCandidates.Count <= MaxNumPossibleSequences)
            {
                //CIS CODE
                foreach (string fc in ambiguousCandidates)
                {
                    if (fc.Length < 5)
                        continue;

                    string substring = fc.Substring(0, 4);
                    NTermDictionary.TryGetValue(fc.Substring(0, 4), out List<string> nFrag);
                    if (nFrag != null && nFrag.AsParallel().Any(seq => seq.Length >= fc.Length && seq.Substring(0, fc.Length).Equals(fc)))
                    {
                        psm.FusionType = FusionType.TL;
                        psm.Candidates.Add(new FusionCandidate(fc) { FusionType = FusionType.TL });
                    }
                    else
                    {
                        //Is it cis?
                        bool cis = false;
                        List<Protein> otherPossibleProteins;
                        if (ProtDictionary.TryGetValue(substring, out List<Protein> possibleProteins))
                        {
                            //possibleProteins.ForEach(prot => originalPossibleProteins.Add(prot));
                            for (int i = 5; i <= fc.Length; i++)
                            {
                                string otherSubstring = fc.Substring(i - 1, fc.Length - i + 1);
                                //get proteins containing both halves
                                otherPossibleProteins = possibleProteins.AsParallel().Where(prot => prot.BaseSequence.Contains(otherSubstring)).ToList();

                                //check if both halves are the correct distance apart and are not overlapping
                                foreach (Protein prot in otherPossibleProteins)
                                {
                                    string seq = prot.BaseSequence;
                                    List<int> indexes = new List<int>();
                                    List<int> otherIndexes = new List<int>();
                                    int index = seq.IndexOf(substring);
                                    int maxStartingIndex = Math.Max(0, index - 25 - fc.Length);
                                    int otherIndex = seq.IndexOf(otherSubstring, maxStartingIndex);

                                    while (index != -1)
                                    {
                                        indexes.Add(index);
                                        index = seq.IndexOf(substring, index + 1);
                                    }
                                    while (otherIndex != -1)
                                    {
                                        otherIndexes.Add(otherIndex);
                                        otherIndex = seq.IndexOf(otherSubstring, otherIndex + 1);
                                    }

                                    int n = 0;
                                    int c = 0;
                                    while (n < indexes.Count && c < otherIndexes.Count)
                                    {
                                        int difference = otherIndexes[c] - (indexes[n] + i);
                                        if (difference > 25)
                                            n++;
                                        else if (difference < -25 - fc.Length)
                                            c++;
                                        else if (difference <= 0 && difference > -fc.Length)
                                        {
                                            int originalC = c;
                                            c++;
                                            while (n < indexes.Count && c < otherIndexes.Count)
                                            {
                                                difference = otherIndexes[c] - (indexes[n] + i);
                                                if (difference > 25)
                                                {
                                                    c = originalC;
                                                    n++;
                                                }
                                                else if (difference <= 0 && difference > -fc.Length)
                                                    c++;
                                                else //it's cis
                                                    break;
                                            }
                                        }
                                        else
                                        {
                                            //it's cis!
                                            psm.Candidates.Add(new FusionCandidate(fc) { FusionType = FusionType.NC });
                                            cis = true;
                                            if (psm.FusionType == FusionType.TS)
                                                psm.FusionType = FusionType.NC;
                                        }
                                        if (cis)
                                            break;
                                    }
                                    if (cis)
                                        break;
                                }
                                if (cis)
                                    break;

                                substring = fc.Substring(0, i);
                                possibleProteins = possibleProteins.AsParallel().Where(prot => prot.BaseSequence.Contains(substring)).ToList();
                                if (possibleProteins.Count == 0)
                                    break;
                            }
                        }
                        if (!cis)
                        {
                            if (ProtDictionary.TryGetValue(fc.Substring(fc.Length - 4, 4), out otherPossibleProteins))
                            {
                                string otherSubstring = fc.Substring(3, fc.Length - 3);
                                otherPossibleProteins = otherPossibleProteins.Where(prot => prot.BaseSequence.Contains(otherSubstring)).ToList();
                                for (int i = 2; i >= 0; i--)
                                {
                                    substring = fc.Substring(0, i + 1);
                                    foreach (Protein prot in otherPossibleProteins)
                                    {
                                        string seq = prot.BaseSequence;
                                        List<int> indexes = new List<int>();
                                        List<int> otherIndexes = new List<int>();
                                        int index = seq.IndexOf(substring);
                                        int maxStartingIndex = Math.Max(0, index - 25 - fc.Length);
                                        int otherIndex = seq.IndexOf(otherSubstring, maxStartingIndex);

                                        while (index != -1)
                                        {
                                            indexes.Add(index);
                                            index = seq.IndexOf(substring, index + 1);
                                        }
                                        while (otherIndex != -1)
                                        {
                                            otherIndexes.Add(otherIndex);
                                            otherIndex = seq.IndexOf(otherSubstring, otherIndex + 1);
                                        }

                                        int n = 0;
                                        int c = 0;
                                        while (n < indexes.Count && c < otherIndexes.Count)
                                        {
                                            int difference = otherIndexes[c] - (indexes[n] + i);
                                            if (difference > 25)
                                                n++;
                                            else if (difference < -25 - fc.Length)
                                                c++;
                                            else if (difference <= 0 && difference > -fc.Length)
                                            {
                                                int originalC = c;
                                                c++;
                                                while (n < indexes.Count && c < otherIndexes.Count)
                                                {
                                                    difference = otherIndexes[c] - (indexes[n] + i);
                                                    if (difference > 25)
                                                    {
                                                        c = originalC;
                                                        n++;
                                                    }
                                                    else if (difference <= 0 && difference > -fc.Length)
                                                        c++;
                                                    else //it's cis
                                                        break;
                                                }
                                            }
                                            else
                                            {
                                                //it's cis!
                                                psm.Candidates.Add(new FusionCandidate(fc) { FusionType = FusionType.NC });
                                                cis = true;
                                                if (psm.FusionType == FusionType.TS)
                                                    psm.FusionType = FusionType.NC;
                                            }
                                            if (cis)
                                                break;
                                        }
                                        if (cis)
                                            break;
                                    }
                                    if (cis)
                                        break;

                                    otherSubstring = fc.Substring(i, fc.Length - i);
                                }
                            }
                            if (!cis)
                                psm.Candidates.Add(new FusionCandidate(fc));
                        }
                    }
                }
            }
            return psm.Candidates.Count > 0;
        }

        //returns false if a full fusion sequence could not be made or was found in the database, making it translated instead of a novel fusion.
        private static bool PossibleCandidate(NeoPsm psm)
        {
            FoundSequences = new Dictionary<string, List<Protein>>(); //used for finding longer fragments than those previously identified. Also populates ParentInfo
            NotFoundSequences = new HashSet<string>(); //don't bother looking for these fragments, since we know they don't exist. Good for multiple homologous putative fusion peptide sequences

            //conduct an initial search of each candidate's full sequence to identify any that are translated
            for (int i = 0; i < psm.Candidates.Count; i++) //foreach fusion peptide sequence that could map to this scan
            {
                string novelSeq = psm.Candidates[i].Seq;
                if (foundParent(novelSeq, ParentInfo.Terminal.C, psm.Candidates[i], false)) //check really quick to see if the whole thing exists as is. If so, assign it as translated. Terminal C was arbitrarily chosen
                {
                    foreach (ParentInfo info in psm.Candidates[i].ParentInfo)
                        foreach (Protein protein in info.TheoreticalProteins)
                            if (protein.BaseSequence.Contains(novelSeq)) //if translated
                                psm.Candidates[i].TranslatedParents.Add(new TranslatedParent(protein.Accession, protein.BaseSequence, protein.BaseSequence.IndexOf(psm.Candidates[i].Seq), psm.Candidates[i].Seq.Length));

                    psm.Candidates[i].FusionType = FusionType.TL;
                    psm.FusionType = psm.Candidates[i].FusionType;
                    for (int j = 0; j < psm.Candidates.Count; j++)
                    {
                        if (j != i)
                        {
                            psm.Candidates.Remove(psm.Candidates[j]);
                            j--;
                            i--;
                        }
                    }
                    return false;
                }
            }
            for (int i = 0; i < psm.Candidates.Count; i++) //foreach fusion peptide sequence that could map to this scan
            {
                //sw.StartFindParents
                if (!isViable(psm.Candidates[i])) //remove this fusion peptide sequence if the parent fragments cannot be found with the given database
                {
                    psm.Candidates.Remove(psm.Candidates[i]);
                    i--;
                }
                else
                {
                    DetermineFusionCandidateType(psm.Candidates[i]); //cis, trans?
                    if (psm.FusionType > psm.Candidates[i].FusionType) //if more likely than previous types, change the psm type (golf scoring)
                    {
                        psm.FusionType = psm.Candidates[i].FusionType;
                    }
                    if (psm.FusionType.Equals(FusionType.TL))//if there's a possible sequence that's present in the database, it is likely correct and is it is not worth it to identify parents of other sequences.
                    {
                        //remove all other candidates
                        for (int j = 0; j < psm.Candidates.Count; j++)
                        {
                            if (j != i)
                            {
                                psm.Candidates.Remove(psm.Candidates[j]);
                                j--;
                                i--;
                            }
                        }
                        return false;
                    }
                }
            }

            return psm.Candidates.Count != 0; //if no candidates are left, we couldn't make the sequence with the database and we'll discard it.
        }

        private static List<string> GetCombinations(double key, double theoreticalMass, double experimentalMass)
        {
            List<string> combinations = new List<string>();
            double closestPeak = double.NaN;
            var ipos = Array.BinarySearch(Keys, key);
            if (ipos < 0)
                ipos = ~ipos;

            if (ipos > 0)
            {
                var downIpos = ipos - 1;
                // Try down
                while (downIpos >= 0)
                {
                    closestPeak = Keys[downIpos];
                    if (NeoMassCalculator.IdenticalMasses(experimentalMass, theoreticalMass + closestPeak, PrecursorMassTolerancePpm))
                    {
                        if (MassDict.TryGetValue(closestPeak, out string[] value))
                            foreach (string frag in value)
                                combinations.Add(frag);
                    }
                    else
                        break;
                    downIpos--;
                }
            }
            if (ipos < Keys.Length)
            {
                var upIpos = ipos;
                // Try here and up
                while (upIpos < Keys.Length)
                {
                    closestPeak = Keys[upIpos];
                    if (NeoMassCalculator.IdenticalMasses(experimentalMass, theoreticalMass + closestPeak, PrecursorMassTolerancePpm))
                    {
                        if (MassDict.TryGetValue(closestPeak, out string[] value))
                            foreach (string frag in value)
                                combinations.Add(frag);
                    }
                    else
                        break;
                    upIpos++;
                }
            }
            return combinations;
        }

        private static bool isViable(FusionCandidate tempCandidate) //returns if sequence could be made from one or two proteins in database and writes fusion type, parents, and junctions to fusionCandidate
        {
            //need to check that each index is viable
            string novelSeq = tempCandidate.Seq;
            //N//
            int nTermParentLength = novelSeq.Length - 1;//length-1, because if the whole thing existed we wouldn't have made it to the else loop. Several edits are made to reflect this
            if (nTermParentLength > 6) //used to speed up search by finding an ideal starting point
            {
                nTermParentLength = 6; //low point of random probability (5 or 7 may also be suitable)
            }

            bool done = false;
            bool foundFirstSearch = false;

            //First pass search
            string testFrag = novelSeq.Substring(0, nTermParentLength);
            if (foundParent(testFrag, ParentInfo.Terminal.N, tempCandidate, foundFirstSearch)) //if found
            {
                foundFirstSearch = true;
                nTermParentLength++;
            }
            else //if not found
            {
                foundFirstSearch = false;
                nTermParentLength--;
            }

            //All other passes
            while (nTermParentLength < novelSeq.Length && nTermParentLength > 0 && !done) //while in range and not done
            {
                testFrag = novelSeq.Substring(0, nTermParentLength);
                if (foundParent(testFrag, ParentInfo.Terminal.N, tempCandidate, foundFirstSearch)) //if found
                {
                    if (!foundFirstSearch)
                    {
                        nTermParentLength--;
                        done = true;
                    }
                    nTermParentLength++;
                }
                else //if not found
                {
                    if (foundFirstSearch)
                    {
                        done = true;
                    }
                    nTermParentLength--;
                }
            }

            //C//
            done = false; //reset tracker
            foundFirstSearch = false;
            int cTermParentLength = novelSeq.Length - 1;
            if (cTermParentLength > 6) //used to speed up search by finding an ideal starting point
            {
                cTermParentLength = 6; //low point of random probability
            }

            testFrag = novelSeq.Substring(novelSeq.Length - cTermParentLength, cTermParentLength);
            //First pass search
            if (foundParent(testFrag, ParentInfo.Terminal.C, tempCandidate, foundFirstSearch)) //if found
            {
                foundFirstSearch = true;
                cTermParentLength++;
            }
            else //if not found
            {
                foundFirstSearch = false;
                cTermParentLength--;
            }

            while (cTermParentLength > 0 && cTermParentLength < novelSeq.Length && !done)
            {
                testFrag = novelSeq.Substring(novelSeq.Length - cTermParentLength, cTermParentLength);
                if (foundParent(testFrag, ParentInfo.Terminal.C, tempCandidate, foundFirstSearch))
                {
                    if (!foundFirstSearch)
                    {
                        cTermParentLength--;
                        done = true;
                    }
                    cTermParentLength++;
                }
                else
                {
                    if (foundFirstSearch)
                    {
                        done = true;
                    }
                    cTermParentLength--;
                }
            }

            if (cTermParentLength + nTermParentLength < novelSeq.Length) //if no overlap
            {
                return false;
            }
            else
            {
                for (int junction = tempCandidate.Seq.Length - cTermParentLength - 1; junction < nTermParentLength; junction++)
                {
                    tempCandidate.addJunctionIndex(junction);
                }
                return true;
            }
        }

        private static bool foundParent(string frag, ParentInfo.Terminal terminal, FusionCandidate candidate, bool foundFirstSearch)
        {
            //localTheoreticals.AsParallel().Where(x => x.Contains(frag)).ToList();
            if (NotFoundSequences.Contains(frag)) //has the fragment been searched but not found before?
                return false;

            List<Protein> matches = new List<Protein>();
            if (FoundSequences.TryGetValue(frag, out matches)) //has the fragment been searched AND found before?
            {
                candidate.ParentInfo.Add(new ParentInfo(matches, terminal, frag));
                return true;
            }

            if (foundFirstSearch) //Has something smaller been found before? Well, then we can just search against those found sequences
            {
                string shorterFrag = terminal.Equals(ParentInfo.Terminal.N) ? frag.Substring(0, frag.Length - 1) : frag.Substring(1, frag.Length - 1);

                foreach (ParentInfo info in candidate.ParentInfo)
                {
                    if (info.ParentType.Equals(terminal) && info.FragFound.Equals(shorterFrag))
                    {
                        List<Protein> tempProtList = new List<Protein>();
                        info.TheoreticalProteins.ForEach(protein => tempProtList.Add(protein));
                        matches = tempProtList.AsParallel().Where(x => x.BaseSequence.Contains(frag)).ToList();
                    }
                }
            }
            else //it hasn't been found before... we need to search against the whole database :(
            {
                matches = TheoreticalProteins.AsParallel().Where(x => x.BaseSequence.Contains(frag)).ToList();
            }
            if (matches != null && matches.Count > 0)
            {
                FoundSequences.Add(frag, matches);
                candidate.ParentInfo.Add(new ParentInfo(matches, terminal, frag));
                return true;
            }
            else
            {
                NotFoundSequences.Add(frag);
                return false;
            }
        }

        private static void DetermineFusionCandidateType(FusionCandidate fusionCandidate)
        {
            if (!fusionCandidate.FusionType.Equals(FusionType.TL))
            {
                string sequence = fusionCandidate.Seq;
                foreach (ParentInfo info in fusionCandidate.ParentInfo)
                {
                    int foundLength = info.FragFound.Length;

                    string compFrag = "";// fusionCandidate.seq.Substring()
                    if (info.ParentType.Equals(ParentInfo.Terminal.N))
                    {
                        compFrag = fusionCandidate.Seq.Substring(foundLength, fusionCandidate.Seq.Length - foundLength);
                    }
                    else
                    {
                        compFrag = fusionCandidate.Seq.Substring(0, fusionCandidate.Seq.Length - foundLength);
                    }

                    foreach (Protein protein in info.TheoreticalProteins)
                    {
                        //get the index(es) of where the found fragment is
                        string subProt = protein.BaseSequence;
                        List<int> originalIndexes = new List<int>();
                        int pastIndex = 0;
                        while (subProt.Contains(info.FragFound))
                        {
                            int newIndex = subProt.IndexOf(info.FragFound);
                            originalIndexes.Add(newIndex + pastIndex);
                            subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                            pastIndex += newIndex + 1;
                        }
                        fusionCandidate.TransParents.Add(new TransParent(protein.Accession, protein.BaseSequence, originalIndexes, foundLength, info.ParentType));

                        //get the index(es) of where the complimentary fragment is (if it's a cis fusion peptide)
                        subProt = protein.BaseSequence;
                        List<int> complementaryIndexes = new List<int>();
                        pastIndex = 0;
                        while (subProt.Contains(compFrag))
                        {
                            int newIndex = subProt.IndexOf(compFrag);
                            complementaryIndexes.Add(newIndex + pastIndex);
                            subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                            pastIndex += newIndex + 1;
                        }
                        if (complementaryIndexes.Count > 0) //if it is not trans
                        {
                            //if it is cis
                            if (info.ParentType.Equals(ParentInfo.Terminal.N))
                                fusionCandidate.CisParents.Add(new CisParent(protein.Accession, protein.BaseSequence, originalIndexes, foundLength, complementaryIndexes, compFrag.Length));
                            else
                                fusionCandidate.CisParents.Add(new CisParent(protein.Accession, protein.BaseSequence, complementaryIndexes, compFrag.Length, originalIndexes, foundLength));
                        }
                    }
                }
            }
            else
            {
                string seq = fusionCandidate.Seq;
                foreach (ParentInfo info in fusionCandidate.ParentInfo)
                {
                    foreach (Protein protein in info.TheoreticalProteins)
                    {
                        if (protein.BaseSequence.Contains(seq)) //if translated
                        {
                            fusionCandidate.TranslatedParents.Add(new TranslatedParent(protein.Accession, protein.BaseSequence, protein.BaseSequence.IndexOf(fusionCandidate.Seq), fusionCandidate.Seq.Length));
                        }
                    }
                }
            }
            foreach (CisParent cisParent in fusionCandidate.CisParents)
            {
                if (cisParent.CisType < fusionCandidate.FusionType)
                {
                    fusionCandidate.FusionType = cisParent.CisType;
                }
            }
        }

        private static void RemoveTranslatablePeptides(List<NeoPsm> psms)
        {
            for (int i = 0; i < psms.Count; i++)
            {
                bool remove = false;
                foreach (FusionCandidate fc in psms[i].Candidates)
                {
                    if (fc.JunctionIndexes.Count > fc.Seq.Length) //delete it
                    {
                        remove = true;
                    }
                }
                if (remove)
                {
                    psms.Remove(psms[i]);
                    i--;
                }
            }
        }

        private static void PopulateSequenceLookUpDictionaries(string databaseFileName, List<Protein> proteins)
        {
            string[] array = databaseFileName.Split('\\');
            string filename = databaseFileName + "_NeoIndex\\NeoIndex_" + array[array.Length - 1] + ".txt";
            //index is ; separated with subsequence;Nsequence;Csequence;protaccession with internal delimited by _
            //the subsequence is removed from Nsequence and Csequence to preserve memory
            if (File.Exists(filename))
            {
                Dictionary<string, Protein> idToSequence = new Dictionary<string, Protein>();
                foreach (Protein prot in proteins)
                    idToSequence.Add(prot.Accession, prot);
                //Load existing index
                string[] index = File.ReadAllLines(filename);
                Parallel.ForEach(index, s =>
                {
                    string[] line = s.Replace("_;", ";").Split(';').ToArray();
                    string key = line[0];
                    List<string> nList = line[1].Split('_').ToList();
                    List<string> cList = line[2].Split('_').ToList();
                    if (nList.Count > 1 || nList[0].Length != 0)
                    {
                        for (int i = 0; i < nList.Count; i++)
                            nList[i] = key + nList[i];

                        List<Protein> prots = new List<Protein>();
                        string[] accessions = line[3].Split('|').ToArray();
                        for (int i = 0; i < accessions.Length - 1; i++)
                            prots.Add(idToSequence[accessions[i]]);
                        lock (ProtDictionary)
                        {
                            ProtDictionary.Add(line[0], prots);
                            NTermDictionary.Add(key, nList);
                        }
                    }
                    if (cList.Count > 1 || cList[0].Length != 0)
                    {
                        string cKey = key[3] + key.Substring(0, 3);
                        for (int i = 0; i < cList.Count; i++)
                            cList[i] = cList[i] + cKey;
                        lock (CTermDictionary)
                            CTermDictionary.Add(cKey, cList);
                    }
                });
            }
            else
            {
                Directory.CreateDirectory(databaseFileName + "_NeoIndex");
                //make a new index
                using (StreamWriter file = new StreamWriter(filename))
                {
                    foreach (char aa1 in AANames)
                    {
                        foreach (char aa2 in AANames)
                        {
                            foreach (char aa3 in AANames)
                            {
                                string threeMer = aa1.ToString() + aa2.ToString() + aa3.ToString();
                                List<string> nEntry = new List<string>();
                                List<string> cEntry = new List<string>();
                                List<Protein> protEntry = new List<Protein>();
                                Parallel.ForEach(proteins, protein =>
                                {
                                    List<string> localNEntry = new List<string>();
                                    List<string> localCEntry = new List<string>();

                                    int index = protein.BaseSequence.IndexOf(threeMer);
                                    while (index != -1)
                                    {
                                        localNEntry.Add(protein.BaseSequence.Length - index > 50 ? protein.BaseSequence.Substring(index, 50) : protein.BaseSequence.Substring(index, protein.BaseSequence.Length - index));
                                        localCEntry.Add(index + 3 > 50 ? protein.BaseSequence.Substring(index - 47, 50) : protein.BaseSequence.Substring(0, index + 3));
                                        index = protein.BaseSequence.IndexOf(threeMer, index + 1);
                                    }
                                    lock (nEntry)
                                    {
                                        foreach (string s in localNEntry)
                                            nEntry.Add(s);
                                        foreach (string s in localCEntry)
                                            cEntry.Add(s);
                                        if (localNEntry.Count != 0)
                                            protEntry.Add(protein);
                                    }
                                });

                                foreach (char aa4 in AANames)
                                {
                                    string nFourMer = threeMer + aa4.ToString();
                                    string cFourMer = aa4.ToString() + threeMer;
                                    List<string> localNEntry = nEntry.AsParallel().Where(seq => seq.Length > 3 && seq.Substring(0, 4).Equals(nFourMer)).ToList();
                                    List<string> LocalCEntry = cEntry.AsParallel().Where(seq => seq.Length > 3 && seq.Substring(seq.Length - 4, 4).Equals(cFourMer)).ToList();
                                    List<Protein> localProtEntry = protEntry.AsParallel().Where(prot => prot.BaseSequence.Contains(nFourMer)).ToList();
                                    NTermDictionary.Add(nFourMer, localNEntry);
                                    CTermDictionary.Add(cFourMer, LocalCEntry);
                                    ProtDictionary.Add(nFourMer, localProtEntry);
                                    file.Write(nFourMer + ";");
                                    localNEntry.ForEach(s => file.Write(s.Substring(4, s.Length - 4) + "_"));
                                    file.Write(";");
                                    LocalCEntry.ForEach(s => file.Write(s.Substring(0, s.Length - 4) + "_"));
                                    file.Write(";");
                                    localProtEntry.ForEach(prot => file.Write(prot.Accession + "|"));
                                    file.WriteLine(";");
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}