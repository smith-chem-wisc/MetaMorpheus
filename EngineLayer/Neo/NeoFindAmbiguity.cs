using Proteomics;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using Chemistry;
using System.Threading.Tasks;

namespace EngineLayer.Neo
{
    public class NeoFindAmbiguity
    {
        private BackgroundWorker worker = null;
        public List<Parent> theoreticalProteins;
        public HashSet<string> notFoundSequences = new HashSet<string>();
        public Dictionary<string, List<Parent>> foundSequences = new Dictionary<string, List<Parent>>();
        private const int maxMissingConsecutivePeaks = 2;
        private const int maxNumPossibleSequences = 2000;
        private const int decimalDigitsForFragmentMassRounding = 3;
        public static Dictionary<double, string[]> massDict = new Dictionary<double, string[]>();
        public static double[] keys;
        public List<List<string>> fragmentSeqs = new List<List<string>>();
        public static double productMassTolerancePpm = 20; //(Ppm)
        public static double precursorMassTolerancePpm = 5; //(Ppm)
        public static List<ProductType> ionsUsed = new List<ProductType> { ProductType.B, ProductType.Y };
        public static Dictionary<double, char> massesToResidues = new Dictionary<double, char>();
        public static List<double> singleAminoAcidMasses = new List<double>();
        public static double maxDifference;
        private readonly double waterMonoisotopicMass = Math.Round(PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass, decimalDigitsForFragmentMassRounding);
        public static HashSet<char> N2Position = new HashSet<char> { 'A', 'E', 'F', 'H', 'I', 'L', 'Q', 'R', 'T', 'V' };
        public static HashSet<char> C1Position = new HashSet<char> { 'F', 'I', 'L', 'M', 'W', 'Y' };
        public static Dictionary<string, List<string>> nTermDictionary = new Dictionary<string, List<string>>();
        public static Dictionary<string, List<string>> cTermDictionary = new Dictionary<string, List<string>>();
        public static Dictionary<string, List<Parent>> protDictionary = new Dictionary<string, List<Parent>>();

        public void FindAmbiguity(List<NeoPsm> candidates, List<Parent> theoreticalProteins, Ms2ScanWithSpecificMass[] spectra, out string error_message)
        {
            error_message = "";
            this.theoreticalProteins = theoreticalProteins;
            for (int i = 0; i < candidates.Count(); i++) //must be mutable while iterating
            {
                NeoPsm psm = candidates[i];

                Ms2ScanWithSpecificMass spectrum = spectra[psm.scanNumber];
                psm.fusionType = FusionCandidate.FusionType.TS; //for some maddening reason, this is not arriving here as trans, but instead translated
                if (IsTooMessy(psm, spectrum, out string e)) //having explosion of combinations when greater than 3 consequtive peaks producing tens of thousands of sequences ids, causes hanging
                {
                    candidates.RemoveAt(i);
                    i--;
                }
                else
                {
                    if (!GeneratePossibleSequences(psm, spectrum, out string error_message1)) //return true if fewer than specified number of ambiguities
                    {
                        candidates.RemoveAt(i);
                        i--;
                    }
                }
                this.worker.ReportProgress((int)((Convert.ToDouble(i > 0 ? i : 0) / candidates.Count()) * 100));
            }
            nTermDictionary.Clear();
            cTermDictionary.Clear();
            protDictionary.Clear();
        }

        public bool IsTooMessy(NeoPsm psm, Ms2ScanWithSpecificMass spectrum, out string error_message) //return true if too messy for confident identification
        {
            error_message = "";
            List<string> baseSequences = new List<string>();
            int currentBestScore = 0;
            for (int index = 0; index < psm.candidates.Count(); index++)
            {
                if (psm.candidates[index].seq.Length < 5)
                {
                    psm.candidates.RemoveAt(index);
                    index--;
                    continue;
                }
                bool badID = false;
                FusionCandidate fc = psm.candidates[index];
                FindIons(fc, psm, spectrum, out string error_message1);
                error_message += error_message1;
                int consecutiveMissedCounter = 0;
                int totalHitCounter = 0;
                if (fc.foundIons.Count(x => x) * 2 < fc.foundIons.Length)
                {
                    psm.candidates.RemoveAt(index);
                    index--;
                    continue;
                }
                for (int b = 0; b < fc.foundIons.Length; b++)
                {
                    if (consecutiveMissedCounter > maxMissingConsecutivePeaks) //if too many permutations possible because of an unmapped region
                        badID = true;
                    else if (!fc.foundIons[b])
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
                        psm.candidates.Remove(psm.candidates[0]);
                        index--;
                    }
                    currentBestScore = totalHitCounter;
                    baseSequences = new List<string> { psm.candidates[index].seq };
                }
                else if (totalHitCounter < currentBestScore || badID || baseSequences.Contains(psm.candidates[index].seq))
                {
                    psm.candidates.Remove(psm.candidates[index]);
                    index--;
                }
            }
            //If there's anything left
            if (psm.candidates.Count() > 0) //It wasn't too messy! Yay!
            {
                List<FusionCandidate> nTermCands = new List<FusionCandidate>();
                foreach (FusionCandidate oldCand in psm.candidates)
                {
                    oldCand.foundIons[1] = false;
                    oldCand.foundIons[2] = false;
                    if (!oldCand.foundIons[3])
                    {

                        FusionCandidate fc = new FusionCandidate(oldCand.seq);
                        FindIons(fc, psm, spectrum, out string error);
                        fc.foundIons[1] = true;
                        fc.foundIons[2] = false;
                        //fc.foundIons[3] = false;
                        if (!oldCand.foundIons[4])
                            fc.foundIons[4] = true;
                        nTermCands.Add(fc);
                        fc = new FusionCandidate(oldCand.seq);
                        FindIons(fc, psm, spectrum, out error);
                        fc.foundIons[1] = false;
                        fc.foundIons[2] = true;
                        //fc.foundIons[3] = false;
                        nTermCands.Add(fc);
                        fc = new FusionCandidate(oldCand.seq);
                        FindIons(fc, psm, spectrum, out error);
                        fc.foundIons[1] = false;
                        fc.foundIons[2] = false;
                        fc.foundIons[3] = true;
                        nTermCands.Add(fc);
                    }
                }
                foreach (FusionCandidate cand in nTermCands)
                    psm.candidates.Add(cand);
                return false;
            }
            else //this might be a fusion peptide, but we won't get enough valuable information from this spectra, so discard it
                return true;
        }

        //use ion hits to know where peaks have been found by morpheus and where there is ambiguity
        public static void FindIons(FusionCandidate fusionCandidate, NeoPsm psm, Ms2ScanWithSpecificMass spectrum, out string error_message)
        {
            error_message = "";
            //double[] nPeaks = psm.nInfo.peakHits; //get peaks
            //double[] cPeaks = psm.cInfo.peakHits;
            fusionCandidate.makeFoundIons();
            string candSeq = fusionCandidate.seq;
            bool[] foundIons = fusionCandidate.foundIons;

            //find which aa have peaks
            for (int i = 0; i < foundIons.Count() - 1; i++)
            {
                //B IONS//
                if (ionsUsed.Contains(ProductType.B))
                {
                    double bTheoMass = NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i)) - NeoConstants.WATER_MONOISOTOPIC_MASS;
                    //foreach (PTM ptm in psm.nInfo.ptms)
                    //{
                    //    if (ptm.index <= i)
                    //        bTheoMass += ptm.mass;
                    //}
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, bTheoMass, productMassTolerancePpm))
                            foundIons[i + 1] = true;
                    }
                }
                //Y IONS//
                if (ionsUsed.Contains(ProductType.Y))
                {
                    double yTheoMass = NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1));
                    //foreach (PTM ptm in psm.cInfo.ptms)
                    //{
                    //    if (ptm.index >= candSeq.Length - 2 - i)
                    //        yTheoMass += ptm.mass;
                    //}
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, yTheoMass, productMassTolerancePpm))
                            foundIons[foundIons.Count() - 1 - i] = true;
                    }
                }
                //C IONS//
                if (ionsUsed.Contains(ProductType.C))
                {
                    double cTheoMass = NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(0, 1 + i)) - NeoConstants.WATER_MONOISOTOPIC_MASS + NeoConstants.nitrogenMonoisotopicMass + 3 * NeoConstants.hydrogenMonoisotopicMass;
                    //foreach (PTM ptm in psm.nInfo.ptms)
                    //{
                    //    if (ptm.index <= i)
                    //        cTheoMass += ptm.mass;
                    //}
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, cTheoMass, productMassTolerancePpm))
                            foundIons[i + 1] = true;
                    }
                }
                //ZDOT IONS//
                if (ionsUsed.Contains(ProductType.Zdot))
                {
                    double zdotTheoMass = NeoMassCalculator.MonoIsoptopicMass(candSeq.Substring(candSeq.Length - 1 - i, i + 1)) - NeoConstants.nitrogenMonoisotopicMass - 2 * NeoConstants.hydrogenMonoisotopicMass;
                    //foreach (PTM ptm in psm.cInfo.ptms)
                    //{
                    //    if (ptm.index >= candSeq.Length - 2 - i)
                    //        zdotTheoMass += ptm.mass;
                    //}
                    foreach (double expPeak in spectrum.TheScan.MassSpectrum.XArray)
                    {
                        if (NeoMassCalculator.IdenticalMasses(expPeak, zdotTheoMass, productMassTolerancePpm))
                            foundIons[foundIons.Count() - 1 - i] = true;
                    }
                }
            }
            //foundIons[0] = true; //AspN always starts with a D
            foundIons[0] = true;//|A|B|C|D|E|F|K where the whole peptide peak is always placed arbitrarily at the n term


            //Look for mass shifts that may be indicative of internal ions to resolve ambiguity
            char[] originalSeq = fusionCandidate.seq.ToCharArray();
            //for (int i = 0; i < foundIons.Length; i++)
            //{
            //    if (!foundIons[i])
            //    {
            //        List<Edge> edges = spectrum.edges[originalSeq[i]];
            //        if (i != 0)
            //            foreach (Edge edge1 in edges)
            //                foreach (Edge edge2 in edge1.lowMassNode.lowerMassEdges)
            //                    if (edge2.aminoAcid == originalSeq[i - 1])
            //                        foundIons[i] = true;

            //        if (i != foundIons.Length - 1)
            //            foreach (Edge edge1 in edges)
            //                foreach (Edge edge2 in edge1.highMassNode.higherMassEdges)
            //                    if (edge2.aminoAcid == originalSeq[i + 1])
            //                        foundIons[i] = true;
            //    }
            //}

            //if (!N2Position.Contains(originalSeq[1]))
            //{
            //    if (N2Position.Contains(originalSeq[1]))
            //        foundIons[1] = false;
            //    else if (N2Position.Contains(originalSeq[2]))
            //        foundIons[2] = false;
            //}
        }

        public static void ReadMassDictionary()
        {
            List<double> tempKeys = new List<double>();

            using (StreamReader masses = new StreamReader(Environment.CurrentDirectory + @"\Data\Dictionary" + maxMissingConsecutivePeaks + ".txt")) //file located in Morpheus folder
            {
                while (masses.Peek() != -1)
                {
                    string line = masses.ReadLine();
                    string[] fields = line.Split('\t');
                    double key = Convert.ToDouble(fields[0]);
                    string[] sequences = fields[1].Split(';');
                    massDict.Add(key, sequences);
                    tempKeys.Add(key);
                }
            }
            keys = new double[tempKeys.Count()];
            for (int i = 0; i < tempKeys.Count(); i++)
            {
                keys[i] = tempKeys[i];
            }
        }

        public static void ReadAminoAcids()
        {
            for (int i = 0; i < Residue.ResidueMonoisotopicMass.Length; i++)
            {
                if (!double.IsNaN(Residue.ResidueMonoisotopicMass[i]) && !singleAminoAcidMasses.Contains(Residue.ResidueMonoisotopicMass[i]))
                {
                    massesToResidues.Add(Residue.ResidueMonoisotopicMass[i], (char)i);
                    singleAminoAcidMasses.Add(Residue.ResidueMonoisotopicMass[i]);
                }
            }
            singleAminoAcidMasses.Sort();
            maxDifference = singleAminoAcidMasses[singleAminoAcidMasses.Count - 1];
        }

        public bool GeneratePossibleSequences(NeoPsm psm, Ms2ScanWithSpecificMass spectrum, out string error_message) //returns false if over the specified number of sequences are generated
        {
            error_message = "";
            List<string> ambiguousCandidates = new List<string>();
            foreach (FusionCandidate fusionCandidate in psm.candidates)
            {
                List<string> nMatchedSequences = new List<string>();
                List<string> cMatchedSequences = new List<string>();
                List<string> partialSequences = new List<string> { "" };
                double theoreticalMass = NeoMassCalculator.MonoIsoptopicMass(fusionCandidate.seq);
                //FindIons(fusionCandidate, psm, spectrum, out string error_message1); //populate the foundIons array
                //error_message += error_message1;

                int globalIndex = 0;
                //C-Terminus
                //Increment through the amino acids to determine all possible sequences
                //once a sequence is no longer findable (is not in database), then stop generating possible sequences for that sequence
                //stop when reaching end of peptide or no more present in database
                while (fusionCandidate.foundIons.Length > globalIndex)
                {
                    int reverseIndex = fusionCandidate.foundIons.Length - globalIndex - 1;
                    int mostRecent = fusionCandidate.foundIons.Length; //most recent Ion found prior to this one (start point)
                    if (fusionCandidate.foundIons[reverseIndex])
                    {
                        for (int i = reverseIndex + 1; i < fusionCandidate.foundIons.Length; i++) //identify start point
                            if (fusionCandidate.foundIons[i])
                            {
                                mostRecent = i; //save most recent hit, exclusive of the current index
                                break;
                            }

                        //get combos
                        string ambiguousFrag = fusionCandidate.seq.Substring(reverseIndex, mostRecent - reverseIndex);
                        double key = NeoMassCalculator.MonoIsoptopicMass(ambiguousFrag);

                        List<string> combinations = GetCombinations(key, theoreticalMass - key, psm.expMass);

                        List<string> tempPartialSequences = new List<string>();
                        foreach (string old in partialSequences)
                            foreach (string neu in combinations)
                                tempPartialSequences.Add(neu + old);

                        partialSequences.Clear();

                        //for (int i = tempPartialSequences.Count - 1; i >= 0; i--)
                        //{
                        //    string seq = tempPartialSequences[i];
                        //    if (seq.Length > 0 && !C1Position.Contains(seq[seq.Length - 1]))
                        //        tempPartialSequences.RemoveAt(i);
                        //}
                        if (tempPartialSequences.Count == 0)
                            break;

                        for (int i = tempPartialSequences.Count - 1; i >= 0; i--)
                        {
                            string partialSeq = tempPartialSequences[i];
                            if (partialSeq.Length > 3)
                            {
                                //can it be made?
                                if (cTermDictionary.TryGetValue(partialSeq.Substring(partialSeq.Length - 4, 4), out List<string> cEntry))
                                {
                                    int n = 4;
                                    while (cEntry.Count != 0 && n <= partialSeq.Length)
                                    {
                                        string currentSeq = partialSeq.Substring(partialSeq.Length - n, n);
                                        if (!cMatchedSequences.Any(x => x.Length >= n && x.Substring(x.Length - n, n).Equals(currentSeq)))
                                            cEntry = cEntry.AsParallel().Where(seq => seq.Length >= n && seq.Substring(seq.Length - n, n).Equals(currentSeq)).ToList();
                                        n++;
                                    }
                                    if (cEntry.Count == 0)
                                    {
                                        string strToAdd = partialSeq.Substring(partialSeq.Length - n + 2, n - 2);
                                        if (!cMatchedSequences.Any(x => x.Equals(strToAdd)))
                                        {
                                            for (int j = cMatchedSequences.Count - 1; j >= 0; j--)
                                                if (strToAdd.Contains(cMatchedSequences[j]))
                                                    cMatchedSequences.RemoveAt(j);
                                            if (!cMatchedSequences.Any(x => x.Contains(strToAdd)))
                                                cMatchedSequences.Add(strToAdd);
                                        }
                                    }
                                    else
                                        partialSequences.Add(partialSeq); //keep adding stuff to it
                                }
                                else { }
                            }
                            else
                                partialSequences.Add(partialSeq); //keep adding stuff to it

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
                while (fusionCandidate.foundIons.Length >= globalIndex)
                {
                    int mostRecent = 0; //most recent Ion found prior to this one (start point)
                    if (fusionCandidate.foundIons.Length == globalIndex || fusionCandidate.foundIons[globalIndex])
                    {
                        for (int i = globalIndex - 1; i >= 0; i--) //identify start point
                            if (fusionCandidate.foundIons[i])
                            {
                                mostRecent = i; //save most recent hit, exclusive of the current index
                                break;
                            }

                        //get combos
                        string ambiguousFrag = fusionCandidate.seq.Substring(mostRecent, globalIndex - mostRecent);
                        double key = NeoMassCalculator.MonoIsoptopicMass(ambiguousFrag);

                        List<string> combinations = GetCombinations(key, theoreticalMass - key, psm.expMass);

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

                                if (nTermDictionary.TryGetValue(partialSeq.Substring(0, 4), out List<string> nEntry))
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
                            if (ambiguousCandidates.Count > maxNumPossibleSequences)
                                break;
                        }
                        double cMass = NeoMassCalculator.MonoIsoptopicMass(c);
                        double totalMass = nMass + cMass - NeoConstants.WATER_MONOISOTOPIC_MASS;
                        if (totalMass + 1 > psm.expMass)
                        {
                            if (NeoMassCalculator.IdenticalMasses(psm.expMass, totalMass, precursorMassTolerancePpm))
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
                                    if (NeoMassCalculator.IdenticalMasses(psm.expMass, totalMass, precursorMassTolerancePpm))
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
                                    else if (totalMass + 10 < psm.expMass)
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


            psm.candidates.Clear();
            if (ambiguousCandidates.Count <= maxNumPossibleSequences)
            {
                //CIS CODE
                foreach (string fc in ambiguousCandidates)
                {
                    if (fc.Length < 5)
                        continue;
                    //if(fc.Equals("TVDDKIVVQKY")||fc.Equals("AEVVVKEGL"))
                    //{ }
                    string substring = fc.Substring(0, 4);
                    //  bool translated = false;
                    List<string> nFrag;
                    nTermDictionary.TryGetValue(fc.Substring(0, 4), out nFrag);
                    if (nFrag != null && nFrag.AsParallel().Any(seq => seq.Length >= fc.Length && seq.Substring(0, fc.Length).Equals(fc)))
                    {
                        psm.fusionType = FusionCandidate.FusionType.TL;
                        psm.candidates.Add(new FusionCandidate(fc) { fusionType = FusionCandidate.FusionType.TL });
                    }
                    else
                    {
                        //Is it cis?
                        bool cis = false;
                        List<Parent> otherPossibleProteins;
                        if (protDictionary.TryGetValue(substring, out List<Parent> possibleProteins))
                        {
                            //possibleProteins.ForEach(prot => originalPossibleProteins.Add(prot));
                            for (int i = 5; i <= fc.Length; i++)
                            {
                                string otherSubstring = fc.Substring(i - 1, fc.Length - i + 1);
                                //get proteins containing both halves
                                otherPossibleProteins = possibleProteins.AsParallel().Where(prot => prot.seq.Contains(otherSubstring)).ToList();

                                //check if both halves are the correct distance apart and are not overlapping
                                foreach (Parent prot in otherPossibleProteins)
                                {
                                    string seq = prot.seq;
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
                                            psm.candidates.Add(new FusionCandidate(fc) { fusionType = FusionCandidate.FusionType.NC });
                                            cis = true;
                                            if (psm.fusionType == FusionCandidate.FusionType.TS)
                                                psm.fusionType = FusionCandidate.FusionType.NC;
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
                                possibleProteins = possibleProteins.AsParallel().Where(prot => prot.seq.Contains(substring)).ToList();
                                if (possibleProteins.Count == 0)
                                    break;
                            }
                        }
                        if (!cis)
                        {
                            if (protDictionary.TryGetValue(fc.Substring(fc.Length - 4, 4), out otherPossibleProteins))
                            {
                                string otherSubstring = fc.Substring(3, fc.Length - 3);
                                otherPossibleProteins = otherPossibleProteins.Where(prot => prot.seq.Contains(otherSubstring)).ToList();
                                for (int i = 2; i >= 0; i--)
                                {
                                    substring = fc.Substring(0, i + 1);
                                    foreach (Parent prot in otherPossibleProteins)
                                    {
                                        string seq = prot.seq;
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
                                                psm.candidates.Add(new FusionCandidate(fc) { fusionType = FusionCandidate.FusionType.NC });
                                                cis = true;
                                                if (psm.fusionType == FusionCandidate.FusionType.TS)
                                                    psm.fusionType = FusionCandidate.FusionType.NC;
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
                                psm.candidates.Add(new FusionCandidate(fc));
                        }
                    }
                }

                //if (psm.fusionType != FusionCandidate.FusionType.TS)
                //    for (int i = psm.candidates.Count - 1; i >= 0; i--)
                //        if (psm.candidates[i].fusionType == FusionCandidate.FusionType.TS)
                //            psm.candidates.RemoveAt(i);
            }
            //Old Code: Simply add all ambiguous strings to psm.candidates
            return psm.candidates.Count > 0;
        }

        //returns false if a full fusion sequence could not be made or was found in the database, making it translated instead of a novel fusion.
        public bool PossibleCandidate(NeoPsm psm)
        {
            foundSequences = new Dictionary<string, List<Parent>>(); //used for finding longer fragments than those previously identified. Also populates ParentInfo
            notFoundSequences = new HashSet<string>(); //don't bother looking for these fragments, since we know they don't exist. Good for multiple homologous putative fusion peptide sequences

            //conduct an initial search of each candidate's full sequence to identify any that are translated
            for (int i = 0; i < psm.candidates.Count(); i++) //foreach fusion peptide sequence that could map to this scan
            {
                string novelSeq = psm.candidates[i].seq;
                if (foundParent(novelSeq, ParentInfo.terminal.C, psm.candidates[i], false)) //check really quick to see if the whole thing exists as is. If so, assign it as translated. Terminal C was arbitrarily chosen
                {
                    foreach (ParentInfo info in psm.candidates[i].parentInfo)
                        foreach (Parent protein in info.theoreticalProteins)
                            if (protein.seq.Contains(novelSeq)) //if translated
                                psm.candidates[i].translatedParents.Add(new TranslatedParent(protein.id, protein.seq, protein.seq.IndexOf(psm.candidates[i].seq), psm.candidates[i].seq.Length));

                    psm.candidates[i].fusionType = FusionCandidate.FusionType.TL;
                    psm.fusionType = psm.candidates[i].fusionType;
                    for (int j = 0; j < psm.candidates.Count(); j++)
                    {
                        if (j != i)
                        {
                            psm.candidates.Remove(psm.candidates[j]);
                            j--;
                            i--;
                        }
                    }
                    return false;
                }
            }
            for (int i = 0; i < psm.candidates.Count(); i++) //foreach fusion peptide sequence that could map to this scan
            {
                //sw.StartFindParents
                if (!isViable(psm.candidates[i])) //remove this fusion peptide sequence if the parent fragments cannot be found with the given database
                {
                    psm.candidates.Remove(psm.candidates[i]);
                    i--;
                }
                else
                {
                    DetermineFusionCandidateType(psm.candidates[i]); //cis, trans?
                    if (psm.fusionType > psm.candidates[i].fusionType) //if more likely than previous types, change the psm type (golf scoring)
                    {
                        psm.fusionType = psm.candidates[i].fusionType;
                    }
                    if (psm.fusionType.Equals(FusionCandidate.FusionType.TL))//if there's a possible sequence that's present in the database, it is likely correct and is it is not worth it to identify parents of other sequences.
                    {
                        //remove all other candidates
                        for (int j = 0; j < psm.candidates.Count(); j++)
                        {
                            if (j != i)
                            {
                                psm.candidates.Remove(psm.candidates[j]);
                                j--;
                                i--;
                            }
                        }
                        return false;
                    }
                }
            }

            return psm.candidates.Count() != 0; //if no candidates are left, we couldn't make the sequence with the database and we'll discard it.

        }

        public List<string> GetCombinations(double key, double theoreticalMass, double experimentalMass)
        {
            List<string> combinations = new List<string>();
            double closestPeak = double.NaN;
            var ipos = Array.BinarySearch(keys, key);
            if (ipos < 0)
                ipos = ~ipos;

            if (ipos > 0)
            {
                var downIpos = ipos - 1;
                // Try down
                while (downIpos >= 0)
                {
                    closestPeak = keys[downIpos];
                    if (NeoMassCalculator.IdenticalMasses(experimentalMass, theoreticalMass + closestPeak, precursorMassTolerancePpm))
                    {
                        if (massDict.TryGetValue(closestPeak, out string[] value))
                            foreach (string frag in value)
                                combinations.Add(frag);
                    }
                    else
                        break;
                    downIpos--;
                }
            }
            if (ipos < keys.Length)
            {
                var upIpos = ipos;
                // Try here and up
                while (upIpos < keys.Length)
                {
                    closestPeak = keys[upIpos];
                    if (NeoMassCalculator.IdenticalMasses(experimentalMass, theoreticalMass + closestPeak, precursorMassTolerancePpm))
                    {
                        if (massDict.TryGetValue(closestPeak, out string[] value))
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


        public bool isViable(FusionCandidate tempCandidate) //returns if sequence could be made from one or two proteins in database and writes fusion type, parents, and junctions to fusionCandidate
        {
            //need to check that each index is viable
            string novelSeq = tempCandidate.seq;
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
            if (foundParent(testFrag, ParentInfo.terminal.N, tempCandidate, foundFirstSearch)) //if found
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
                if (foundParent(testFrag, ParentInfo.terminal.N, tempCandidate, foundFirstSearch)) //if found
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
            if (foundParent(testFrag, ParentInfo.terminal.C, tempCandidate, foundFirstSearch)) //if found
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
                if (foundParent(testFrag, ParentInfo.terminal.C, tempCandidate, foundFirstSearch))
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
            {
                /*       foreach (ParentInfo PI in nTermPI) //append B searches 
                           {
                               tempCandidate.parentInfo.Add(PI);
                           }*/
                //Add Result
            }
            if (cTermParentLength + nTermParentLength < novelSeq.Length) //if no overlap
            {
                return false;
            }
            else
            {
                for (int junction = tempCandidate.seq.Length - cTermParentLength - 1; junction < nTermParentLength; junction++)
                {
                    tempCandidate.addJunctionIndex(junction);
                }
                return true;
            }
        }

        public bool foundParent(string frag, ParentInfo.terminal terminal, FusionCandidate candidate, bool foundFirstSearch)
        {
            //localTheoreticals.AsParallel().Where(x => x.Contains(frag)).ToList();
            if (notFoundSequences.Contains(frag)) //has the fragment been searched but not found before?
                return false;

            List<Parent> matches = new List<Parent>();
            if (foundSequences.TryGetValue(frag, out matches)) //has the fragment been searched AND found before?
            {
                candidate.parentInfo.Add(new ParentInfo(matches, terminal, frag));
                return true;
            }

            if (foundFirstSearch) //Has something smaller been found before? Well, then we can just search against those found sequences
            {
                string shorterFrag = terminal.Equals(ParentInfo.terminal.N) ? frag.Substring(0, frag.Length - 1) : frag.Substring(1, frag.Length - 1);

                foreach (ParentInfo info in candidate.parentInfo)
                {
                    if (info.parentType.Equals(terminal) && info.fragFound.Equals(shorterFrag))
                    {
                        List<Parent> tempProtList = new List<Parent>();
                        info.theoreticalProteins.ForEach(protein => tempProtList.Add(protein));
                        matches = tempProtList.AsParallel().Where(x => x.seq.Contains(frag)).ToList();
                    }
                }
            }
            else //it hasn't been found before... we need to search against the whole database :(
            {
                matches = theoreticalProteins.AsParallel().Where(x => x.seq.Contains(frag)).ToList();
            }
            if (matches != null && matches.Count() > 0)
            {
                foundSequences.Add(frag, matches);
                candidate.parentInfo.Add(new ParentInfo(matches, terminal, frag));
                return true;
            }
            else
            {
                notFoundSequences.Add(frag);
                return false;
            }
        }

        public void DetermineFusionCandidateType(FusionCandidate fusionCandidate)
        {
            if (!fusionCandidate.fusionType.Equals(FusionCandidate.FusionType.TL))
            {
                string sequence = fusionCandidate.seq;
                foreach (ParentInfo info in fusionCandidate.parentInfo)
                {
                    int foundLength = info.fragFound.Length;

                    string compFrag = "";// fusionCandidate.seq.Substring()
                    if (info.parentType.Equals(ParentInfo.terminal.N))
                    {
                        compFrag = fusionCandidate.seq.Substring(foundLength, fusionCandidate.seq.Length - foundLength);
                    }
                    else
                    {
                        compFrag = fusionCandidate.seq.Substring(0, fusionCandidate.seq.Length - foundLength);
                    }

                    foreach (Parent protein in info.theoreticalProteins)
                    {
                        //get the index(es) of where the found fragment is
                        string subProt = protein.seq;
                        List<int> originalIndexes = new List<int>();
                        int pastIndex = 0;
                        while (subProt.Contains(info.fragFound))
                        {
                            int newIndex = subProt.IndexOf(info.fragFound);
                            originalIndexes.Add(newIndex + pastIndex);
                            subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                            pastIndex += newIndex + 1;
                        }
                        fusionCandidate.transParents.Add(new TransParent(protein.id, protein.seq, originalIndexes, foundLength, info.parentType));

                        //get the index(es) of where the complimentary fragment is (if it's a cis fusion peptide)
                        subProt = protein.seq;
                        List<int> complementaryIndexes = new List<int>();
                        pastIndex = 0;
                        while (subProt.Contains(compFrag))
                        {
                            int newIndex = subProt.IndexOf(compFrag);
                            complementaryIndexes.Add(newIndex + pastIndex);
                            subProt = subProt.Substring(newIndex + 1, subProt.Length - newIndex - 1); //need to remove old match
                            pastIndex += newIndex + 1;
                        }
                        if (complementaryIndexes.Count() > 0) //if it is not trans
                        {
                            //if it is cis
                            if (info.parentType.Equals(ParentInfo.terminal.N))
                                fusionCandidate.cisParents.Add(new CisParent(protein.id, protein.seq, originalIndexes, foundLength, complementaryIndexes, compFrag.Length));
                            else
                                fusionCandidate.cisParents.Add(new CisParent(protein.id, protein.seq, complementaryIndexes, compFrag.Length, originalIndexes, foundLength));
                        }
                    }
                }
            }
            else
            {
                string seq = fusionCandidate.seq;
                foreach (ParentInfo info in fusionCandidate.parentInfo)
                {
                    foreach (Parent protein in info.theoreticalProteins)
                    {
                        if (protein.seq.Contains(seq)) //if translated
                        {
                            fusionCandidate.translatedParents.Add(new TranslatedParent(protein.id, protein.seq, protein.seq.IndexOf(fusionCandidate.seq), fusionCandidate.seq.Length));
                        }
                    }
                }
            }
            foreach (CisParent cisParent in fusionCandidate.cisParents)
            {
                if (cisParent.cisType < fusionCandidate.fusionType)
                {
                    fusionCandidate.fusionType = cisParent.cisType;
                }
            }

        }

        public void RemoveTranslatablePeptides(List<NeoPsm> psms)
        {
            for (int i = 0; i < psms.Count(); i++)
            {
                bool remove = false;
                foreach (FusionCandidate fc in psms[i].candidates)
                {
                    if (fc.junctionIndexes.Count() > fc.seq.Length) //delete it
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
    }
}

