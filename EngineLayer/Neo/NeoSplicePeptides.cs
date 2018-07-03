using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Neo
{
    public static class NeoSplicePeptides
    {
        public static readonly int IonsUsedMassVer = 2;
        public static double FixedModMass = 0.0;
        //carbamidomethyl is defaulted at 0.

        public static List<NeoPsm> SplicePeptides(List<NeoPsm> psms)
        {
            List<NeoPsm> candidates = new List<NeoPsm>();

            //int counter = 0;
            foreach (NeoPsm psm in psms)
            {
                //preliminary filters can be removed if MassMatch calls to IonCrop are set to true.
                string B = IonCrop(psm.NInfo.Seq, psm.ExpMass, 0, ProductType.B, true); //used as a preliminary filter to prevent longer searches from seq ID's that are larger than the precursor mass
                string Y = IonCrop(psm.CInfo.Seq, psm.ExpMass, 0, ProductType.Y, true); //used as a preliminary filter to prevent longer searches from seq ID's that are larger than the precursor mass
                for (int y = 0; y < Y.Length - IonsUsedMassVer; y++) //foreach y aa removed
                {
                    for (int b = 0; b < B.Length - IonsUsedMassVer; b++) //foreach b aa removed
                    {
                        MassMatch(B, Y, psm, b, y); //return a FusionCandidate
                    }
                }
                if (psm.Candidates.Count > 0) //match was found
                    candidates.Add(psm);
                else
                {
                    AutoFill(psm);
                    if (psm.Candidates.Count > 0)
                    {
                        candidates.Add(psm);
                    }
                }
            }
            //DetermineCommonAminoAcids();
            return candidates;
        }

        //method was originally written recursively, but large peptides result in stackoverflow exceptions
        public static void MassMatch(string B, string Y, NeoPsm psm, int BIndex, int YIndex) //this is the workhorse of SpliceFragments
        {
            double experimentalMass = psm.ExpMass;
            string bFrag = IonCrop(B, experimentalMass, BIndex, ProductType.B, false); //returns a B ion sequence that has a mass smaller than the experimental mass by cleaving C term AA
            if (bFrag.Length < IonsUsedMassVer)//If the number of AA from the N-term peptide is less than desired amount, start over loop and remove a single aa from the C-term
                return;
            string yFrag = IonCrop(Y, experimentalMass, YIndex, ProductType.Y, false); //returns a Y ion sequence that has a mass smaller than the experimental mass by cleaving N term AA
            if (yFrag.Length < IonsUsedMassVer)//If the number of AA from the C-term peptide is less than desired amount, end recursion.
                return;
            double theoreticalMass = NeoMassCalculator.MonoIsoptopicMass(bFrag) + NeoMassCalculator.MonoIsoptopicMass(yFrag) - NeoConstants.WATER_MONOISOTOPIC_MASS + FixedModMass; //water added once in b and once in y

            if (NeoMassCalculator.IdenticalMasses(experimentalMass, theoreticalMass, NeoFindAmbiguity.PrecursorMassTolerancePpm))//if match
            {
                string novelSeq = bFrag + yFrag;
                if (!psm.Candidates.Any(x => x.Seq.Equals(novelSeq))) //if fusion sequence was not previously assigned to this psm
                {
                    psm.Candidates.Add(new FusionCandidate(novelSeq));
                }
            }
        }

        public static void AutoFill(NeoPsm psm)
        {
            double experimentalMass = psm.ExpMass;
            for (int i = 0; i < 2; i++)
            {
                string ionSequence;
                ProductType ion;
                if (i == 0)
                {
                    ionSequence = psm.NInfo.Seq;
                    ion = ProductType.B;
                }
                else
                {
                    ionSequence = psm.CInfo.Seq;
                    ion = ProductType.Y;
                }

                int fragNumber = 0;
                string ionFrag = ionSequence;
                while (ionFrag.Length > 1)
                {
                    if (ion == ProductType.B)
                    {
                        ionFrag = ionSequence.Substring(0, (ionSequence.Length - fragNumber));
                        if (ionFrag.Substring(ionSequence.Length - fragNumber - 1, 1) == "]") //if end of a PTM annotation
                        {
                            while (ionFrag.Substring(ionSequence.Length - fragNumber - 1, 1) != "[")
                            {
                                fragNumber++;
                                ionFrag = ionSequence.Substring(0, (ionSequence.Length - fragNumber));
                            }
                            fragNumber += 2; //removes "(" and the AA the PTM was attached to
                            ionFrag = ionSequence.Substring(0, (ionSequence.Length - fragNumber));
                        }
                    }
                    else //Ion==Y
                    {
                        ionFrag = ionSequence.Substring((0 + fragNumber), (ionSequence.Length - fragNumber));
                        if (ionFrag.Substring(0, 1) == "[") //if start of a PTM annotation
                        {
                            while (ionFrag.Substring(0, 1) != "]")
                            {
                                fragNumber++;
                                ionFrag = ionSequence.Substring((0 + fragNumber), (ionSequence.Length - fragNumber));
                            }
                            fragNumber++; //removes ")"
                            ionFrag = ionSequence.Substring((0 + fragNumber), (ionSequence.Length - fragNumber));
                        }
                    }

                    double theoreticalMass = NeoMassCalculator.MonoIsoptopicMass(ionFrag) - NeoConstants.WATER_MONOISOTOPIC_MASS;
                    if (theoreticalMass < experimentalMass)
                    {
                        double key = experimentalMass - theoreticalMass;
                        List<string> combinations = new List<string>();

                        double closestPeak = double.NaN;
                        var ipos = Array.BinarySearch(NeoFindAmbiguity.Keys, key);
                        if (ipos < 0)
                            ipos = ~ipos;

                        if (ipos > 0)
                        {
                            var downIpos = ipos - 1;
                            // Try down
                            while (downIpos >= 0)
                            {
                                closestPeak = NeoFindAmbiguity.Keys[downIpos];
                                if (NeoMassCalculator.IdenticalMasses(experimentalMass, closestPeak + theoreticalMass, NeoFindAmbiguity.PrecursorMassTolerancePpm))
                                    foreach (string frag in NeoFindAmbiguity.MassDict[closestPeak])
                                        combinations.Add(frag);
                                else
                                    break;
                                downIpos--;
                            }
                        }
                        if (ipos < NeoFindAmbiguity.Keys.Length)
                        {
                            var upIpos = ipos;
                            // Try here and up
                            while (upIpos < NeoFindAmbiguity.Keys.Length)
                            {
                                closestPeak = NeoFindAmbiguity.Keys[upIpos];
                                if (NeoMassCalculator.IdenticalMasses(experimentalMass, closestPeak + theoreticalMass, NeoFindAmbiguity.PrecursorMassTolerancePpm))
                                    foreach (string frag in NeoFindAmbiguity.MassDict[closestPeak])
                                        combinations.Add(frag);
                                else
                                    break;
                                upIpos++;
                            }
                        }
                        if (combinations.Count != 0)
                        {
                            foreach (string s in combinations)
                            {
                                if (ion == ProductType.B)
                                    psm.Candidates.Add(new FusionCandidate(ionFrag + s));
                                else
                                    psm.Candidates.Add(new FusionCandidate(s + ionFrag));
                            }
                            break;
                        }
                    }
                    fragNumber++;
                }
            }
        }

        //This method removes the number of amino acids specified by FragNumber from the respecitve terminus specified by ion of IonSequence
        //If checkToRemoveExtraAA is true, additional AA will be removed to achieve a theoretical mass less than the Experimental mass
        public static string IonCrop(string ionSequence, double experimentalMass, int fragNumber, ProductType ion, bool checkToRemoveExtraAA)
        {
            do
            {
                string ionFrag;
                if (ion == ProductType.B)
                {
                    ionFrag = ionSequence.Substring(0, (ionSequence.Length - fragNumber));
                    if (ionFrag.Substring(ionSequence.Length - fragNumber - 1, 1) == "]") //if end of a PTM annotation
                    {
                        while (ionFrag.Substring(ionSequence.Length - fragNumber - 1, 1) != "[")
                        {
                            fragNumber++;
                            ionFrag = ionSequence.Substring(0, (ionSequence.Length - fragNumber));
                        }
                        fragNumber += 2; //removes "(" and the AA the PTM was attached to
                        ionFrag = ionSequence.Substring(0, (ionSequence.Length - fragNumber));
                    }
                }
                else //Ion==Y
                {
                    ionFrag = ionSequence.Substring((0 + fragNumber), (ionSequence.Length - fragNumber));
                    if (ionFrag.Substring(0, 1) == "[") //if start of a PTM annotation
                    {
                        while (ionFrag.Substring(0, 1) != "]")
                        {
                            fragNumber++;
                            ionFrag = ionSequence.Substring((0 + fragNumber), (ionSequence.Length - fragNumber));
                        }
                        fragNumber++; //removes ")"
                        ionFrag = ionSequence.Substring((0 + fragNumber), (ionSequence.Length - fragNumber));
                    }
                }
                if (!checkToRemoveExtraAA)
                    return ionFrag;
                else
                {
                    double IonMass = NeoMassCalculator.MonoIsoptopicMass(ionFrag);
                    if (NeoMassCalculator.IdenticalMasses(experimentalMass, IonMass, NeoFindAmbiguity.PrecursorMassTolerancePpm))
                    {
                        //Assume this ID is correct
                        return "";
                    }
                    else if (IonMass < experimentalMass) //end if the mass of the fragment is lower than the experimental
                        return ionFrag;
                    else //call the function again to remove another amino acid.
                        fragNumber++;
                }
            }
            while (true);
        }
    }
}