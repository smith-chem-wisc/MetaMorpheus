using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Neo
{
    public static class NeoSplicePeptides
    {
        public static readonly int ionsUsedMassVer = 2;
        public static double fixedModMass = 0.0; //carbamidomethyl is defaulted at 0.


        public static List<NeoPsm> SplicePeptides(List<NeoPsm> psms)
        {
            List<NeoPsm> candidates = new List<NeoPsm>();

            //int counter = 0;
            foreach (NeoPsm psm in psms)
            {
                //this.worker.ReportProgress(Convert.ToInt16((Convert.ToDouble(counter) / Convert.ToDouble(psms.Count())) * 100));
                //counter++;
                //preliminary filters can be removed if MassMatch calls to IonCrop are set to true.
                if (psm.cInfo.seq.Contains("AWRRSC"))
                { }
                string B = IonCrop(psm.nInfo.seq, psm.expMass, 0, ProductType.B, true); //used as a preliminary filter to prevent longer searches from seq ID's that are larger than the precursor mass

                if(psm.cInfo.seq.Contains("AWRRSC"))
                { }
                string Y = IonCrop(psm.cInfo.seq, psm.expMass, 0, ProductType.Y, true); //used as a preliminary filter to prevent longer searches from seq ID's that are larger than the precursor mass
                for (int y = 0; y < Y.Length - ionsUsedMassVer; y++) //foreach y aa removed
                {
                    for (int b = 0; b < B.Length - ionsUsedMassVer; b++) //foreach b aa removed
                    {
                        MassMatch(B, Y, psm, b, y); //return a FusionCandidate
                    }
                }
                if (psm.candidates.Count > 0) //match was found
                    candidates.Add(psm);
                else
                {
                    AutoFill(psm);
                    if (psm.candidates.Count > 0)
                        candidates.Add(psm);
                }
            }
            //DetermineCommonAminoAcids();
            return candidates;
        }

        //method was originally written recursively, but large peptides result in stackoverflow exceptions
        public static void MassMatch(string B, string Y, NeoPsm psm, int BIndex, int YIndex) //this is the workhorse of SpliceFragments
        {
            double experimentalMass = psm.expMass;
            string bFrag = IonCrop(B, experimentalMass, BIndex, ProductType.B, false); //returns a B ion sequence that has a mass smaller than the experimental mass by cleaving C term AA
            //BIndex = B.Length - BFrag.Length; //added 11/8/16 Useful first pass to record how many AA have been cleaved from C term 
            if (bFrag.Length < ionsUsedMassVer)//If the number of AA from the N-term peptide is less than desired amount, start over loop and remove a single aa from the C-term
                return;
            string yFrag = IonCrop(Y, experimentalMass, YIndex, ProductType.Y, false); //returns a Y ion sequence that has a mass smaller than the experimental mass by cleaving N term AA
            //YIndex = Y.Length - YFrag.Length; //added 11/8/16 Useful first pass to record how many AA have been cleaved from N term
            if (yFrag.Length < ionsUsedMassVer)//If the number of AA from the C-term peptide is less than desired amount, end recursion.
                return;
            double theoreticalMass = NeoMassCalculator.MonoIsoptopicMass(bFrag) + NeoMassCalculator.MonoIsoptopicMass(yFrag) - NeoConstants.WATER_MONOISOTOPIC_MASS + fixedModMass; //water added once in b and once in y

            //add PTM masses
            //foreach (PTM ptm in psm.nInfo.ptms)
            //{
            //    if (ptm.index < bFrag.Length)
            //    {
            //        theoreticalMass += ptm.mass;
            //    }
            //}
            //foreach (PTM ptm in psm.cInfo.ptms)
            //{
            //    if (Y.Length - ptm.index < yFrag.Length)
            //    {
            //        theoreticalMass += ptm.mass;
            //    }
            //}
            if (NeoMassCalculator.IdenticalMasses(experimentalMass, theoreticalMass, NeoFindAmbiguity.precursorMassTolerancePpm))//if match
            {
                string novelSeq = bFrag + yFrag;
                if (!psm.candidates.Any(x => x.seq.Equals(novelSeq))) //if fusion sequence was not previously assigned to this psm
                    psm.candidates.Add(new FusionCandidate(novelSeq));
            }
        }

        public static void AutoFill(NeoPsm psm)
        {
            double experimentalMass = psm.expMass;
            for (int i = 0; i < 2; i++)
            {
                string ionSequence;
                ProductType ion;
                if (i == 0)
                {
                    ionSequence = psm.nInfo.seq;
                    ion = ProductType.B;
                }
                else
                {
                    ionSequence = psm.cInfo.seq;
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
                        var ipos = Array.BinarySearch(NeoFindAmbiguity.keys, key);
                        if (ipos < 0)
                            ipos = ~ipos;

                        if (ipos > 0)
                        {
                            var downIpos = ipos - 1;
                            // Try down
                            while (downIpos >= 0)
                            {
                                closestPeak = NeoFindAmbiguity.keys[downIpos];
                                if (NeoMassCalculator.IdenticalMasses(experimentalMass, closestPeak + theoreticalMass, NeoFindAmbiguity.precursorMassTolerancePpm))
                                    foreach (string frag in NeoFindAmbiguity.massDict[closestPeak])
                                        combinations.Add(frag);
                                else
                                    break;
                                downIpos--;
                            }
                        }
                        if (ipos < NeoFindAmbiguity.keys.Length)
                        {
                            var upIpos = ipos;
                            // Try here and up
                            while (upIpos < NeoFindAmbiguity.keys.Length)
                            {
                                closestPeak = NeoFindAmbiguity.keys[upIpos];
                                if (NeoMassCalculator.IdenticalMasses(experimentalMass, closestPeak + theoreticalMass, NeoFindAmbiguity.precursorMassTolerancePpm))
                                    foreach (string frag in NeoFindAmbiguity.massDict[closestPeak])
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
                                    psm.candidates.Add(new FusionCandidate(ionFrag + s));
                                else
                                    psm.candidates.Add(new FusionCandidate(s + ionFrag));
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
            if(ionSequence.Contains("AWRRSC"))
            { }
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
                    if (NeoMassCalculator.IdenticalMasses(experimentalMass, IonMass, NeoFindAmbiguity.precursorMassTolerancePpm))
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
