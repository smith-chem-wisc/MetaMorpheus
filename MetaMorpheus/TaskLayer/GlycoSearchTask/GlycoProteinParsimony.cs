using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using EngineLayer.GlycoSearch;
using Readers;
using ThermoFisher.CommonCore.Data;

namespace TaskLayer
{
    public class GlycoProteinParsimony
    {
        //id: Accession, ProtienPos, ModId.islocalized, minQValue, maxProb

        public GlycoProteinParsimony(string proteinAccess, int proteinPos, char aminoAcid, bool isLocalized, double minQValue, LocalizationLevel bestLocalizeLevel, double maxProb)
        {
            ProteinAccession = proteinAccess;

            ProteinPos = proteinPos;

            AminoAcid = aminoAcid;

            IsLocalized = isLocalized;

            MinQValue = minQValue;

            BestLocalizeLevel = bestLocalizeLevel;

            MaxProbability = maxProb;
        }
        public string ProteinAccession { get; }

        public int ProteinPos { get; } // One-based position in the protein sequence

        public char AminoAcid { get; }

        public bool IsLocalized { get; set; }

        public double MinQValue { get; set; }

        public LocalizationLevel BestLocalizeLevel { get; set; }

        public double MaxProbability { get; set; }

        public static Dictionary<(string proteinAccession, string proteinPosition, int modId), GlycoProteinParsimony> ProteinLevelGlycoParsimony(List<GlycoSpectralMatch> allPsmsGly)
        {
            //key: proPosId
            Dictionary<(string proteinAccession, string proteinPosition, int modId), GlycoProteinParsimony> localizedMod = new Dictionary<(string proteinAccession, string proteinPosition, int modId), GlycoProteinParsimony>();

            foreach (var gsm in allPsmsGly)
            {
                if (gsm.IsContaminant || gsm.IsDecoy)
                {
                    continue;
                }

                if ((!gsm.LocalizedGlycan.IsNullOrEmpty()) && gsm.LocalizedGlycan.Count > 0)
                {
                    foreach (var local in gsm.LocalizedGlycan)
                    {
                        int proteinPos = local.SiteIndex + gsm.OneBasedStartResidue.Value - 2;

                        (string,string,int) proPosId = new (gsm.Accession, proteinPos.ToString(), local.ModId);

                        double prob = -1;
                        if (gsm.ModSitePairProbDict != null && gsm.ModSitePairProbDict.ContainsKey(local))
                        {
                            prob = local.Probability;
                        }


                        if (!localizedMod.ContainsKey(proPosId))
                        {
                            GlycoProteinParsimony gpp = new GlycoProteinParsimony(gsm.Accession, proteinPos, gsm.BaseSequence[local.SiteIndex -2], local.Confident, gsm.FdrInfo.QValue, gsm.LocalizationLevel, prob);
                            localizedMod.Add(proPosId, gpp);
                        }
                        else
                        {
                            bool islocalized = (local.Confident || localizedMod[proPosId].IsLocalized);
                            double minQValue = localizedMod[proPosId].MinQValue > gsm.FdrInfo.QValue ? gsm.FdrInfo.QValue : localizedMod[proPosId].MinQValue;
                            double maxProb = localizedMod[proPosId].MaxProbability > prob ? localizedMod[proPosId].MaxProbability : prob;
                            var localLevel = localizedMod[proPosId].BestLocalizeLevel < gsm.LocalizationLevel ? localizedMod[proPosId].BestLocalizeLevel : gsm.LocalizationLevel;

                            localizedMod[proPosId].IsLocalized = islocalized;
                            localizedMod[proPosId].MinQValue = minQValue;
                            localizedMod[proPosId].BestLocalizeLevel = localLevel;
                            localizedMod[proPosId].MaxProbability = maxProb;
                        }
                    }
                }
            }
                
            return localizedMod;
        }

    }
}
