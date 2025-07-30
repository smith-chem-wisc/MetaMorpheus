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
        //id: Accession, ProtienPos, GlycanId.islocalized, minQValue, maxProb

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

        public static Dictionary<(string proteinAccession, string proteinPosition, int glycanId), GlycoProteinParsimony> ProteinLevelGlycoParsimony(List<GlycoSpectralMatch> allPsmsGly)
        {
            //key: proPosId
            Dictionary<(string proteinAccession, string proteinPosition, int glycanId), GlycoProteinParsimony> localizedGlycan = new Dictionary<(string proteinAccession, string proteinPosition, int glycanId), GlycoProteinParsimony>();

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
                        if (gsm.SiteSpeciLocalProb != null && gsm.SiteSpeciLocalProb.ContainsKey(local.SiteIndex))
                        {
                            prob = gsm.SiteSpeciLocalProb[local.SiteIndex].Where(p => p.Item1 == local.ModId).FirstOrDefault().Item2;
                        }


                        if (!localizedGlycan.ContainsKey(proPosId))
                        {
                            GlycoProteinParsimony gpp = new GlycoProteinParsimony(gsm.Accession, proteinPos, gsm.BaseSequence[local.SiteIndex -2], local.Confident, gsm.FdrInfo.QValue, gsm.LocalizationLevel, prob);
                            localizedGlycan.Add(proPosId, gpp);
                        }
                        else
                        {
                            bool islocalized = (local.Confident || localizedGlycan[proPosId].IsLocalized);
                            double minQValue = localizedGlycan[proPosId].MinQValue > gsm.FdrInfo.QValue ? gsm.FdrInfo.QValue : localizedGlycan[proPosId].MinQValue;
                            double maxProb = localizedGlycan[proPosId].MaxProbability > prob ? localizedGlycan[proPosId].MaxProbability : prob;
                            var localLevel = localizedGlycan[proPosId].BestLocalizeLevel < gsm.LocalizationLevel ? localizedGlycan[proPosId].BestLocalizeLevel : gsm.LocalizationLevel;

                            localizedGlycan[proPosId].IsLocalized = islocalized;
                            localizedGlycan[proPosId].MinQValue = minQValue;
                            localizedGlycan[proPosId].BestLocalizeLevel = localLevel;
                            localizedGlycan[proPosId].MaxProbability = maxProb;
                        }
                    }
                }
            }

            return localizedGlycan;
        }

    }
}
