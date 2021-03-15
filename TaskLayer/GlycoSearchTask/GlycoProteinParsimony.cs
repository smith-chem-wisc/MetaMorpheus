using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using EngineLayer;
using EngineLayer.GlycoSearch;

namespace TaskLayer
{
    public class GlycoProteinParsimony
    {
        //id: ProteinAccession, ProtienPos, GlycanId.islocalized, minQValue, maxProb

        public GlycoProteinParsimony(string proteinAccess, int proteinPos, char aminoAcid, string glycanCom, bool isLocalized, double minQValue, LocalizationLevel bestLocalizeLevel, double maxProb)
        {
            ProteinAccession = proteinAccess;

            ProteinPos = proteinPos;

            AminoAcid = aminoAcid;

            GlycanCom = glycanCom;

            IsLocalized = isLocalized;

            MinQValue = minQValue;

            BestLocalizeLevel = bestLocalizeLevel;

            MaxProbability = maxProb;
        }
        public string ProteinAccession { get; }

        public int ProteinPos { get; }

        public char AminoAcid { get; }

        public string GlycanCom { get; }

        public bool IsLocalized { get; set; }

        public double MinQValue { get; set; }

        public LocalizationLevel BestLocalizeLevel { get; set; }

        public double MaxProbability { get; set; }

        public static Dictionary<string, GlycoProteinParsimony> ProteinLevelGlycoParsimony(List<GlycoSpectralMatch> allPsmsGly)
        {
            //key: proPosId
            Dictionary<string, GlycoProteinParsimony> localizedGlycan = new Dictionary<string, GlycoProteinParsimony>();

            foreach (var gsm in allPsmsGly)
            {
                if (gsm.IsContaminant || gsm.IsDecoy)
                {
                    continue;
                }

                if (gsm.LocalizedGlycan!=null && gsm.LocalizedGlycan.Count > 0)
                {
                    foreach (var local in gsm.LocalizedGlycan)
                    {
                        string glycanCom = "";
                        if (gsm.GlycanType == GlycoType.OGlycoPep)
                        {
                            glycanCom = GlycanBox.GlobalOGlycans[local.GlycanID].Composition;
                        }
                        else if (gsm.GlycanType == GlycoType.NGlycoPep)
                        {
                            glycanCom = GlycanBox.GlobalNGlycans[local.GlycanID].Composition;
                        }
                        else
                        {
                            glycanCom = GlycanBox.GlobalMixedGlycans[local.GlycanID].Composition;
                        }

                        int proteinPos = local.ModSite + gsm.OneBasedStartResidueInProtein.Value - 2;

                        string proPosId = gsm.ProteinAccession + "-" + proteinPos.ToString() + "-" + glycanCom;

                        double prob = -1;
                        if (gsm.SiteSpeciLocalProb != null && gsm.SiteSpeciLocalProb.ContainsKey(local.ModSite))
                        {
                            prob = gsm.SiteSpeciLocalProb[local.ModSite].Where(p => p.Item1 == local.GlycanID).FirstOrDefault().Item2;
                        }


                        if (!localizedGlycan.ContainsKey(proPosId))
                        {
                            GlycoProteinParsimony gpp = new GlycoProteinParsimony(gsm.ProteinAccession, proteinPos, gsm.BaseSequence[local.ModSite - 2], glycanCom, local.IsLocalized, gsm.FdrInfo.QValue, gsm.LocalizationLevel, prob);
                            localizedGlycan.Add(proPosId, gpp);
                        }
                        else
                        {
                            bool islocalized = (local.IsLocalized || localizedGlycan[proPosId].IsLocalized);
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
