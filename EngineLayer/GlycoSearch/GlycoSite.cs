using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using Proteomics.Fragmentation;

namespace EngineLayer.GlycoSearch
{
    public class GlycoSite
    { 
        public GlycoSite(int modSite, int glycoID, bool isLocalized)
        {
            ModSite = modSite;
            GlycanID = glycoID;
            IsLocalized = isLocalized;
        }

        //mod site, glycanID, isLocalized
        public int ModSite { get; set; }

        public int GlycanID { get; set; }

        public bool IsLocalized { get; set; }

        public static List<GlycoSite> GetLocalizedGlycan(List<Route> routes, out LocalizationLevel localizationLevel)
        {
            List<GlycoSite> localizedGlycan = new List<GlycoSite>();

            //Dictionary<string, int>: modsite-id, count
            Dictionary<string, int> seenModSite = new Dictionary<string, int>();

            foreach (var route in routes)
            {
                foreach (var gs in route.Mods)
                {
                    var k = gs.ModSite.ToString() + "-" + gs.GlycanID.ToString();
                    if (seenModSite.ContainsKey(k))
                    {
                        seenModSite[k] += 1;
                    }
                    else
                    {
                        seenModSite.Add(k, 1);
                    }
                }
            }

            localizationLevel = LocalizationLevel.Level3;
            if (routes.Count == 1)
            {
                localizationLevel = LocalizationLevel.Level1;
            }
            else if (routes.Count > 1)
            {
                //When there is one glyco site appear in all routes, then the glyco site is localized.
                if (seenModSite.Values.Where(p => p == routes.Count).Count() > 0)
                {
                    localizationLevel = LocalizationLevel.Level2;
                }
                else
                {
                    localizationLevel = LocalizationLevel.Level3;
                }
            }

            foreach (var seenMod in seenModSite)
            {
                if (seenMod.Value == routes.Count)
                {
                    localizedGlycan.Add(new GlycoSite(int.Parse(seenMod.Key.Split('-')[0]), int.Parse(seenMod.Key.Split('-')[1]), true));
                }
                else
                {
                    localizedGlycan.Add(new GlycoSite(int.Parse(seenMod.Key.Split('-')[0]), int.Parse(seenMod.Key.Split('-')[1]), false));
                }
            }

            return localizedGlycan;
        }

        //Glyco Localization
        public static void GlycoLocalizationCalculation(List<GlycoSpectralMatch> gsms, CommonParameters CommonParameters)
        {
            foreach (var glycoSpectralMatch in gsms)
            {
                if (glycoSpectralMatch.LocalizationGraphs == null)
                {
                    continue;
                }

                if (glycoSpectralMatch.LocalizationGraphs != null)
                {
                    //TO DO: The is_HCD_only_data is not totally correct. 
                    if (glycoSpectralMatch.GlycanType == GlycoType.OGlycoPep &&
                        !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType) && 
                        !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType))
                    {
                        glycoSpectralMatch.LocalizationLevel = LocalizationLevel.Level3;
                        if (glycoSpectralMatch.LocalizationGraphs.Count == 1 && glycoSpectralMatch.LocalizationGraphs.First().ModPos.Length == 1)
                        {
                            glycoSpectralMatch.LocalizationLevel = LocalizationLevel.Level1b;
                        }

                    }
                    else
                    {
                        List<Route> localizationCandidates = new List<Route>();

                        for (int i = 0; i < glycoSpectralMatch.LocalizationGraphs.Count; i++)
                        {
                            var allPathWithMaxScore = LocalizationGraph.GetAllHighestScorePaths(glycoSpectralMatch.LocalizationGraphs[i].array, glycoSpectralMatch.LocalizationGraphs[i].ChildModBoxes);

                            foreach (var path in allPathWithMaxScore)
                            {
                                var local = LocalizationGraph.GetLocalizedPath(glycoSpectralMatch.LocalizationGraphs[i], path);
                                local.ModBoxId = glycoSpectralMatch.LocalizationGraphs[i].ModBoxId;
                                localizationCandidates.Add(local);
                            }
                        }

                        glycoSpectralMatch.Routes = localizationCandidates;
                    }
                }

                if (glycoSpectralMatch.Routes != null)
                {
                    LocalizationLevel localLevel;
                    glycoSpectralMatch.LocalizedGlycan = GlycoSite.GetLocalizedGlycan(glycoSpectralMatch.Routes, out localLevel);
                    glycoSpectralMatch.LocalizationLevel = localLevel;

                    //Localization PValue.
                    if (localLevel == LocalizationLevel.Level1 || localLevel == LocalizationLevel.Level2)
                    {
                        List<Route> allRoutes = new List<Route>();
                        foreach (var graph in glycoSpectralMatch.LocalizationGraphs)
                        {
                            allRoutes.AddRange(LocalizationGraph.GetAllPaths_CalP(graph, glycoSpectralMatch.ScanInfo_p, glycoSpectralMatch.Thero_n));
                        }
                        glycoSpectralMatch.SiteSpeciLocalProb = LocalizationGraph.CalSiteSpecificLocalizationProbability(allRoutes, glycoSpectralMatch.LocalizationGraphs.First().ModPos);
                    }
                }

                CorrectLocalizationLevel(glycoSpectralMatch);

            }
        }

        //Correct Localization Level based on site specific probability. If LocalizationLevel = 1, and there are site probability lower than 0.75, Correct the level to 1b.
        public static void CorrectLocalizationLevel(GlycoSpectralMatch gsm)
        {
            if (gsm.SiteSpeciLocalProb == null || gsm.LocalizationLevel != LocalizationLevel.Level1)
            {
                return;
            }

            if (gsm.LocalizationGraphs.First().ModPos.Length == 1 && gsm.LocalizationGraphs.First().TotalScore == 0)
            {
                gsm.LocalizationLevel = LocalizationLevel.Level1b;
                return;
            }


            for (int i = 0; i < gsm.LocalizedGlycan.Count; i++)
            {
                var g = gsm.LocalizedGlycan[i];
                if (gsm.SiteSpeciLocalProb[g.ModSite].Where(p => p.Item1 == g.GlycanID).First().Item2 < 0.75)
                {
                    gsm.LocalizationLevel = LocalizationLevel.Level1b;
                    return;
                }

                if (!gsm.Routes.First().Mods[i].IsLocalized)
                {
                    gsm.LocalizationLevel = LocalizationLevel.Level1b;
                    return;
                }
            }
        }

        //Deprecated function. 
        //For N-Glycopeptide localzation Level.
        public static void NGlycoCorrectLocalizationLevel(List<GlycoSpectralMatch> ngsms)
        {
            foreach (var gsm in ngsms)
            {
                //Correct nglycoSpectrumMatch localizationLevel to Level1b based on matched fragment ions.
                if (gsm.ModPos.Count() == 1)
                {
                    gsm.LocalizationLevel = LocalizationLevel.Level1b;
                    //correct localization.
                    if (gsm.MatchedFragmentIons.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.b || p.NeutralTheoreticalProduct.ProductType == ProductType.y) && p.NeutralTheoreticalProduct.NeutralLoss > 0).Count() > 0
                        || gsm.MatchedFragmentIons.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.c && p.NeutralTheoreticalProduct.AminoAcidPosition >= gsm.ModPos[0] - 1) || (p.NeutralTheoreticalProduct.ProductType == ProductType.zDot && p.NeutralTheoreticalProduct.AminoAcidPosition <= gsm.ModPos[0] - 1)).Count() > 0)
                    {
                        gsm.LocalizationLevel = LocalizationLevel.Level1;
                    }
                    else
                    {
                        foreach (var childIons in gsm.ChildMatchedFragmentIons)
                        {
                            if (childIons.Value.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.b || p.NeutralTheoreticalProduct.ProductType == ProductType.y) && p.NeutralTheoreticalProduct.NeutralLoss > 0).Count() > 0
                                || childIons.Value.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.c && p.NeutralTheoreticalProduct.AminoAcidPosition >= gsm.ModPos[0] - 1) || (p.NeutralTheoreticalProduct.ProductType == ProductType.zDot && p.NeutralTheoreticalProduct.AminoAcidPosition <= gsm.ModPos[0] - 1)).Count() > 0)
                            {
                                gsm.LocalizationLevel = LocalizationLevel.Level1;
                            }
                        }
                    }
                }
            }
        }
    }
}
