using EngineLayer.GlycoSearch;
using EngineLayer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public static class WriteGlycoFile
    {
        public static void WritePsmGlycoToTsv(List<GlycoSpectralMatch> gsms, string filePath, bool writeGlycoPsms)
        {
            if (gsms.Count == 0)
            {
                return;
            }

            string header = GlycoSpectralMatch.GetTabSepHeaderSingle();
            if (writeGlycoPsms)
            {
                header += (GlycoSpectralMatch.GetTabSeperatedHeaderGlyco());
            }

            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(header);
                if (writeGlycoPsms)
                {
                    foreach (var gsm in gsms)
                    {
                        output.WriteLine(gsm.SingleToString() + "\t" + gsm.GlycoToString());
                    }
                }
                else
                {
                    foreach (var gsm in gsms)
                    {
                        output.WriteLine(gsm.SingleToString());
                    }
                }

            }
        }

        //The function is to summarize localized glycan by protein site.
        public static void WriteSeenProteinGlycoLocalization(Dictionary<(string proteinAccession, string proteinPosition, int glycanId), GlycoProteinParsimony> glycoProteinParsimony, string outputPath)
        {
            if (glycoProteinParsimony.Keys.Count == 0)
            { return; }
            var writtenFile = Path.Combine(outputPath);
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("Protein Accession\tModification Site\tAminoAcid\tLocalized Glycans\tLocalized\tLowest Qvalue\tBest Localization Level\tMax Site Specific Probability");
                foreach (var protein in glycoProteinParsimony.OrderBy(i => i.Key.proteinAccession))
                {
                    if (protein.Value != null)
                    {
                        output.WriteLine(
                                protein.Key.proteinAccession + "\t" +
                                protein.Key.proteinPosition + "\t" +
                                protein.Value.AminoAcid + "\t" + 
                                (protein.Key.glycanId >= 0 ? GlycanBox.GlobalOGlycans[protein.Key.glycanId].Composition:
                                GlycanBox.GlobalNGlycans[protein.Key.glycanId].Composition) + "\t" +
                                protein.Value.IsLocalized + "\t" +
                                protein.Value.MinQValue.ToString("0.000") + "\t" +
                                protein.Value.BestLocalizeLevel + "\t" +
                                protein.Value.MaxProbability.ToString("0.000"));

                    }
                }
            }
        }

        /// <summary>
        /// To summarize localized glycosylation of each protein site. The filter parameter is MinQValue <= 0.01 and IsLocalized = true.
        /// </summary>
        /// <param name="glycoProteinParsimony"></param>
        /// <param name="outputPath"></param>
        public static void WriteProteinGlycoLocalization(Dictionary<(string proteinAccession, string proteinPosition, int glycanId), GlycoProteinParsimony> glycoProteinParsimony, string outputPath)
        {
            if (glycoProteinParsimony.Count == 0)
            { return; }

            Dictionary<string, HashSet<string>> localizedglycans = new Dictionary<string, HashSet<string>>();
            foreach (var item in glycoProteinParsimony.Where(p => p.Value.IsLocalized && p.Value.MinQValue <= 0.01))
            {

                var key = item.Key.proteinAccession + "#" + item.Key.proteinPosition;
                if (localizedglycans.ContainsKey(key))
                {
                    localizedglycans[key].Add(item.Key.glycanId.ToString());
                }
                else
                {
                    localizedglycans[key] = new HashSet<string>();
                    localizedglycans[key].Add(item.Key.glycanId.ToString());
                }

            }

            var writtenFile = Path.Combine(outputPath);
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("Protein Accession\tModification Site\tLocalized Glycan Number\tLocalized Glycans");
                foreach (var local in localizedglycans.OrderBy(p => p.Key))
                {
                    var x = local.Key.Split('#');
                    output.WriteLine(
                        x[0] + "\t" +
                        x[1] + "\t" +
                        local.Value.Count() + "\t" +
                        String.Join(",", local.Value.Select(p => int.Parse(p) >= 0? GlycanBox.GlobalOGlycans[int.Parse(p)].Composition 
                            : GlycanBox.GlobalNGlycans[int.Parse(p)].Composition))
                        );
                }
            }
        }
    }
}
