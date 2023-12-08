using EngineLayer.CrosslinkSearch;
using EngineLayer.GlycoSearch;
using EngineLayer;
using Proteomics.Fragmentation;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Serialization;

namespace TaskLayer
{
    public static class WriteGlycoFile
    {
        public static void WritePsmGlycoToTsv(List<GlycoSpectralMatch> items, string filePath, bool writeGlycoPsms)
        {
            if (items.Count == 0)
            {
                return;
            }

            string header = GlycoSpectralMatch.GetTabSepHeaderSingle();
            if (writeGlycoPsms)
            {
                header += ("\t" + GlycoSpectralMatch.GetTabSeperatedHeaderGlyco());
            }

            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(header);
                if (writeGlycoPsms)
                {
                    foreach (var heh in items)
                    {
                        output.WriteLine(heh.SingleToString() + "\t" + heh.GlycoToString());
                    }
                }
                else
                {
                    foreach (var heh in items)
                    {
                        output.WriteLine(heh.SingleToString());
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
                foreach (var item in glycoProteinParsimony.OrderBy(i => i.Key.proteinAccession))
                {
                    if (item.Value != null)
                    {
                        output.WriteLine(
                                item.Key.proteinAccession + "\t" +
                                item.Key.proteinPosition + "\t" +
                                item.Value.AminoAcid + "\t" +
                                GlycanBox.GlobalOGlycans[item.Key.glycanId].Composition + "\t" +
                                item.Value.IsLocalized + "\t" +
                                item.Value.MinQValue.ToString("0.000") + "\t" +
                                item.Value.BestLocalizeLevel + "\t" +
                                item.Value.MaxProbability.ToString("0.000"));

                    }
                }
            }
        }

        //The function is to summarize localized glycosylation of each protein site. 
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
                        String.Join(",", local.Value.Select(p => GlycanBox.GlobalOGlycans[int.Parse(p)].Composition))
                        );
                }
            }
        }
    }
}
