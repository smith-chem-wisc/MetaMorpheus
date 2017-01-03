using IndexSearchAndAnalyze;
using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace FragmentGeneration
{
    public class Writing
    {
        public static void WriteTree(BinTreeStructure myTreeStructure, string output_folder, string fileName)
        {
            using (StreamWriter output = new StreamWriter(Path.Combine(output_folder, fileName + ".mytsv")))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tNtermLocFrac\tCtermLocFrac\tUniprot");
                foreach (Bin bin in myTreeStructure.finalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.NlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }
        }

        public static void WriteToTabDelimitedTextFileWithDecoys(List<NewPsmWithFDR> items, string output_folder, string fileName)
        {
            Console.WriteLine("Writing psms");
            using (StreamWriter output = new StreamWriter(Path.Combine(output_folder, fileName + ".psmtsv")))
            {
                output.WriteLine("Spectrum File\tScan Number\tRetention Time\tPrecursor MZ\tPrecursor Charge\tPrecursor Intensity\tExperimental Peaks\tTotal Intensity\tPrecursor Mass\tScoreFromSearch\tPreviousAminoAcid\tSequence\tNextAminoAcid\tnumVariableMods\tStart Residue\tEnd Residue\tPeptide\tMissed Cleavages\tPeptide Mass\tProtein\tMass Diff(Da)\tMatched Fragments\tMatched Counts\tLocalized Scores\tImprovement\tImprovment Residue\tImprovement Terminus\tDecoy\tCumulative Target\tCumulative Decoy\tQ-value");
                for (int i = 0; i < items.Count; i++)
                    output.WriteLine(items[i].ToString());
            }
        }
    }
}