using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;

namespace EngineLayer.Neo
{
    public static class AggregateSearchFiles
    {
        public static void Combine(string primaryFilePath, string secondaryFilePath)
        {
            string header = "";
            string[] primaryLines = (System.IO.File.ReadAllLines(primaryFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(secondaryFilePath));
            List<PsmTsvLine> aggregatedLines = AggregateTargetDecoy(primaryLines, secondaryLines);
            CalculateFDR(aggregatedLines);

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"\\bison\share\users\Zach\NeoPaper\e001318\Parameter_Benchmarking\500\2017-12-18-13-52-20\Task1-Neo_Cis\180116NormalNeoCis_Aggregated_MinusOnePointFive.txt"))
            {
                file.WriteLine(header);
                foreach (PsmTsvLine psm in aggregatedLines)
                    file.WriteLine(psm.ToString());
            }
        }

        private static List<PsmTsvLine> AggregateTargetDecoy(string[] primaryLines, string[] secondaryLines)
        {
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);
            List<PsmTsvLine> aggregatedLines = new List<PsmTsvLine>();
            int p = 0;
            int s = 0;
            while (p < primaryPsms.Count || s < secondaryPsms.Count)
            {
                PsmTsvLine psmP = primaryPsms[p];
                PsmTsvLine psmS = secondaryPsms[s];
                if (psmP.scanNumber < psmS.scanNumber)
                {
                    aggregatedLines.Add(psmP);
                    p++;
                }
                else if (psmP.scanNumber > psmS.scanNumber)
                {
                    aggregatedLines.Add(psmS);
                    s++;
                }
                else
                {
                    if (psmP.score > psmS.score - 0.001 && psmP.score < psmS.score + 0.001)
                    {
                        aggregatedLines.Add(psmP.AggregateLine(psmS));
                        p++;
                        s++;
                    }
                    else if (psmP.score > psmS.score)
                    {
                        aggregatedLines.Add(psmP);
                        p++;
                    }
                    else
                    {
                        aggregatedLines.Add(psmS);
                        s++;
                    }
                }
            }
            return aggregatedLines;
        }

        private static List<NeoPsm> 

        private static void CalculateFDR(List<PsmTsvLine> aggregatedLines)
        {
            aggregatedLines = aggregatedLines.OrderByDescending(x => x.score).ToList();
            int targets = 0;
            int decoys = 0;
            foreach (PsmTsvLine line in aggregatedLines)
            {
                if (line.DCT.Contains("T") || line.DCT.Contains("C"))
                    targets++;
                else
                    decoys++;

                line.target = targets.ToString();
                line.decoy = decoys.ToString();
                line.q = ((1.0d * decoys) / (targets)).ToString();
            }
        }
    }
}
