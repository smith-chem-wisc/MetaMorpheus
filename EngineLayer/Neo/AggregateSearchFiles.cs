using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;

namespace EngineLayer.Neo
{
    public static class AggregateSearchFiles
    {

        public static int scanNumberIndex = -1;
        public static int scoreIndex = -1;
        public static int baseIndex = -1;
        public static int fullIndex = -1;
        public static int accessionIndex = -1;
        public static int proteinIndex = -1;
        public static int geneIndex = -1;
        public static int DCTIndex = -1;
        public static int targetIndex = -1;
        public static int decoyIndex = -1;
        public static int qIndex = -1;

        public static readonly string scanNumberHeader = "Scan Number";
        public static readonly string scoreHeader = "Score";
        public static readonly string baseHeader = "Base Sequence";
        public static readonly string fullHeader = "Full Sequence";
        public static readonly string accessionHeader = "Protein Accession";
        public static readonly string proteinHeader = "Protein Name";
        public static readonly string geneHeader = "Gene Name";
        public static readonly string DCTHeader = "Decoy / Contaminant / Target";
        public static readonly string targetHeader = "Cumulative Target";
        public static readonly string decoyHeader = "Cumulative Decoy";
        public static readonly string qHeader = "QValue";

        public static void Combine(string primaryFilePath, string secondaryFilePath)
        {
            string header = "";
            string[] primaryLines = (System.IO.File.ReadAllLines(primaryFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(secondaryFilePath));
            List<PsmTsvLine> aggregatedLines = ParsePsmTsv(primaryLines, secondaryLines);
            CalculateFDR(aggregatedLines);

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"\\bison\share\users\Zach\NeoPaper\e001318\Parameter_Benchmarking\500\2017-12-18-13-52-20\Task1-Neo_Cis\180116NormalNeoCis_Aggregated_MinusOnePointFive.txt"))
            {
                file.WriteLine(header);
                foreach (PsmTsvLine psm in aggregatedLines)
                    file.WriteLine(psm.ToString());
            }
        }

        private static List<PsmTsvLine> ParsePsmTsv(string[] primaryLines, string[] secondaryLines)
        {
            List<PsmTsvLine> primaryPsms = new List<PsmTsvLine>();
            List<PsmTsvLine> secondaryPsms = new List<PsmTsvLine>();

            string[] headerArray = primaryLines[0].Split('\t');
            for (int i = 0; i < headerArray.Length; i++)
            {
                string currentHeader = headerArray[i];
                if (currentHeader.Equals(scanNumberHeader))
                    scanNumberIndex = i;
                else if (currentHeader.Equals(scoreHeader))
                    scoreIndex = i;
                else if (currentHeader.Equals(baseHeader))
                    baseIndex = i;
                else if (currentHeader.Equals(fullHeader))
                    fullIndex = i;
                else if (currentHeader.Equals(accessionHeader))
                    accessionIndex = i;
                else if (currentHeader.Equals(proteinHeader))
                    proteinIndex = i;
                else if (currentHeader.Equals(geneHeader))
                    geneIndex = i;
                else if (currentHeader.Equals(DCTHeader))
                    DCTIndex = i;

            }
            for (int i = 1; i < primaryLines.Length; i++)
            {
                string[] lineArray = primaryLines[i].Split('\t');
                primaryPsms.Add(new PsmTsvLine(lineArray, Convert.ToInt32(lineArray[scanNumberIndex]), Convert.ToDouble(lineArray[scoreIndex]), lineArray[baseIndex], lineArray[fullIndex], lineArray[accessionIndex], lineArray[proteinIndex], lineArray[geneIndex], lineArray[DCTIndex], lineArray[targetIndex], lineArray[decoyIndex], lineArray[qIndex]));
            }
            for (int i = 1; i < secondaryLines.Length; i++)
            {
                string[] lineArray = secondaryLines[i].Split('\t');
                secondaryPsms.Add(new PsmTsvLine(lineArray, Convert.ToInt32(lineArray[scanNumberIndex]), Convert.ToDouble(lineArray[scoreIndex]), lineArray[baseIndex], lineArray[fullIndex], lineArray[accessionIndex], lineArray[proteinIndex], lineArray[geneIndex], lineArray[DCTIndex], lineArray[targetIndex], lineArray[decoyIndex], lineArray[qIndex]));
            }
            primaryPsms = primaryPsms.OrderBy(x => x.scanNumber).ToList();
            secondaryPsms = secondaryPsms.OrderBy(x => x.scanNumber).ToList();

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
