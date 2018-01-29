using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;

namespace EngineLayer.Neo
{
    public static class AggregateSearchFiles
    {
        public static void Combine(string primaryFilePath, string secondaryFilePath, string outputFilePath)
        {
            string header = "";
            string[] primaryLines = (System.IO.File.ReadAllLines(@primaryFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(@secondaryFilePath));
            List<PsmTsvLine> aggregatedLines = AggregateDifferentDatabaseSearches(primaryLines, secondaryLines);
            CalculateFDR(aggregatedLines);

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@outputFilePath+"_TargetsAndDecoys.psmtsv"))
            {
                file.WriteLine(header);
                foreach (PsmTsvLine psm in aggregatedLines)
                    file.WriteLine(psm.ToString());
            }

            List<PsmTsvLine> targets = AssignFDRToTarget(primaryLines, secondaryLines);
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@outputFilePath + "_Targets.psmtsv"))
            {
                file.WriteLine(header);
                foreach (PsmTsvLine psm in targets)
                    file.WriteLine(psm.ToString());
            }
        }

        private static List<PsmTsvLine> AggregateDifferentDatabaseSearches(string[] primaryLines, string[] secondaryLines)
        {
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);
            List<PsmTsvLine> aggregatedLines = new List<PsmTsvLine>();
            int p = 0;
            int s = 0;
            while (p < primaryPsms.Count && s < secondaryPsms.Count)
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
                        aggregatedLines.Add(psmP.AggregateLine(psmS));
                    else if (psmP.score > psmS.score)
                        aggregatedLines.Add(psmP);
                    else                  
                        aggregatedLines.Add(psmS);
                    p++;
                    s++;
                }
            }
            while (p < primaryPsms.Count)
            {
                aggregatedLines.Add(primaryPsms[p]);
                p++;
            }
            while (s < secondaryPsms.Count)
            {
                aggregatedLines.Add(primaryPsms[s]);
                s++;
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

        private static List<PsmTsvLine> AssignFDRToTarget(string[] primaryLines, string[] secondaryLines)
        {
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);
            primaryPsms = primaryPsms.OrderByDescending(x => x.score).ToList();
            secondaryPsms = secondaryPsms.OrderByDescending(x => x.score).ToList();
            int p = 0;
            int s = 0;
            int target = 0;
            int decoy = 0;
            double qMax = 0;
            PsmTsvLine decoyLine = secondaryPsms[s];

            while (p < primaryPsms.Count)
            {
                PsmTsvLine targetLine = primaryPsms[p];
                if (targetLine.score > decoyLine.score || s == secondaryPsms.Count)
                {
                    target++;
                    targetLine.target = target.ToString();
                    targetLine.decoy = decoy.ToString();
                    double qValue = (1.0d * decoy / target);
                    qMax = (qMax > qValue) ? qMax : qValue;
                    targetLine.q = qMax.ToString();
                    p++;
                }
                else
                {
                    decoy++;
                    s++;
                }
            }
            return primaryPsms;
        }

        public static void RecursiveNeoAggregation(string standardFilePath, string neoResultFilePath)
        {
            double qThreshold = 0;
            double oldQThreshold = 0;
            double scoreDifferenceThreshold = 0;
            double oldScoreDifferenceThreshold = 0;
            int numSplicedHighScore = -1;
            int numSplicedScore = 0;
            bool increaseQ = true;
            bool increaseScoreDifference = true;

            string[] primaryLines = (System.IO.File.ReadAllLines(@standardFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(@neoResultFilePath));
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            primaryPsms.ForEach(x => x.neoType = PsmTsvLine.NeoType.Normal);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);
            List<PsmTsvLine> aggregatedLines = new List<PsmTsvLine>();
            do
            {
                if (!(numSplicedScore > numSplicedHighScore))
                    increaseScoreDifference = false;
                do
                {
                    if (!(numSplicedScore > numSplicedHighScore))
                        increaseQ = false;
                    else
                    {
                        numSplicedHighScore = numSplicedScore;
                        oldQThreshold = qThreshold;
                        oldScoreDifferenceThreshold = scoreDifferenceThreshold;
                    }

                    qThreshold = UpdateQThreshold(primaryPsms, qThreshold, increaseQ);
                    int p = 0;
                    int s = 0;
                    while (p < primaryPsms.Count && s < secondaryPsms.Count)
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
                            if (psmP.score > psmS.score - scoreDifferenceThreshold)
                            {
                                aggregatedLines.Add(psmP);
                            }
                            else
                            {
                                psmS.neoType = (Convert.ToDouble(psmP.q) < qThreshold) ? PsmTsvLine.NeoType.DecoySpliced : PsmTsvLine.NeoType.Spliced;
                                aggregatedLines.Add(psmS);
                            }
                            p++;
                            s++;
                        }
                    }
                    while (p < primaryPsms.Count)
                    {
                        aggregatedLines.Add(primaryPsms[p]);
                        p++;
                    }
                    while (s < secondaryPsms.Count)
                    {
                        aggregatedLines.Add(primaryPsms[s]);
                        s++;
                    }
                } while (numSplicedScore > numSplicedHighScore || increaseQ);
                scoreDifferenceThreshold = UpdateScoreDifferenceThreshold(scoreDifferenceThreshold, increaseScoreDifference);
            } while (numSplicedScore > numSplicedHighScore||increaseScoreDifference);
        }


        public static double UpdateQThreshold(List<PsmTsvLine> primaryLines, double qThreshold, bool increaseQ)
        {
            List<double> qValues = primaryLines.Select(x => Convert.ToDouble(x.q)).ToList();
            if (increaseQ)
            {
                qValues = qValues.OrderBy(x => x).ToList();
                foreach (double qValue in qValues)
                    if (qValue > qThreshold)
                        return qValue;
            }
            else
            {
                qValues = qValues.OrderByDescending(x => x).ToList();
                foreach (double qValue in qValues)
                    if (qValue < qThreshold)
                        return qValue;
            }
            return qThreshold;
        }

        public static double UpdateScoreDifferenceThreshold(double scoreDifferenceThreshold, bool increaseScoreDifference)
        {
            if (increaseScoreDifference)
                return scoreDifferenceThreshold += 0.01;
            else
                return scoreDifferenceThreshold -= 0.01;
        }
    }
}
