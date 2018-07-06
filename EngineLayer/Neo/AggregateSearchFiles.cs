using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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

            using (StreamWriter file = new StreamWriter(@outputFilePath + "_TargetsAndDecoys.psmtsv"))
            {
                file.WriteLine(header);
                foreach (PsmTsvLine psm in aggregatedLines)
                    file.WriteLine(psm.ToString());
            }

            List<PsmTsvLine> targets = AssignFDRToTarget(primaryLines, secondaryLines);
            using (StreamWriter file = new StreamWriter(@outputFilePath + "_Targets.psmtsv"))
            {
                file.WriteLine(header);
                foreach (PsmTsvLine psm in targets)
                    file.WriteLine(psm.ToString());
            }
        }

        public static void RecursiveNeoAggregation(string standardFilePath, string neoResultFilePath, string outputFolder, string identifier)
        {
            //This method determines the optimum cutoff for gold standard identification and the minimum score difference required for a splice to outscore a normal
            double qThreshold = 0;
            double oldQThreshold = 0;
            double scoreDifferenceThreshold = 0;
            double oldScoreDifferenceThreshold = 0;
            int numSplicedHighScoreQ = -1; //highest number of splice assignments at a 1% local FDR
            int numSplicedHighScoreThreshold = -1; //highest number of splice assignments at a 1% local FDR
            int numSplicedScore = 0; //current number of splic assignments at a 1% local FDR
            bool increaseQ = true;
            bool increaseScoreDifference = true;

            string[] primaryLines = (System.IO.File.ReadAllLines(@standardFilePath));
            string[] secondaryLines = (System.IO.File.ReadAllLines(@neoResultFilePath));
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            primaryPsms.ForEach(x => x.NeoType = NeoType.Normal);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);

            do //determine if score difference should be changed
            {
                if (numSplicedScore <= numSplicedHighScoreThreshold) //check second time around if move score Threhold the other way
                {
                    increaseScoreDifference = false;
                }
                else
                {
                    numSplicedHighScoreThreshold = numSplicedScore; //update highscore
                    oldScoreDifferenceThreshold = scoreDifferenceThreshold; //update score difference
                }
                scoreDifferenceThreshold = UpdateScoreDifferenceThreshold(scoreDifferenceThreshold, increaseScoreDifference); //update score Threshold

                do //determine gold standards to use
                {
                    if (numSplicedScore <= numSplicedHighScoreQ) //check second time around if move qValue the other way
                    {
                        increaseQ = false;
                    }
                    else
                    {
                        numSplicedHighScoreQ = numSplicedScore; //update highscore
                        oldQThreshold = qThreshold; //updateQ
                    }

                    qThreshold = UpdateQThreshold(primaryPsms, qThreshold, increaseQ); //get qValue
                    List<PsmTsvLine> aggregatedLines = Percolate(primaryPsms, secondaryPsms, qThreshold, scoreDifferenceThreshold);
                    numSplicedScore = CalculateNumberOfConfidentSpliced(aggregatedLines);
                } while (numSplicedScore > numSplicedHighScoreQ || increaseQ); //do again the otherway if done increasing
                List<PsmTsvLine> oldAggregatedLines = Percolate(primaryPsms, secondaryPsms, oldQThreshold, scoreDifferenceThreshold);
                numSplicedScore = CalculateNumberOfConfidentSpliced(oldAggregatedLines);
                increaseQ = true;
                qThreshold = oldQThreshold;
            } while (numSplicedScore > numSplicedHighScoreThreshold || increaseScoreDifference);
            List<PsmTsvLine> finalAggregatedLines = Percolate(primaryPsms, secondaryPsms, oldQThreshold, oldScoreDifferenceThreshold);
            using (StreamWriter file = new StreamWriter(Path.Combine(outputFolder, identifier)))
            {
                file.WriteLine(primaryLines[0]); //header
                foreach (PsmTsvLine line in finalAggregatedLines)
                {
                    file.WriteLine(line.ToString());
                }
            }
            using (StreamWriter file = new StreamWriter(Path.Combine(outputFolder, "PercolatorInfo_" + identifier)))
            {
                file.WriteLine("Maxmimum q-Value of Gold Standards: " + oldQThreshold);
                file.WriteLine("Minimum Score Difference for Splice Selection Over Normal: " + oldScoreDifferenceThreshold);
            }
        }

        public static double UpdateQThreshold(List<PsmTsvLine> primaryLines, double qThreshold, bool increaseQ)
        {
            List<double> qValues = primaryLines.Select(x => Convert.ToDouble(x.Q)).ToList(); //grab all q values
            if (increaseQ) //get next highest qValue
            {
                qValues = qValues.OrderBy(x => x).ToList();
                foreach (double qValue in qValues)
                {
                    if (qValue > qThreshold)
                    {
                        return qValue;
                    }
                }
            }
            else //get next lowest qValue
            {
                qValues = qValues.OrderByDescending(x => x).ToList();
                foreach (double qValue in qValues)
                {
                    if (qValue < qThreshold)
                    {
                        return qValue;
                    }
                }
            }
            return qThreshold; //if nothing, return the same qValue
        }

        public static double UpdateScoreDifferenceThreshold(double scoreDifferenceThreshold, bool increaseScoreDifference)
        {
            return increaseScoreDifference ? scoreDifferenceThreshold += 0.01 : scoreDifferenceThreshold -= 0.01;
        }

        public static List<PsmTsvLine> Percolate(List<PsmTsvLine> primaryPsms, List<PsmTsvLine> secondaryPsms, double qThreshold, double scoreDifferenceThreshold)
        {
            List<PsmTsvLine> aggregatedLines = new List<PsmTsvLine>();
            int p = 0;
            int s = 0;
            //While loop generates a combined list of primary (normal) and secondary (spliced) psms based on scoreDifferenceThreshold
            while (p < primaryPsms.Count && s < secondaryPsms.Count)
            {
                PsmTsvLine psmP = primaryPsms[p];
                PsmTsvLine psmS = secondaryPsms[s];
                if (psmP.ScanNumber < psmS.ScanNumber)
                {
                    aggregatedLines.Add(psmP);
                    p++;
                }
                else if (psmP.ScanNumber > psmS.ScanNumber)
                {
                    aggregatedLines.Add(psmS);
                    s++;
                }
                else
                {
                    if (psmP.Score > psmS.Score - scoreDifferenceThreshold)
                    {
                        aggregatedLines.Add(psmP);
                    }
                    else
                    {
                        psmS.NeoType = (Convert.ToDouble(psmP.Q) < qThreshold) ? NeoType.DecoySpliced : NeoType.Spliced; //if the beaten score belonged to a gold standard, this is a false discovery
                        aggregatedLines.Add(psmS);
                    }
                    p++;
                    s++;
                }
            }

            // Wrap up any leftover psms without scanCounts
            for (; p < primaryPsms.Count; p++)
            {
                aggregatedLines.Add(primaryPsms[p]);
            }
            for (; s < secondaryPsms.Count; s++)
            {
                aggregatedLines.Add(secondaryPsms[s]);
            }

            return aggregatedLines;
        }

        public static int CalculateNumberOfConfidentSpliced(List<PsmTsvLine> aggregatedLines)
        {
            int numConfidentSpliced = 1; //prevent divide-by-zero error, subtract later
            int numDecoySpliced = 0;
            aggregatedLines = aggregatedLines.OrderByDescending(x => x.Score).ToList();
            foreach (PsmTsvLine line in aggregatedLines)
            {
                if ((1.0 * numDecoySpliced) / numConfidentSpliced >= 0.05)
                {
                    break;
                }
                if (line.NeoType.Equals(NeoType.Spliced))
                {
                    numConfidentSpliced++;
                }
                else if (line.NeoType.Equals(NeoType.DecoySpliced))
                {
                    numDecoySpliced++;
                }
            }
            return numConfidentSpliced;
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
                if (psmP.ScanNumber < psmS.ScanNumber)
                {
                    aggregatedLines.Add(psmP);
                    p++;
                }
                else if (psmP.ScanNumber > psmS.ScanNumber)
                {
                    aggregatedLines.Add(psmS);
                    s++;
                }
                else
                {
                    if (psmP.Score > psmS.Score - 0.001 && psmP.Score < psmS.Score + 0.001)
                    {
                        aggregatedLines.Add(psmP.AggregateLine(psmS));
                    }
                    else if (psmP.Score > psmS.Score)
                    {
                        aggregatedLines.Add(psmP);
                    }
                    else
                    {
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
            return aggregatedLines;
        }

        private static void CalculateFDR(List<PsmTsvLine> aggregatedLines)
        {
            aggregatedLines = aggregatedLines.OrderByDescending(x => x.Score).ToList();
            int targets = 0;
            int decoys = 0;
            foreach (PsmTsvLine line in aggregatedLines)
            {
                if (line.DCT.Contains("T") || line.DCT.Contains("C"))
                {
                    targets++;
                }
                else
                {
                    decoys++;
                }

                line.Target = targets.ToString();
                line.Decoy = decoys.ToString();
                line.Q = ((1.0d * decoys) / (targets)).ToString();
            }
        }

        private static List<PsmTsvLine> AssignFDRToTarget(string[] primaryLines, string[] secondaryLines)
        {
            List<PsmTsvLine> primaryPsms = ImportPsmtsv.ImportLinesToAggregate(primaryLines);
            List<PsmTsvLine> secondaryPsms = ImportPsmtsv.ImportLinesToAggregate(secondaryLines);
            primaryPsms = primaryPsms.OrderByDescending(x => x.Score).ToList();
            secondaryPsms = secondaryPsms.OrderByDescending(x => x.Score).ToList();
            int p = 0;
            int s = 0;
            int target = 0;
            int decoy = 0;
            double qMax = 0;
            PsmTsvLine decoyLine = secondaryPsms[s];

            while (p < primaryPsms.Count)
            {
                PsmTsvLine targetLine = primaryPsms[p];
                if (targetLine.Score > decoyLine.Score || s == secondaryPsms.Count)
                {
                    target++;
                    targetLine.Target = target.ToString();
                    targetLine.Decoy = decoy.ToString();
                    double qValue = (1.0d * decoy / target);
                    qMax = (qMax > qValue) ? qMax : qValue;
                    targetLine.Q = qMax.ToString();
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
    }
}