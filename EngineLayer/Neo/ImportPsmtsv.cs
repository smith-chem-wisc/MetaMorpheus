using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.Neo
{
    public static class ImportPsmtsv
    {
        public static readonly string ScanNumberHeader = "Scan Number";
        public static readonly string ScoreHeader = "Score";
        public static readonly string BaseHeader = "Base Sequence";
        public static readonly string FullHeader = "Full Sequence";
        public static readonly string AccessionHeader = "Protein Accession";
        public static readonly string ProteinHeader = "Protein Name";
        public static readonly string GeneHeader = "Gene Name";
        public static readonly string DCTHeader = "Decoy/Contaminant/Target";
        public static readonly string TargetHeader = "Cumulative Target";
        public static readonly string DecoyHeader = "Cumulative Decoy";
        public static readonly string QHeader = "QValue";
        public static readonly string ScanPrecursorMassHeader = "Precursor Mass";
        public static readonly string MatchedIonsHeader = "Matched Ion Masses";
        public static readonly string MatchedionCountsHeader = "Matched Ion Counts";
        public static int ScanNumberIndex = -1;
        public static int ScoreIndex = -1;
        public static int BaseIndex = -1;
        public static int FullIndex = -1;
        public static int AccessionIndex = -1;
        public static int ProteinIndex = -1;
        public static int GeneIndex = -1;
        public static int DCTIndex = -1;
        public static int TargetIndex = -1;
        public static int DecoyIndex = -1;
        public static int QIndex = -1;
        public static int ScanPrecursorMassIndex = -1;
        public static int MatchedIonsIndex = -1;
        public static int MatchedIonCountsIndex = -1;

        public static void ParseHeader(string header)
        {
            string[] headerArray = header.Split('\t');
            for (int i = 0; i < headerArray.Length; i++)
            {
                string currentHeader = headerArray[i];
                if (currentHeader.Equals(ScanNumberHeader))
                {
                    ScanNumberIndex = i;
                }
                else if (currentHeader.Equals(ScoreHeader))
                {
                    ScoreIndex = i;
                }
                else if (currentHeader.Equals(BaseHeader))
                {
                    BaseIndex = i;
                    FullIndex = i; //workaround for open mass searches generating thousands of combinations; eventually patch
                }
                else if (currentHeader.Equals(AccessionHeader))
                {
                    AccessionIndex = i;
                }
                else if (currentHeader.Equals(ProteinHeader))
                {
                    ProteinIndex = i;
                }
                else if (currentHeader.Equals(GeneHeader))
                {
                    GeneIndex = i;
                }
                else if (currentHeader.Equals(DCTHeader))
                {
                    DCTIndex = i;
                }
                else if (currentHeader.Equals(TargetHeader))
                {
                    TargetIndex = i;
                }
                else if (currentHeader.Equals(DecoyHeader))
                {
                    DecoyIndex = i;
                }
                else if (currentHeader.Equals(QHeader))
                {
                    QIndex = i;
                }
                else if (currentHeader.Equals(ScanPrecursorMassHeader))
                {
                    ScanPrecursorMassIndex = i;
                }
                else if (currentHeader.Equals(MatchedIonsHeader))
                {
                    MatchedIonsIndex = i;
                }
                else if (currentHeader.Equals(MatchedionCountsHeader))
                {
                    MatchedIonCountsIndex = i;
                }
            }
        }

        public static List<PsmTsvLine> ImportLinesToAggregate(string[] lines)
        {
            List<PsmTsvLine> results = new List<PsmTsvLine>();

            ParseHeader(lines[0]);
            for (int i = 1; i < lines.Length; i++)
            {
                string[] lineArray = lines[i].Split('\t');
                results.Add(new PsmTsvLine(lineArray, Convert.ToInt32(lineArray[ScanNumberIndex]), Convert.ToDouble(lineArray[ScoreIndex]), lineArray[BaseIndex], lineArray[FullIndex], lineArray[AccessionIndex], lineArray[ProteinIndex], lineArray[GeneIndex], lineArray[DCTIndex], lineArray[TargetIndex], lineArray[DecoyIndex], lineArray[QIndex]));
            }
            return results.OrderBy(x => x.ScanNumber).ToList();
        }

        public static List<NeoPsm> ImportNeoPsms(string nFileName, string cFileName)
        {
            string[] nInput = File.ReadAllLines(nFileName);
            string[] cInput = File.ReadAllLines(cFileName);
            List<NeoPsm> psms = new List<NeoPsm>();
            ParseHeader(nInput[0]);

            List<InitialID> nAssignment = new List<InitialID>();
            List<InitialID> cAssignment = new List<InitialID>();

            for (int i = 1; i < nInput.Length; i++)
            {
                string[] line = nInput[i].Split('\t').ToArray();
                InitialID id = new InitialID(Convert.ToInt32(line[ScanNumberIndex]), Convert.ToDouble(line[ScanPrecursorMassIndex]), line[AccessionIndex], line[FullIndex], line[MatchedIonsIndex], line[ScoreIndex]);
                nAssignment.Add(id);
            }

            for (int i = 1; i < cInput.Length; i++)
            {
                string[] line = cInput[i].Split('\t').ToArray();
                InitialID id = new InitialID(Convert.ToInt32(line[ScanNumberIndex]), Convert.ToDouble(line[ScanPrecursorMassIndex]), line[AccessionIndex], line[FullIndex], line[MatchedIonsIndex], line[ScoreIndex]);
                cAssignment.Add(id);
            }
            //sort by scan number
            List<InitialID> nAssignmentSorted = nAssignment.OrderBy(o => o.ScanNumber).ToList();
            List<InitialID> cAssignmentSorted = cAssignment.OrderBy(o => o.ScanNumber).ToList();

            //remove scans not found in both files
            double maxCount = nAssignmentSorted.Count;
            for (int i = 0; i < maxCount; i++)
            {
                if (i < cAssignmentSorted.Count)
                {
                    if (nAssignmentSorted[i].ScanNumber.Equals(cAssignmentSorted[i].ScanNumber))
                    {
                        NeoPsm psm = new NeoPsm(nAssignmentSorted[i].ScanNumber, nAssignmentSorted[i].ExpMass, nAssignmentSorted[i], cAssignmentSorted[i]);
                        psms.Add(psm);
                        continue;
                    }
                    else if (nAssignmentSorted[i].ScanNumber < cAssignmentSorted[i].ScanNumber) //no information was found for the b scan using y ions, so remove it
                    {
                        nAssignmentSorted.Remove(nAssignmentSorted[i]);
                        maxCount--;
                        i--;
                    }
                    else  //no information was found for the y scan using b ions, so remove it
                    {
                        cAssignmentSorted.Remove(cAssignmentSorted[i]);
                        i--;
                    }
                }
            }

            return psms;
        }
    }
}