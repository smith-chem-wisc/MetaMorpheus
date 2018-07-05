﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.Neo
{
    public static class ImportPsmtsv
    {
        public static readonly string scanNumberHeader = "Scan Number";
        public static readonly string scoreHeader = "Score";
        public static readonly string baseHeader = "Base Sequence";
        public static readonly string fullHeader = "Full Sequence";
        public static readonly string accessionHeader = "Protein Accession";
        public static readonly string proteinHeader = "Protein Name";
        public static readonly string geneHeader = "Gene Name";
        public static readonly string DCTHeader = "Decoy/Contaminant/Target";
        public static readonly string targetHeader = "Cumulative Target";
        public static readonly string decoyHeader = "Cumulative Decoy";
        public static readonly string qHeader = "QValue";
        public static readonly string scanPrecursorMassHeader = "Precursor Mass";
        public static readonly string matchedIonsHeader = "Matched Ion Masses";
        public static readonly string matchedionCountsHeader = "Matched Ion Counts";
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
        public static int scanPrecursorMassIndex = -1;
        public static int matchedIonsIndex = -1;
        public static int matchedIonCountsIndex = -1;

        public static void ParseHeader(string header)
        {
            string[] headerArray = header.Split('\t');
            for (int i = 0; i < headerArray.Length; i++)
            {
                string currentHeader = headerArray[i];
                if (currentHeader.Equals(scanNumberHeader))
                    scanNumberIndex = i;
                else if (currentHeader.Equals(scoreHeader))
                    scoreIndex = i;
                else if (currentHeader.Equals(baseHeader))
                {
                    baseIndex = i;
                    fullIndex = i; //workaround for open mass searches generating thousands of combinations; eventually patch
                }
                else if (currentHeader.Equals(accessionHeader))
                    accessionIndex = i;
                else if (currentHeader.Equals(proteinHeader))
                    proteinIndex = i;
                else if (currentHeader.Equals(geneHeader))
                    geneIndex = i;
                else if (currentHeader.Equals(DCTHeader))
                    DCTIndex = i;
                else if (currentHeader.Equals(targetHeader))
                    targetIndex = i;
                else if (currentHeader.Equals(decoyHeader))
                    decoyIndex = i;
                else if (currentHeader.Equals(qHeader))
                    qIndex = i;
                else if (currentHeader.Equals(scanPrecursorMassHeader))
                    scanPrecursorMassIndex = i;
                else if (currentHeader.Equals(matchedIonsHeader))
                    matchedIonsIndex = i;
                else if (currentHeader.Equals(matchedionCountsHeader))
                    matchedIonCountsIndex = i;
            }
        }

        public static List<PsmTsvLine> ImportLinesToAggregate(string[] lines)
        {
            List<PsmTsvLine> results = new List<PsmTsvLine>();

            ParseHeader(lines[0]);
            for (int i = 1; i < lines.Length; i++)
            {
                string[] lineArray = lines[i].Split('\t');
                results.Add(new PsmTsvLine(lineArray, Convert.ToInt32(lineArray[scanNumberIndex]), Convert.ToDouble(lineArray[scoreIndex]), lineArray[baseIndex], lineArray[fullIndex], lineArray[accessionIndex], lineArray[proteinIndex], lineArray[geneIndex], lineArray[DCTIndex], lineArray[targetIndex], lineArray[decoyIndex], lineArray[qIndex]));
            }
            return results.OrderBy(x => x.scanNumber).ToList();
        }

        public static List<NeoPsm> ImportNeoPsms(string nFileName, string cFileName)
        {
            string[] nInput = File.ReadAllLines(nFileName);
            string[] cInput = File.ReadAllLines(cFileName);
            List<NeoPsm> psms = new List<NeoPsm>();
            ParseHeader(nInput[0]);

            List<InitialID> nAssignment = new List<InitialID>();
            List<InitialID> cAssignment = new List<InitialID>();

            for (int i = 1; i < nInput.Count(); i++)
            {
                string[] line = nInput[i].Split('\t').ToArray();
                InitialID id = new InitialID(Convert.ToInt32(line[scanNumberIndex]), Convert.ToDouble(line[scanPrecursorMassIndex]), line[accessionIndex], line[fullIndex], line[matchedIonsIndex], line[scoreIndex]);
                nAssignment.Add(id);
            }

            for (int i = 1; i < cInput.Count(); i++)
            {
                string[] line = cInput[i].Split('\t').ToArray();
                InitialID id = new InitialID(Convert.ToInt32(line[scanNumberIndex]), Convert.ToDouble(line[scanPrecursorMassIndex]), line[accessionIndex], line[fullIndex], line[matchedIonsIndex], line[scoreIndex]);
                cAssignment.Add(id);
            }
            //sort by scan number
            List<InitialID> nAssignmentSorted = nAssignment.OrderBy(o => o.scanNumber).ToList();
            List<InitialID> cAssignmentSorted = cAssignment.OrderBy(o => o.scanNumber).ToList();

            //remove scans not found in both files
            double maxCount = nAssignmentSorted.Count();
            for (int i = 0; i < maxCount; i++)
            {
                if (i < cAssignmentSorted.Count())
                {
                    if (nAssignmentSorted[i].scanNumber.Equals(cAssignmentSorted[i].scanNumber))
                    {
                        NeoPsm psm = new NeoPsm(nAssignmentSorted[i].scanNumber, nAssignmentSorted[i].expMass, nAssignmentSorted[i], cAssignmentSorted[i]);
                        psms.Add(psm);
                        continue;
                    }
                    else if (nAssignmentSorted[i].scanNumber < cAssignmentSorted[i].scanNumber) //no information was found for the b scan using y ions, so remove it
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