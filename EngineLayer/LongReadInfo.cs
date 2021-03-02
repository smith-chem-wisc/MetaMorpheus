using Proteomics;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class LongReadInfo
    {
        // protein inference weight for proteins w/ no corresponding ORF calling data.
        // this would be used if UniProt proteins are added for genes that weren't observed in the transcriptomics data
        public static readonly double ProteinInferenceWeightForNoTranscriptomicsData = 0.5;

        public readonly string ProteinAccession;
        public readonly double CPM;
        public readonly double ProteinInferenceWeight;

        public LongReadInfo(string protein, double cpm)
        {
            this.ProteinAccession = protein;
            this.CPM = cpm;
            ProteinInferenceWeight = CalculateProteinInferenceWeight();
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(CPM.ToString("F2"));

            return sb.ToString();
        }

        private double CalculateProteinInferenceWeight()
        {
            if (CPM < 0.5)
            {
                return 0.1;
            }
            else if (CPM < 1)
            {
                return 0.3;
            }
            else if (CPM < 5)
            {
                return 0.5;
            }
            else if (CPM < 10)
            {
                return 0.7;
            }

            return 1;
        }
    }
}
