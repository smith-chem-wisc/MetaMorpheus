using System;
using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class InitialID : NeoPsm
    {
        public InitialID(int scan, double expMass, string id, string seq, string peaks, string scorem)
            : base(scan, expMass)
        {
            Seq = CleanSeq(seq);
            ID = id;
            Score = Score;
            try
            {
                SetPeakHits(peaks);
            }
            catch
            {
                double[] temp = new double[] { 0, 0 };
                PeakHits = temp;
            }
        }

        public string ID { get; set; }
        public double[] PeakHits { get; set; }
        public Boolean Target { get; set; }
        public double QValue { get; set; }
        public string Seq { get; set; }
        public string Score { get; set; }

        public void SetPeakHits(string peakHits)
        {
            peakHits = peakHits.Replace("];[", ","); //used for EThcD to convert partition into addition
            peakHits = peakHits.Replace('"', '[');
            peakHits = peakHits.Replace("[", "");
            peakHits = peakHits.Replace("]", "");
            peakHits = peakHits.Replace(" ", "");
            peakHits = peakHits.Replace(";", "");
            string[] strArray = peakHits.Split(','); //split as strings
            List<double> dbList = new List<double>(); //double list to sort multi-ion input
            foreach (string str in strArray)
            {
                dbList.Add(Convert.ToDouble(str));
            }
            dbList.Sort();
            double[] dbArray = new double[strArray.Length]; //new double array
            for (int i = 0; i < strArray.Length; i++)
            {
                dbArray[i] = dbList[i];
            }
            PeakHits = dbArray;
        }

        //removes all modifications (in parenthesis) and returns the base sequence
        public string CleanSeq(string seq)
        {
            bool ModificationOn = false;
            string ModificationName = "";
            string cleanedSeq = "";
            string[] seqArray = seq.Split(' ');
            foreach (char amino_acid in seqArray[0]) //if there are synonymous peaks, then the sequences must be identical or possess ambiguities that will be caught later
            {
                if (amino_acid == ']') //only occurs at end of mod
                {
                    ModificationOn = false;
                }
                if (ModificationOn == true) //only occurs if "(" already found
                {
                    ModificationName += amino_acid;
                }
                if (amino_acid == '[') //start collecting PTM name
                {
                    ModificationOn = true;
                }
                if (ModificationOn == false && amino_acid != ']')
                {
                    cleanedSeq += amino_acid;
                }
            }
            return cleanedSeq;
        }
    }
}