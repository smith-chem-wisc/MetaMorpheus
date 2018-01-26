using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Neo
{
    public class InitialID : NeoPsm
    {
        public string id { get; set; }
        public double[] peakHits { get; set; }
        public Boolean target { get; set; }
        public double qValue { get; set; }
        public string seq { get; set; }
        public string score { get; set; }

        public InitialID(int scan, double expMass, string id, string seq, string peaks, string scorem)
            : base(scan, expMass)
        {
            this.seq = CleanSeq(seq);
            this.id = id;
            this.score = score;
            try
            {
                SetPeakHits(peaks);
            }
            catch
            {
                double[] temp = new double[] { 0, 0 };
                this.peakHits = temp;
            }
        }

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
            this.peakHits = dbArray;
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
                    //ModificationOn = false;
                    //double modMass = MassCalculator.getPTMMass(ModificationName, out string e);
                    //error_message += e;
                    //PTM ptm = new PTM(ModificationName, aaIndex, modMass);
                    //this.ptms.Add(ptm);
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
