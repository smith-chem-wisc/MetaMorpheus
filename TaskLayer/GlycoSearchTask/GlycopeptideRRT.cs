using EngineLayer;
using EngineLayer.GlycoSearch;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics.LinearRegression;
using FlashLFQ;
using MathNet.Numerics.LinearAlgebra;
using SharpLearning.Containers.Arithmetic;

namespace TaskLayer
{
    public class GlycopeptideRRT
    {
        public GlycopeptideRRT(GlycoSpectralMatch[] gsms)
        {
            GSMs = gsms;
        }

        public GlycoSpectralMatch[] GSMs { get; }

        public byte[] Kind
        {
            get
            {
                return GlycoSpectralMatch.GetKind(GSMs.First());
            }
        }

        public string BaseSeq
        {
            get
            {
                return GSMs.First().BaseSequence;
            }
        }

        public string FullSeq
        {
            get
            {
                return GSMs.First().FullSequence;
            }
        }


        public int SameBaseSeqCount {get; set;}

        public double PredictedRetentionTime { get; set; }

        public double OriginalRetentionTime { get; set; }

        public double Intensity { get; set; }

        public double[] _x { get; set; }


        //Retention time prediction.
        public static List<GlycopeptideRRT> ConstructGlycopeptideRRT(List<GlycoSpectralMatch> allGsms, FlashLfqResults flashLfqResults)
        {
            var peaks = flashLfqResults.Peaks.SelectMany(p => p.Value).ToList();

            List<GlycopeptideRRT> glycopeptideRRTs = new List<GlycopeptideRRT>();

            foreach (var gsms in allGsms.GroupBy(p=>p.FullSequence).Select(p =>p.ToArray()))
            {
                GlycopeptideRRT glycopeptideRRT = new GlycopeptideRRT(gsms);

                glycopeptideRRTs.Add(glycopeptideRRT);
            }

            foreach (var gp in glycopeptideRRTs)
            {
                //peaks are from FlashLFQ result. They contains the Corrected retention time and intensity.
                var apex = peaks.Where(p => p.Identifications.Select(x => x.ModifiedSequence).Contains(gp.FullSeq)).First().Apex;
                if (apex != null)
                {
                    gp.OriginalRetentionTime = apex.IndexedPeak.RetentionTime;
                    gp.Intensity = apex.Intensity;
                }
                else
                {
                    gp.OriginalRetentionTime = gp.GSMs.First().ScanRetentionTime;
                }
         
            }

            return glycopeptideRRTs;
        }

        //Apply Weighted multiple linear regression for relative retention time prediction of glycopeptide.
        //The method is from "Retention Time Prediction for Glycopeptides in Reversed-Phase Chromatography for Glycoproteomic Applications"
        //To understand the code, please read the paper first.
        //Basicly, the method create matrix and used MathNet.MultipleRegression for calculation.
        public static string PredictRRT(List<GlycopeptideRRT> glycopeptideRRTs)
        {           
            //kind is an arrary contains all the sugar. To build the matrix, we first need to figure out which sugar has been identified in the glycopeptides.
            bool[] ExistSugar = new bool[glycopeptideRRTs.First().Kind.Length];
            List<int> indexOfExistSugar = new List<int>();
            for (int i = 0; i < ExistSugar.Length; i++)
            {
                var exist = glycopeptideRRTs.Sum(p => p.Kind[i]) > 0;
                ExistSugar[i] = exist;
            }
            for (int i = 0; i < ExistSugar.Length; i++)
            {
                if (ExistSugar[i])
                {
                    indexOfExistSugar.Add(i);
                }
            }

            //Group the glycopeptide based on BaseSequence.
            var gps_byId = glycopeptideRRTs.GroupBy(p => p.BaseSeq).Where(p => p.Count() >= 2);

            List<double[]> x = new List<double[]>();

            List<double> y = new List<double>();
            List<double> intensity = new List<double>();
            int ind = 0;
            //here we decide only the baseSequence with more than 3 glycans will be used for RRT.The size of the matrix is decided of baseSequence count. 
            int peptideCount = gps_byId.Count();
            while (ind < peptideCount)
            {
                var gps = gps_byId.ElementAt(ind);

                foreach (var gp in gps)
                {
                    double[] _x = new double[peptideCount + indexOfExistSugar.Count()];                   
                    _x[ind] = 1;
                    //var kind = GlycoSpectralMatch.GetKind(gp.GSMs.First());
                    for (int i = 0; i < indexOfExistSugar.Count(); i++)
                    {
                        _x[peptideCount + i] = gp.Kind[indexOfExistSugar[i]];                
                    }
                    x.Add(_x);
                    gp._x = _x;
                    y.Add(gp.OriginalRetentionTime);
                    intensity.Add(gp.Intensity);
                }

                ind++;
            }

            //double[] _p = MultipleRegression.QR(x.ToArray(), y.ToArray(), intercept: false);

            //var w = new[] { new[] { 1, 0, 0 }, new[] { 0, 1, 0 }, new[] { 0, 0, 1 } };
            double[][] w = new double[x.Count][];
            var maxIntensity = intensity.Max();
            for (int i = 0; i < x.Count(); i++)
            {
                w[i] = new double[x.Count];

                double intense = 0.01;
                if (intensity[i] > 0)
                {
                    intense = intensity[i] / maxIntensity;
                }
                w[i][i] = intense;
            }

            Matrix<double> xm = Matrix<double>.Build.Dense(x.Count(), x.First().Length, (i, j) => x[i][j]);
            Matrix<double> wm = Matrix<double>.Build.Dense(x.Count(), w.First().Length, (i, j) => w[i][j]);
            Vector<double> yv = Vector<double>.Build.Dense(y.Count(), (i) => y[i]);

            string predictMessage = "Successfully Predicted Relative Retention Time.";

            try
            {
                var p = WeightedRegression.Weighted(xm, yv, wm);

                foreach (var gp in glycopeptideRRTs)
                {
                    double z = -1.0;
                    if (gp._x != null)
                    {
                        for (int i = 0; i < gp._x.Length; i++)
                        {
                            z += p[i] * gp._x[i];
                        }
                    }
                    gp.PredictedRetentionTime = z;
                }
                
            }
            catch (Exception)
            {
                predictMessage = "Failed to Predict Relative Retention Time.";
                //throw;
            }

            return predictMessage;

        }

        public static void WriteGlycopeptideRRT(List<GlycopeptideRRT> glycopeptideRRTs, string filepath)
        {
            using (StreamWriter output = new StreamWriter(filepath))
            {
                //for (int i = 0; i < x.Count(); i++)
                foreach (var gp in glycopeptideRRTs)
                {         
                    if (gp._x != null)
                    {
                        string line = "";
                        for (int j = 0; j < gp._x.Length; j++)
                        {
                            line += gp._x[j] + "\t";
                        }
                        line += gp.OriginalRetentionTime + "\t";
                        line += gp.PredictedRetentionTime;
                        output.WriteLine(line);
                    }
                }  
            }
        }

        public static double[] RetentionTimeRegression(List<GlycoSpectralMatch> glycoSpectralMatches)
        {
            double[][] xs = new double[glycoSpectralMatches.Count][];

            double[] ys = glycoSpectralMatches.Select(p => p.ScanRetentionTime).ToArray();

            for (int i = 0; i < glycoSpectralMatches.Count; i++)
            {
                var glycanKind = GlycanBox.OGlycanBoxes[glycoSpectralMatches[i].Routes.First().ModBoxId].Kind;

                xs[i] = new double[1 + glycanKind.Length];

                xs[i][0] = glycoSpectralMatches[i].PredictedHydrophobicity;

                for (int j = 1; j <= glycanKind.Length; j++)
                {
                    xs[i][j] = glycanKind[j - 1];
                }
            }

            using (StreamWriter output = new StreamWriter(@"E:\MassData\Glycan\Nick_2019_StcE\Rep1\_temp2\2020-05-19-17-12-20\Task1-GlycoSearchTask\RetentionTime.csv"))
            {
                for (int i = 0; i < xs.Length; i++)
                {
                    string line = "";
                    for (int j = 0; j < xs[0].Length; j++)
                    {
                        line += xs[i][j].ToString() + "\t";
                    }
                    line += ys[i].ToString();
                    output.WriteLine(line);
                }
            }

            double[] p = MultipleRegression.QR(xs, ys, intercept: false);

            //double[] p = Fit.MultiDim(xs, ys, intercept: true);

            return p;
        }

        public static void RetentionTimePrediction(List<GlycoSpectralMatch> glycoSpectralMatches, double[] regression)
        {
            foreach (var gsm in glycoSpectralMatches)
            {
                var glycanKind = GlycanBox.OGlycanBoxes[gsm.Routes.First().ModBoxId].Kind;

                var xs = new double[1 + glycanKind.Length];

                xs[0] = gsm.PredictedHydrophobicity;

                for (int j = 1; j <= glycanKind.Length; j++)
                {
                    xs[j] = glycanKind[j - 1];
                }

                double prt = regression[0];

                for (int j = 1; j <= xs.Length; j++)
                {
                    prt += regression[j] * xs[j - 1];
                }

                gsm.PredictedRT = prt;
            }
        }
    }
}
