using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.ML;
using Microsoft.ML.Data;
using Microsoft.ML.Trainers;
using System.IO;

namespace EngineLayer.GlycoSearch
{
    public class GlycanIndicator
    {
        public GlycanIndicator(bool isNGlycan, double score, bool nglycanMotifExist, double r138to144, double i_168, double i_186, double i_366, double i_274, double i_292, bool pepnn, bool pepnh, double glycanScore)
        {
            IsNGlycan = isNGlycan;
            Score = score;
            NGlycanMotifExist = nglycanMotifExist;
            R138to144 = r138to144;
            I_168 = i_168;
            I_186 = i_186;
            I_274 = i_274;
            I_366 = i_366;
            I_292 = i_292;
            PepNN = pepnn;
            PepNH = pepnh;
            GlycanScore = glycanScore;
        }

        public GlycanIndicator(GlycoSpectralMatch gsm)
        {
            IsNGlycan = gsm.NGlycan != null;
            Score = gsm.Score;
            NGlycanMotifExist = NGlycanMotifExist;
            R138to144 = gsm.R138to144;
            I_168 = gsm.OxoniumIonIntensity[7];
            I_186 = gsm.OxoniumIonIntensity[8];
            I_366 = gsm.OxoniumIonIntensity[14];
            I_274 = gsm.OxoniumIonIntensity[10];
            I_292 = gsm.OxoniumIonIntensity[12];
            PepNN = gsm.PepNN;
            PepNH = gsm.PepNH;
            GlycanScore = gsm.GlycanScore;
        }

        public bool IsNGlycan { get; set; }
        public double Score { get; set; }
        public bool NGlycanMotifExist { get; set; }

        public double R138to144 { get; set; }

        public double I_168 { get; set; }

        public double I_186 { get; set; }

        public double I_366 { get; set; }

        public double I_274 { get; set; }

        public double I_292 { get; set; }

        public bool PepNN { get; set; }

        public bool PepNH { get; set; }

        public double GlycanScore { get; set; }

        public static List<GlycanIndicator> ReadGsms(string filepath)
        {
            List<GlycanIndicator> glycanIndicators = new List<GlycanIndicator>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filepath);
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Could not read file: " + e.Message);
            }

            char[] split = new char[] { '\t' };
            Dictionary<string, int> parsedHeader = new Dictionary<string, int>();

            int lineCount = 0;
            string line;


            while (reader.Peek() > 0)
            {
                lineCount++;
                line = reader.ReadLine();
                if (lineCount == 1)
                {
                    var spls = line.Split(split);
                    foreach (var spl in spls)
                    {
                        if (spl.Length > 0)
                        {
                            parsedHeader.Add(spl, Array.IndexOf(spls, spl));
                        }
                    }

                    continue;
                }

                GlycanIndicator glycanIndicator = GenerateGlycanIndicator(line, split, parsedHeader);
                glycanIndicators.Add(glycanIndicator);
            }

            return glycanIndicators;
        }

        private static GlycanIndicator GenerateGlycanIndicator(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split);

            return new GlycanIndicator(
                (spl[parsedHeader["IsNGlycoPeptide"]] == "TRUE" ? true : false),
                double.Parse(spl[parsedHeader["Score"]]),
                (spl[parsedHeader["N-Glycan motif Check"]] == "TRUE" ? true : false),
                double.Parse(spl[parsedHeader["R138/144"]]),
                double.Parse(spl[parsedHeader["I_168"]]),
                double.Parse(spl[parsedHeader["I_186"]]),
                double.Parse(spl[parsedHeader["I_366"]]),
                double.Parse(spl[parsedHeader["I_274"]]),
                double.Parse(spl[parsedHeader["I_292"]]),
                (spl[parsedHeader["PepNN"]] == "TRUE" ? true : false),
                (spl[parsedHeader["PepNH"]] == "TRUE" ? true : false),
                double.Parse(spl[parsedHeader["NumOfContateYion"]])
                );
        }
    }

    // Example with label and 50 feature values. A data set is a collection of such examples.
    public class DataPoint
    {
        public bool Label { get; set; } //Ture if is OGlycopeptide.
        [VectorType(11)]
        public float[] Features { get; set; }
        public static IEnumerable<DataPoint> GenerateDataPoints(List<GlycanIndicator> glycanIndicators)
        {
            foreach (var gi in glycanIndicators)
            {
                yield return new DataPoint
                {
                    Label = gi.IsNGlycan,

                    Features = new float[]
                    {
                        (float)gi.Score,
                        gi.NGlycanMotifExist? 0:1,
                        (float)gi.R138to144,
                        (float)gi.I_168,
                        (float)gi.I_186,
                        (float)gi.I_366,
                        (float)gi.I_274,
                        (float)gi.I_292,
                        gi.PepNN ? 0:1,
                        gi.PepNH ? 0:1,
                        (float)gi.GlycanScore
                    }
                };

            }        
        }
    }


    // Class used to capture predictions.
    public class Prediction
    {
        // Original label.
        public bool Label { get; set; }
        // Predicted label from the trainer.
        public bool PredictedLabel { get; set; }
    }
}
