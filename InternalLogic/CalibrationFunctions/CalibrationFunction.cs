using System;
using System.Collections.Generic;
using System.IO;

namespace InternalLogic
{
    public abstract class CalibrationFunction
    {
        internal string name;

        public abstract void Train(IEnumerable<LabeledDataPoint> trainingList);

        public abstract double Predict(double[] t);

        public double getMSE(IEnumerable<LabeledDataPoint> pointList)
        {
            double mse = 0;
            int count = 0;
            foreach (LabeledDataPoint p in pointList)
            {
                mse += Math.Pow(Predict(p.inputs) - p.output, 2);
                count++;
            }
            return count == 0 ? 0 : mse / count;
        }

        internal void writePredictedLables(List<LabeledDataPoint> trainList1, string v)
        {
            var fullFileName = Path.Combine(@"PredictedLabels", v + "newLabels" + ".dat");
            Directory.CreateDirectory(Path.GetDirectoryName(fullFileName));
            using (StreamWriter file = new StreamWriter(fullFileName))
            {
                file.WriteLine("PredictedLabel");
                foreach (LabeledDataPoint d in trainList1)
                    file.WriteLine(Predict(d.inputs));
            }
        }
    }
}