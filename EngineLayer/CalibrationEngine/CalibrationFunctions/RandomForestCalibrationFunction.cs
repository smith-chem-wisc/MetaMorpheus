using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Calibration
{
    public class RandomForestCalibrationFunction : CalibrationFunction
    {
        #region Private Fields

        private readonly RegressionTree[] RegressionTrees;
        private readonly bool[] useFeature;
        private readonly Random rand;

        #endregion Private Fields

        #region Public Constructors

        public RandomForestCalibrationFunction(int numTrees, int doNotSplitIfUnderThis, bool[] useFeature, Random rand)
        {
            RegressionTrees = new RegressionTree[numTrees];
            this.useFeature = useFeature;
            this.rand = rand;
            for (int i = 0; i < numTrees; i++)
                RegressionTrees[i] = new RegressionTree(doNotSplitIfUnderThis, 0, useFeature);
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return "RandomForestCalibrationFunction" + string.Join(",", useFeature);
        }

        public override void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
        {
            List<LabeledDataPoint> trainingListHere = trainingList.ToList();
            Parallel.For(0, RegressionTrees.Length, i =>
            {
                List<LabeledDataPoint> subsampledTrainingPoints = new List<LabeledDataPoint>();
                for (int j = 0; j < trainingListHere.Count; j++)
                {
                    int index = rand.Next(trainingListHere.Count);
                    subsampledTrainingPoints.Add(trainingListHere[index]);
                }
                RegressionTrees[i].Train(subsampledTrainingPoints);
            });
        }

        public override double Predict(double[] t)
        {
            return RegressionTrees.Select(b => b.Predict(t)).Average();
        }

        #endregion Public Methods
    }
}