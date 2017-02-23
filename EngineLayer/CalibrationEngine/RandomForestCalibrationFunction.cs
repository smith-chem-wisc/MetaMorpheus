using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Calibration
{
    internal class RandomForestCalibrationFunction : CalibrationFunction
    {
        #region Private Fields

        private readonly RegressionTree[] RegressionTrees;
        private readonly bool[] useFeature;

        #endregion Private Fields

        #region Public Constructors

        public RandomForestCalibrationFunction(int numTrees, int splitLimit, bool[] useFeature)
        {
            RegressionTrees = new RegressionTree[numTrees];
            this.useFeature = useFeature;
            for (int i = 0; i < numTrees; i++)
                RegressionTrees[i] = new RegressionTree(splitLimit, 0, useFeature);
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return "RandomForestCalibrationFunction" + string.Join(",", useFeature);
        }

        #endregion Public Methods

        #region Internal Methods

        internal override void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
        {
            var rand = new Random(1);
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

        internal override double Predict(double[] t)
        {
            return RegressionTrees.Select(b => b.Predict(t)).Average();
        }

        #endregion Internal Methods
    }
}