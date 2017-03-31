using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Threading;

namespace EngineLayer.Calibration
{
    internal class RandomForestCalibrationFunction : CalibrationFunction
    {

        #region Private Fields

        private readonly RegressionTree[] RegressionTrees;
        private readonly bool[] useFeature;
        private readonly Thread taskThread;

        #endregion Private Fields

        #region Public Constructors

        public RandomForestCalibrationFunction(int numTrees, int splitLimit, bool[] useFeature, Thread taskThread)
        {
            RegressionTrees = new RegressionTree[numTrees];
            this.useFeature = useFeature;
            for (int i = 0; i < numTrees; i++)
            {
                RegressionTrees[i] = new RegressionTree(splitLimit, 0, useFeature);
            }
            this.taskThread = taskThread;
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
            var rand = new Random();
            List<LabeledDataPoint> trainingListHere = trainingList.ToList();

            Thread currentThread = Thread.CurrentThread;
            Parallel.For(0, RegressionTrees.Length, (i, loopState) =>
            {
                if (currentThread.ThreadState == ThreadState.Aborted || taskThread != null && taskThread.ThreadState == ThreadState.Aborted)
                    loopState.Stop();

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