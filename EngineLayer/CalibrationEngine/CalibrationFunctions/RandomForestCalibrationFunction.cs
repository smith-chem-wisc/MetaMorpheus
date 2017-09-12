using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Calibration
{
    public class RandomForestCalibrationFunction : ILearner<double>
    {
        #region Private Fields

        private readonly RegressionTree[] RegressionTrees;
        #endregion Private Fields

        #region Public Constructors

        public RandomForestCalibrationFunction(int numTrees, int doNotSplitIfUnderThis)
        {
            RegressionTrees = new RegressionTree[numTrees];
            for (int i = 0; i < numTrees; i++)
                RegressionTrees[i] = new RegressionTree(doNotSplitIfUnderThis, 0);
        }

        #endregion Public Constructors

        #region Public Methods

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            Random rand = new Random(0);
            List<LabeledDataPoint> trainingListHere = new List<LabeledDataPoint>();
            for (int i = 0; i < targets.Length; i++)
                trainingListHere.Add(new LabeledDataPoint(observations.Row(i), targets[i]));
            IPredictorModel<double>[] RegressionTreesPMs = new IPredictorModel<double>[RegressionTrees.Length];
            Parallel.For(0, RegressionTrees.Length, i =>
            {
                List<LabeledDataPoint> subsampledTrainingPoints = new List<LabeledDataPoint>();
                for (int j = 0; j < trainingListHere.Count; j++)
                {
                    int index = rand.Next(trainingListHere.Count);
                    subsampledTrainingPoints.Add(trainingListHere[index]);
                }

                double[] trainList2Concat = subsampledTrainingPoints.SelectMany(b => b.Inputs).ToArray();
                F64Matrix trainList2Matrix = new F64Matrix(trainList2Concat, subsampledTrainingPoints.Count, trainList2Concat.Length / subsampledTrainingPoints.Count);
                double[] trainList2Targets = subsampledTrainingPoints.Select(b => b.Label).ToArray();

                RegressionTreesPMs[i] = RegressionTrees[i].Learn(trainList2Matrix, trainList2Targets);
            });

            return new RandomForestCalibrationPredictorModel(RegressionTreesPMs);
        }

        #endregion Public Methods
    }

    internal class RandomForestCalibrationPredictorModel : IPredictorModel<double>
    {
        #region Private Fields

        private readonly IPredictorModel<double>[] RegressionTrees;

        #endregion Private Fields

        #region Public Constructors

        public RandomForestCalibrationPredictorModel(IPredictorModel<double>[] RegressionTrees)
        {
            this.RegressionTrees = RegressionTrees;
        }

        #endregion Public Constructors

        #region Public Methods

        public double[] GetRawVariableImportance()
        {
            throw new NotImplementedException();
        }

        public Dictionary<string, double> GetVariableImportance(Dictionary<string, int> featureNameToIndex)
        {
            throw new NotImplementedException();
        }

        public double Predict(double[] observation)
        {
            return RegressionTrees.Select(b => b.Predict(observation)).Average();
        }

        #endregion Public Methods
    }
}