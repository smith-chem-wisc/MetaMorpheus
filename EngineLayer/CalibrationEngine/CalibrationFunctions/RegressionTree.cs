using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Calibration
{
    internal class RegressionTree : CalibrationFunction
    {

        #region Private Fields

        private readonly int level;
        private readonly bool[] useFeature;
        private RegressionTree leftTree;
        private RegressionTree rightTree;

        private double output = double.NaN;

        private double bestValue;
        private int bestI = -1;
        private int doNotSplitIfUnderThis;

        #endregion Private Fields

        #region Public Constructors

        public RegressionTree(int doNotSplitIfUnderThis, int level, bool[] useFeature)
        {
            this.doNotSplitIfUnderThis = doNotSplitIfUnderThis;
            this.level = level;
            this.useFeature = useFeature;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal override double Predict(double[] t)
        {
            if (double.IsNaN(output))
            {
                if (t[bestI] < bestValue)
                    return leftTree.Predict(t);
                else
                    return rightTree.Predict(t);
            }
            else
                return output;
        }

        internal override void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
        {
            var trainingPoints = trainingList.ToList();
            var averageOutputs = trainingPoints.Select(b => b.label).Average();
            if (trainingPoints.Count() < doNotSplitIfUnderThis)
            {
                output = averageOutputs;
                return;
            }
            var bestSumSquaredErrors = trainingPoints.Select(b => Math.Pow(averageOutputs - b.label, 2)).Sum();

            var prunedTrainingPoints = trainingPoints;

            // For every variable, try to find the best split
            for (int i = 0; i < useFeature.Length; i++)
                if (useFeature[i])
                {
                    prunedTrainingPoints.Sort(Comparer<LabeledDataPoint>.Create((x, y) => x.inputs[i].CompareTo(y.inputs[i])));
                    int num_splits = Math.Min(15, prunedTrainingPoints.Count - 1);
                    for (double j = 0; j < num_splits; j++)
                    {
                        double quantile = prunedTrainingPoints.Select(b => b.inputs[i]).Quantile((j + 1) / (num_splits + 1));
                        if (double.IsNaN(quantile))
                            break;
                        if (quantile == prunedTrainingPoints.First().inputs[i] || quantile == prunedTrainingPoints.Last().inputs[i])
                            continue;
                        double averageFirst = prunedTrainingPoints.TakeWhile(b => b.inputs[i] < quantile).Select(b => b.label).Average();
                        double averageLast = prunedTrainingPoints.SkipWhile(b => b.inputs[i] < quantile).Select(b => b.label).Average();
                        var sumSquaredErrors = prunedTrainingPoints.TakeWhile(b => b.inputs[i] < quantile).Select(b => Math.Pow(averageFirst - b.label, 2)).Sum() +
                                               prunedTrainingPoints.SkipWhile(b => b.inputs[i] < quantile).Select(b => Math.Pow(averageLast - b.label, 2)).Sum();
                        if (sumSquaredErrors < bestSumSquaredErrors)
                        {
                            bestSumSquaredErrors = sumSquaredErrors;
                            bestValue = quantile;
                            bestI = i;
                        }
                    }
                }

            bool reachedBottom = false;
            if (bestI == -1)
                reachedBottom = true;

            if (reachedBottom)
                output = trainingPoints.Select(b => b.label).Average();
            else
            {
                trainingPoints.Sort(Comparer<LabeledDataPoint>.Create((x, y) => x.inputs[bestI].CompareTo(y.inputs[bestI])));
                leftTree = new RegressionTree(doNotSplitIfUnderThis, level + 1, useFeature);
                leftTree.Train(trainingPoints.TakeWhile(b => b.inputs[bestI] < bestValue).ToList());
                rightTree = new RegressionTree(doNotSplitIfUnderThis, level + 1, useFeature);
                rightTree.Train(trainingPoints.SkipWhile(b => b.inputs[bestI] < bestValue).ToList());
            }
        }

        #endregion Internal Methods

    }
}