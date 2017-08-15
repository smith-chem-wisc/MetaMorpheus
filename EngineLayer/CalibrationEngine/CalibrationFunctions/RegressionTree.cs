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
        private readonly int doNotSplitIfUnderThis;
        private RegressionTree leftTree;
        private RegressionTree rightTree;

        private double output = double.NaN;

        private double bestValue;
        private int bestI = -1;

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
            var averageOutputs = trainingPoints.Select(b => b.Label).Average();
            if (trainingPoints.Count < doNotSplitIfUnderThis)
            {
                output = averageOutputs;
                return;
            }
            var bestSumSquaredErrors = trainingPoints.Select(b => Math.Pow(averageOutputs - b.Label, 2)).Sum();

            //if (level == 0)
            //{
            //    Console.WriteLine("useFeature = " + string.Join(",", useFeature));
            //    Console.WriteLine("averageOutputs = " + averageOutputs);
            //    Console.WriteLine("trainingPoints.Count = " + trainingPoints.Count);
            //    Console.WriteLine("bestSumSquaredErrors = " + bestSumSquaredErrors);
            //}
            var prunedTrainingPoints = trainingPoints;

            // For every variable, try to find the best split
            for (int i = 0; i < useFeature.Length; i++)
                if (useFeature[i])
                {
                    //if (level == 0)
                    //    Console.WriteLine(" i = " + i);
                    prunedTrainingPoints.Sort(Comparer<LabeledDataPoint>.Create((x, y) => x.Inputs[i].CompareTo(y.Inputs[i])));
                    int num_splits = Math.Min(15, prunedTrainingPoints.Count - 1);
                    for (double j = 0; j < num_splits; j++)
                    {
                        double quantile = prunedTrainingPoints.Select(b => b.Inputs[i]).Quantile((j + 1) / (num_splits + 1));
                        if (double.IsNaN(quantile))
                            break;
                        if (quantile == prunedTrainingPoints.First().Inputs[i] || quantile == prunedTrainingPoints.Last().Inputs[i])
                            continue;
                        double averageFirst = prunedTrainingPoints.TakeWhile(b => b.Inputs[i] < quantile).Select(b => b.Label).Average();
                        double averageLast = prunedTrainingPoints.SkipWhile(b => b.Inputs[i] < quantile).Select(b => b.Label).Average();

                        var sumSquaredErrors = prunedTrainingPoints.TakeWhile(b => b.Inputs[i] < quantile).Select(b => Math.Pow(averageFirst - b.Label, 2)).Sum() +
                                               prunedTrainingPoints.SkipWhile(b => b.Inputs[i] < quantile).Select(b => Math.Pow(averageLast - b.Label, 2)).Sum();
                        //if (level == 0)
                        //{
                        //    Console.WriteLine("  j = " + j);
                        //    Console.WriteLine("   quantile = " + quantile);
                        //    Console.WriteLine("   averageFirst = " + averageFirst);
                        //    Console.WriteLine("   averageLast = " + averageLast);
                        //    Console.WriteLine("   sumSquaredErrors = " + sumSquaredErrors);
                        //}

                        if (sumSquaredErrors < bestSumSquaredErrors)
                        {
                            bestSumSquaredErrors = sumSquaredErrors;
                            bestValue = quantile;
                            bestI = i;
                            //if (level == 0)
                            //{
                            //    Console.WriteLine("  replacing!");
                            //    Console.WriteLine("   bestI = " + bestI);
                            //    Console.WriteLine("   bestSumSquaredErrors = " + bestSumSquaredErrors);
                            //    Console.WriteLine("   bestValue = " + bestValue);
                            //}
                        }
                    }
                }

            bool reachedBottom = false;
            if (bestI == -1)
                reachedBottom = true;

            if (reachedBottom)
                output = trainingPoints.Select(b => b.Label).Average();
            else
            {
                trainingPoints.Sort(Comparer<LabeledDataPoint>.Create((x, y) => x.Inputs[bestI].CompareTo(y.Inputs[bestI])));
                leftTree = new RegressionTree(doNotSplitIfUnderThis, level + 1, useFeature);
                leftTree.Train(trainingPoints.TakeWhile(b => b.Inputs[bestI] < bestValue).ToList());
                rightTree = new RegressionTree(doNotSplitIfUnderThis, level + 1, useFeature);
                rightTree.Train(trainingPoints.SkipWhile(b => b.Inputs[bestI] < bestValue).ToList());
            }
        }

        #endregion Internal Methods
    }
}