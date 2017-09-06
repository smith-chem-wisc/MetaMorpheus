using MathNet.Numerics.Statistics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Calibration
{
    internal class RegressionTree : ILearner<double>
    {
        #region Private Fields

        private readonly int level;
        private readonly bool[] useFeature;
        private readonly int doNotSplitIfUnderThis;

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

        #region Public Methods

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            //Console.WriteLine("        In Train");

            //Console.WriteLine("  All points in actual trainingList");
            //Console.WriteLine(string.Join(Environment.NewLine, trainingPoints.OrderBy(b => b.Inputs[0]).Select(b => "      " + b.Inputs[0] + " " + b.Label)));

            //Console.WriteLine(string.Join(Environment.NewLine, trainingPoints.Select(b => "          " + b.Inputs[0] + " , " + b.Label)));

            var averageOutputs = targets.Average();

            //Console.WriteLine("          Average outputs: " + averageOutputs);

            if (targets.Length < doNotSplitIfUnderThis)
            {
                output = averageOutputs;
                //Console.WriteLine("          count too low! output = " + output + " level= " + level);
                return new RegressionTreePredictorModel(output);
            }

            double bestSumSquaredErrors = 0;
            var prunedTrainingPoints = new List<IHasInputsAndOutputs>();
            for (int i = 0; i < targets.Length; i++)
            {
                bestSumSquaredErrors += Math.Pow(averageOutputs - targets[i], 2);
                prunedTrainingPoints.Add(new LabeledDataPoint(observations.Row(i), targets[i]));
            }

            // For every variable, try to find the best split
            for (int i = 0; i < useFeature.Length; i++)
                if (useFeature[i])
                {
                    prunedTrainingPoints.Sort(Comparer<IHasInputsAndOutputs>.Create((x, y) => x.Inputs[i].CompareTo(y.Inputs[i])));
                    int num_splits = Math.Min(15, prunedTrainingPoints.Count - 1);
                    for (double j = 0; j < num_splits; j++)
                    {
                        double quantile = prunedTrainingPoints.Select(b => b.Inputs[i]).Quantile((j + 1) / (num_splits + 1));

                        //Console.WriteLine("quantile: " + quantile + " level = " + level);

                        if (double.IsNaN(quantile))
                            break;
                        if (quantile == prunedTrainingPoints.First().Inputs[i] || quantile == prunedTrainingPoints.Last().Inputs[i])
                            continue;
                        double averageFirst = prunedTrainingPoints.TakeWhile(b => b.Inputs[i] < quantile).Select(b => b.Label).Average();
                        double averageLast = prunedTrainingPoints.SkipWhile(b => b.Inputs[i] < quantile).Select(b => b.Label).Average();

                        var sumSquaredErrors = prunedTrainingPoints.TakeWhile(b => b.Inputs[i] < quantile).Select(b => Math.Pow(averageFirst - b.Label, 2)).Sum() +
                                               prunedTrainingPoints.SkipWhile(b => b.Inputs[i] < quantile).Select(b => Math.Pow(averageLast - b.Label, 2)).Sum();

                        if (sumSquaredErrors < bestSumSquaredErrors)
                        {
                            bestSumSquaredErrors = sumSquaredErrors;
                            bestValue = quantile;
                            bestI = i;
                            //Console.WriteLine("bestValue: " + bestValue + " level = " + level);
                        }
                    }
                }

            bool reachedBottom = false;
            if (bestI == -1)
                reachedBottom = true;

            if (reachedBottom)
            {
                output = targets.Average();
                //Console.WriteLine("          count too low! output = " + output + " level= " + level);
                return new RegressionTreePredictorModel(output);
            }
            else
            {
                prunedTrainingPoints.Sort(Comparer<IHasInputsAndOutputs>.Create((x, y) => x.Inputs[bestI].CompareTo(y.Inputs[bestI])));

                IPredictorModel<double> leftTreePM;

                IPredictorModel<double> rightTreePM;

                {
                    var leftTree = new RegressionTree(doNotSplitIfUnderThis, level + 1, useFeature);
                    var forLeft = prunedTrainingPoints.TakeWhile(b => b.Inputs[bestI] < bestValue).ToList();
                    var theLeftArray = forLeft.SelectMany(b => b.Inputs).ToArray();
                    leftTreePM = leftTree.Learn(new F64Matrix(theLeftArray, forLeft.Count, theLeftArray.Length / forLeft.Count), forLeft.Select(b => b.Label).ToArray());
                }
                {
                    var rightTree = new RegressionTree(doNotSplitIfUnderThis, level + 1, useFeature);
                    var forRight = prunedTrainingPoints.SkipWhile(b => b.Inputs[bestI] < bestValue).ToList();
                    var theRightArray = forRight.SelectMany(b => b.Inputs).ToArray();
                    rightTreePM = rightTree.Learn(new F64Matrix(theRightArray, forRight.Count, theRightArray.Length / forRight.Count), forRight.Select(b => b.Label).ToArray());
                }
                return new RegressionTreePredictorModel(bestValue, bestI, leftTreePM as RegressionTreePredictorModel, rightTreePM as RegressionTreePredictorModel);
            }
        }

        #endregion Public Methods

        //Console.WriteLine("        Done With Train");
    }

    internal class LabeledDataPoint : IHasInputsAndOutputs
    {
        #region Public Constructors

        public LabeledDataPoint(double[] v1, double v2)
        {
            this.Inputs = v1;
            this.Label = v2;
        }

        #endregion Public Constructors

        #region Public Properties

        public double[] Inputs { get; private set; }

        public double Label { get; private set; }

        #endregion Public Properties
    }

    internal class RegressionTreePredictorModel : IPredictorModel<double>
    {
        #region Private Fields

        private double output = double.NaN;
        private double bestValue;
        private int bestI;
        private RegressionTreePredictorModel leftTree;
        private RegressionTreePredictorModel rightTree;

        #endregion Private Fields

        #region Public Constructors

        public RegressionTreePredictorModel(double bestValue, int bestI, RegressionTreePredictorModel leftTree, RegressionTreePredictorModel rightTree)
        {
            this.bestValue = bestValue;
            this.bestI = bestI;
            this.leftTree = leftTree;
            this.rightTree = rightTree;
        }

        public RegressionTreePredictorModel(double output)
        {
            this.output = output;
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
            if (double.IsNaN(output))
            {
                if (observation[bestI] < bestValue)
                    return leftTree.Predict(observation);
                else
                    return rightTree.Predict(observation);
            }
            else
                return output;
        }

        #endregion Public Methods
    }
}