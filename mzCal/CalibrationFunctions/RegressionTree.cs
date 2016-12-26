using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    internal class RegressionTree
    {
        private RegressionTree leftTree;
        private RegressionTree rightTree;

        private double output = double.NaN;

        private double bestValue;
        private int bestI = -1;

        internal void Train(List<LabeledDataPoint> trainingPoints, int splitLimit, int level)
        {
            //if (trainingPoints.Count > 10000)
            //    Console.WriteLine(" trainingPoints.Count = " + trainingPoints.Count);
            var averageOutputs = trainingPoints.Select(b => b.output).Average();
            if (trainingPoints.Count() < splitLimit)
            {
                output = averageOutputs;
                return;
            }
            var bestSumSquaredErrors = trainingPoints.Select(b => Math.Pow(averageOutputs - b.output, 2)).Sum();
            //int bestCount = int.MaxValue;

            //Random rand = new Random();
            //var prunedTrainingPoints = trainingPoints.Where(b => rand.NextDouble() > 2.0 / 3).ToList();

            var prunedTrainingPoints = trainingPoints;

            // For every variable, try to find the best split
            for (int i = 0; i < trainingPoints[0].inputs.Length; i++)
            {
                prunedTrainingPoints.Sort(Comparer<LabeledDataPoint>.Create((x, y) => x.inputs[i].CompareTo(y.inputs[i])));

                // Now that have sorted array, find split points:
                int num_splits = Math.Min(15, prunedTrainingPoints.Count - 1);
                for (double j = 0; j < num_splits; j++)
                {
                    //Console.WriteLine(" " + j + "th quantile");
                    double quantile = prunedTrainingPoints.Select(b => b.inputs[i]).Quantile((j + 1) / (num_splits + 1));
                    if (quantile == prunedTrainingPoints.First().inputs[i] || quantile == prunedTrainingPoints.Last().inputs[i])
                        continue;
                    //Console.WriteLine(" quantile = " + quantile);
                    double averageFirst = prunedTrainingPoints.TakeWhile(b => b.inputs[i] < quantile).Select(b => b.output).Average();
                    double averageLast = prunedTrainingPoints.SkipWhile(b => b.inputs[i] < quantile).Select(b => b.output).Average();
                    var sumSquaredErrors = prunedTrainingPoints.TakeWhile(b => b.inputs[i] < quantile).Select(b => Math.Pow(averageFirst - b.output, 2)).Sum() +
                                           prunedTrainingPoints.SkipWhile(b => b.inputs[i] < quantile).Select(b => Math.Pow(averageLast - b.output, 2)).Sum();
                    if (sumSquaredErrors < bestSumSquaredErrors)
                    {
                        //Console.WriteLine(" error is better!");
                        //Console.WriteLine(" sumSquaredErrors = " + sumSquaredErrors);
                        //Console.WriteLine(" bestSumSquaredErrors = " + bestSumSquaredErrors);
                        bestSumSquaredErrors = sumSquaredErrors;
                        bestValue = quantile;
                        bestI = i;
                        //bestCount = prunedTrainingPoints.TakeWhile(b => b.inputs[i] < quantile).Count();
                    }
                }
            }

            bool reachedBottom = false;
            if (bestI == -1)
                reachedBottom = true;

            if (reachedBottom)
            {
                if (trainingPoints.Count > 10)
                {
                    //Console.WriteLine(trainingPoints.Count);
                    //for (int i = 0; i < trainingPoints[0].inputs.Count(); i++)
                    //Console.WriteLine(string.Join(", ", trainingPoints.Select(b => b.inputs[i])));
                    //Console.WriteLine(string.Join(", ", trainingPoints.Select(b => b.output)));
                }
                output = trainingPoints.Select(b => b.output).Average();
            }
            else
            {
                //if (level <= 2)
                //    Console.WriteLine(" level " + level + " bestI = " + bestI);
                //Console.WriteLine(" Splitting");
                //Console.WriteLine(" bestI = " + bestI);
                //Console.WriteLine(" bestValue = " + bestValue);
                trainingPoints.Sort(Comparer<LabeledDataPoint>.Create((x, y) => x.inputs[bestI].CompareTo(y.inputs[bestI])));
                leftTree = new RegressionTree();
                leftTree.Train(trainingPoints.TakeWhile(b => b.inputs[bestI] < bestValue).ToList(), splitLimit, level + 1);
                rightTree = new RegressionTree();
                rightTree.Train(trainingPoints.SkipWhile(b => b.inputs[bestI] < bestValue).ToList(), splitLimit, level + 1);
            }
        }

        public RegressionTree()
        {
        }

        internal double predict(double[] inputs)
        {
            if (double.IsNaN(output))
            {
                if (inputs[bestI] < bestValue)
                    return leftTree.predict(inputs);
                else
                    return rightTree.predict(inputs);
            }
            else
            {
                return output;
            }
        }
    }
}