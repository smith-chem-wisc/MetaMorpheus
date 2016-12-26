using MathNet.Numerics.Statistics;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    internal class MyBinaryTree
    {
        private MyBinaryTree left;
        private MyBinaryTree right;
        private double threshold = double.NaN;
        private double value = double.NaN;

        public MyBinaryTree(List<LabeledDataPoint> trainingList)
        {
            var selected = trainingList.Select(b => b.inputs[2]).ToList();
            if (selected.Min() == selected.Max())
                value = trainingList.Select(b => b.output).Average();
            else
            {
                threshold = selected.Median();
                if (threshold == selected.Min() || threshold == selected.Max())
                    value = trainingList.Select(b => b.output).Average();
                else
                {
                    left = new MyBinaryTree(trainingList.Where(b => b.inputs[2] < threshold).ToList());
                    right = new MyBinaryTree(trainingList.Where(b => b.inputs[2] >= threshold).ToList());
                }
            }
        }

        internal double GetValue(double v)
        {
            if (double.IsNaN(value))
            {
                if (v < threshold)
                    return left.GetValue(v);
                else
                    return right.GetValue(v);
            }
            else
                return value;
        }
    }
}