//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace mzCal
//{
//    public class RandomForest : CalibrationFunction
//    {
//        RegressionTree[] RegressionTrees;
//        bool[] useFeature;
//        int numFeatures;

//        public RandomForest(List<LabeledDataPoint> trainingList, int numTrees, int splitLimit, bool[] useFeature)
//        {
//            this.useFeature = useFeature;
//            var rand = new Random();

//            RegressionTrees = new RegressionTree[numTrees];

//            for (int i = 0; i < numTrees; i++)
//                RegressionTrees[i] = new RegressionTree();

//            numFeatures = useFeature.Where(b => b == true).Count();

//            Parallel.For(0, numTrees, i =>
//            {
//                List<LabeledDataPoint> subsampledTrainingPoints = new List<LabeledDataPoint>();
//                for (int j = 0; j < trainingList.Count; j++)
//                {
//                    int index = rand.Next(trainingList.Count);
//                    var yeesh = new LabeledDataPoint(IndexMap(trainingList[index].inputs), trainingList[index].output);
//                    subsampledTrainingPoints.Add(yeesh);
//                }
//                RegressionTrees[i].Train(subsampledTrainingPoints, splitLimit, 0);
//            });
//        }

//        public override double Predict(double[] input)
//        {
//            return RegressionTrees.Select(b => b.predict(IndexMap(input))).Average();
//        }

//        private double[] IndexMap(double[] input)
//        {
//            double[] output = new double[numFeatures];
//            int featInd = 0;
//            for (int k = 0; k < useFeature.Length; k++)
//                if (useFeature[k])
//                {
//                    output[featInd] = input[k];
//                    featInd++;
//                }
//            return output;
//        }
//    }
//}