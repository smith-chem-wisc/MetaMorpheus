using MathNet.Numerics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System.Collections.Generic;
using System.Linq;

namespace TaskLayer
{
    public class SeparateMzLearner : ILearner<double>
    {
        #region Private Fields

        private ILearner<double> ok;

        #endregion Private Fields

        #region Public Constructors

        public SeparateMzLearner(ILearner<double> ok)
        {
            this.ok = ok;
        }

        #endregion Public Constructors

        #region Public Methods

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            List<double> targetCoeffs = new List<double>();
            var dps = new List<double>();
            double prevRT = observations[0, 1];
            Dictionary<double, List<(double[], double)>> dict = new Dictionary<double, List<(double[], double)>>();
            for (int i = 0; i < observations.RowCount; i++)
            {
                double[] row = observations.Row(i);
                if (dict.ContainsKey(row[1]))
                    dict[row[1]].Add((row, targets[i]));
                else
                    dict[row[1]] = new List<(double[], double)> { (row, targets[i]) };
            }

            foreach (var kvp in dict.Where(b => b.Value.Count > 1))
            {
                // Train only on the MZ value!!!! Since already considering a specific scan

                var f = Fit.LinearMultiDimFunc(kvp.Value.Select(b => new[] { b.Item1[0] }).ToArray(), kvp.Value.Select(b => b.Item2).ToArray());
                foreach (var dp in kvp.Value)
                {
                    dps.AddRange(dp.Item1.Skip(1).ToArray());
                    targetCoeffs.Add(f(new double[] { 1 }));
                }
            }

            // Skip the intensities for now!!!

            F64Matrix kook = new F64Matrix(dps.ToArray(), targetCoeffs.Count, dps.Count / targetCoeffs.Count);
            IPredictorModel<double> predictMzCoeff = ok.Learn(kook, targetCoeffs.ToArray());
            return new MySeparatePredictorModel(predictMzCoeff);
        }

        public override string ToString()
        {
            return "SeparateMzLearner " + ok.ToString();
        }

        #endregion Public Methods
    }

    internal class PredictMzCoeff : IPredictorModel<double>
    {
        #region Private Fields

        private F64Matrix observations;
        private double[] targets;

        #endregion Private Fields

        #region Public Constructors

        public PredictMzCoeff(F64Matrix observations, double[] targets)
        {
            this.observations = observations;
            this.targets = targets;
        }

        #endregion Public Constructors

        #region Public Methods

        public double[] GetRawVariableImportance()
        {
            throw new System.NotImplementedException();
        }

        public Dictionary<string, double> GetVariableImportance(Dictionary<string, int> featureNameToIndex)
        {
            throw new System.NotImplementedException();
        }

        public double Predict(double[] observation)
        {
            throw new System.NotImplementedException();
        }

        #endregion Public Methods
    }

    internal class MySeparatePredictorModel : IPredictorModel<double>
    {
        #region Private Fields

        private IPredictorModel<double> predictMzCoeff;

        #endregion Private Fields

        #region Public Constructors

        public MySeparatePredictorModel(IPredictorModel<double> predictMzCoeff)
        {
            this.predictMzCoeff = predictMzCoeff;
        }

        #endregion Public Constructors

        #region Public Methods

        public double[] GetRawVariableImportance()
        {
            throw new System.NotImplementedException();
        }

        public Dictionary<string, double> GetVariableImportance(Dictionary<string, int> featureNameToIndex)
        {
            throw new System.NotImplementedException();
        }

        public double Predict(double[] observation)
        {
            return predictMzCoeff.Predict(observation.Skip(1).ToArray()) * observation[0];
        }

        #endregion Public Methods
    }
}