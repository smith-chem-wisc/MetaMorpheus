using System;

namespace InternalLogicCalibration
{
    public class TransformFunction
    {
        #region Internal Fields

        internal int numOutputs;

        #endregion Internal Fields

        #region Private Fields

        private readonly Func<double[], double[]> tf;

        #endregion Private Fields

        #region Public Constructors

        public TransformFunction(Func<double[], double[]> tf, int numOutputs)
        {
            this.tf = tf;
            this.numOutputs = numOutputs;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal double[] Transform(double[] t)
        {
            return tf(t);
        }

        #endregion Internal Methods
    }
}