using System;

namespace MetaMorpheusLogic
{
    public class TransformFunction
    {
        internal int numOutputs;
        private Func<double[], double[]> tf;
        internal string name;

        internal double[] Transform(double[] t)
        {
            return tf(t);
        }

        public TransformFunction(Func<double[], double[]> tf, int numOutputs, string name)
        {
            this.tf = tf;
            this.numOutputs = numOutputs;
            this.name = name;
        }
    }
}