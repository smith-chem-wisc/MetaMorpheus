using System;
using System.Collections.Generic;

namespace MetaMorpheus
{
    internal class MyDoubleEqualityComparer : IEqualityComparer<double>
    {
        private double v;

        public MyDoubleEqualityComparer(double v)
        {
            this.v = v;
        }

        public bool Equals(double x, double y)
        {
            return Math.Abs(x - y) <= v;
        }

        public int GetHashCode(double obj)
        {
            return obj.GetHashCode();
        }
    }
}