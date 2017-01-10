using System;

namespace IndexSearchAndAnalyze
{
    public class FunctionSearchMode : SearchMode
    {
        private Func<double, double, bool> p;

        public FunctionSearchMode(string FileNameAddition, Func<double, double, bool> p) : base(FileNameAddition)
        {
            this.FileNameAddition = FileNameAddition;
            this.p = p;
        }

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            return p(scanPrecursorMass, peptideMass);
        }
    }
}