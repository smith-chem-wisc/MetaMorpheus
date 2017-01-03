using System;

namespace IndexSearchAndAnalyze
{
    public class SearchMode
    {
        private Func<double, bool> p;

        public SearchMode(string FileNameAddition, Func<double, bool> p)
        {
            this.FileNameAddition = FileNameAddition;
            this.p = p;
        }

        public string FileNameAddition { get; internal set; }

        internal bool Accepts(double dd)
        {
            return p(dd);
        }
    }
}