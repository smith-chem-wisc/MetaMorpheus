namespace IndexSearchAndAnalyze
{
    public abstract class SearchMode
    {
        public SearchMode(string fileNameAddition)
        {
            FileNameAddition = fileNameAddition;
        }

        public string FileNameAddition { get; internal set; }

        public abstract bool Accepts(double scanPrecursorMass, double peptideMass);
    }
}