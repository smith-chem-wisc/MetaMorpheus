namespace EngineLayer.Neo
{
    public class TranslatedParent : Parent
    {
        #region Public Constructors

        public TranslatedParent(string id, string seq, int start, int length) : base(id, seq)
        {
            this.start = start;
            this.peptideLength = length;
        }

        #endregion Public Constructors

        #region Public Properties

        public int start { get; set; }
        public int peptideLength { get; set; }
        public FusionCandidate.FusionType translated { get; }

        #endregion Public Properties
    }
}