namespace EngineLayer.Neo
{
    public class TranslatedParent : Parent
    {
        public TranslatedParent(string id, string seq, int start, int length)
            : base(id, seq)
        {
            Start = start;
            PeptideLength = length;
        }

        public int Start { get; set; }
        public int PeptideLength { get; set; }
        public FusionType Translated { get; }
    }
}