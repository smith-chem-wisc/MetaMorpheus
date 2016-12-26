namespace MetaMorpheus
{
    public enum CleavageSpecificity
    {
        None,
        Semi,
        SemiN,
        SemiC,
        Full
    }

    internal static class CleavageSpecificityExtensions
    {
        public static int GetMinNumberTermini(this CleavageSpecificity cleavageSpecifity)
        {
            switch (cleavageSpecifity)
            {
                case CleavageSpecificity.None:
                    return 0;

                case CleavageSpecificity.Semi:
                case CleavageSpecificity.SemiN:
                case CleavageSpecificity.SemiC:
                    return 1;

                case CleavageSpecificity.Full:
                default:
                    return 2;
            }
        }
    }
}