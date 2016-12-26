using System.Collections.Generic;

namespace MetaMorpheus
{
    public class SpectraLimits
    {
        public static int Count
        {
            get
            {
                return spectraLimits.Count;
            }
        }

        private static HashSet<int> spectraLimits = new HashSet<int>();
        public static bool limit = false;

        public static void Add(int v)
        {
            spectraLimits.Add(v);
            limit = true;
        }

        public static bool Contains(int spectrumNumber)
        {
            return spectraLimits.Contains(spectrumNumber);
        }
    }
}