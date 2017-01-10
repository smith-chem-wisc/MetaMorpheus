using System.Collections.Generic;

namespace OldInternalLogic
{
    public class SpecificProteinSelection
    {
        public static bool enabled = false;

        private static HashSet<string> proteinLimits = new HashSet<string>();

        public static bool ConsiderProtein(string accession)
        {
            foreach (var ye in proteinLimits)
                if (accession.Contains(ye))
                    return true;
            return false;
        }

        public static void Add(string v)
        {
            proteinLimits.Add(v);
            enabled = true;
        }

        public static int Count
        {
            get
            {
                return proteinLimits.Count;
            }
        }
    }
}