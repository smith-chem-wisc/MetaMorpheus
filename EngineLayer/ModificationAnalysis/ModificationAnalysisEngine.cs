using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ModificationAnalysisEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly int searchModesCount;
        private readonly List<Psm>[] newPsms;

        #endregion Private Fields

        #region Public Constructors

        public ModificationAnalysisEngine(List<Psm>[] newPsms, int searchModesCount, List<string> nestedIds) : base(nestedIds)
        {
            this.searchModesCount = searchModesCount;
            this.newPsms = newPsms;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ModificationAnalysisResults myAnalysisResults = new ModificationAnalysisResults(this);

            Dictionary<string, int>[] nonLocMods = new Dictionary<string, int>[searchModesCount];
            Dictionary<string, int>[] allModsOnProteins = new Dictionary<string, int>[searchModesCount];
            Dictionary<string, int>[] allmodsSeenAndLocalized = new Dictionary<string, int>[searchModesCount];

            for (int j = 0; j < searchModesCount; j++)
            {
                Status("Running modification analysis...", nestedIds);

                Dictionary<string, int> modsSeen = new Dictionary<string, int>();

                foreach (var groupOfPsmsWithSameBaseSequence in newPsms[j].Where(b => b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.BaseSequence != null && b.ModsIdentified != null && b.FullSequence==null).GroupBy(b => (b.BaseSequence)))
                {
                    Dictionary<string, int> set = new Dictionary<string, int>();
                    foreach (var psm in groupOfPsmsWithSameBaseSequence)
                    {
                        foreach (var modWithCount in psm.ModsIdentified)
                        {
                            if (set.ContainsKey(modWithCount.Key))
                            {
                                if (set[modWithCount.Key] < modWithCount.Value)
                                    set[modWithCount.Key] = modWithCount.Value;
                            }
                            else set[modWithCount.Key] = modWithCount.Value;
                        }
                    }

                    foreach (var modSeen in set)
                    {
                        if (modsSeen.ContainsKey(modSeen.Key))
                            modsSeen[modSeen.Key] += modSeen.Value;
                        else
                            modsSeen.Add(modSeen.Key, modSeen.Value);
                    }
                }
                nonLocMods[j] = modsSeen;

                HashSet<(string, string, int)> modsOnProteins = new HashSet<(string, string, int)>();
                HashSet<(string, string, int)> modsSeenAndLocalized = new HashSet<(string, string, int)>();

                foreach (var psm in newPsms[j].Where(b => b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.ProteinAccesion != null && b.OneBasedEndResidueInProtein != null && b.OneBasedStartResidueInProtein != null))
                {
                    var singlePeptide = psm.CompactPeptides.First().Value.Item2.First();
                    foreach (var modInProtein in singlePeptide.Protein.OneBasedPossibleLocalizedModifications.Where(b => b.Key >= singlePeptide.OneBasedStartResidueInProtein && b.Key <= singlePeptide.OneBasedEndResidueInProtein))
                    {
                        foreach (var huh in modInProtein.Value)
                        {
                            modsOnProteins.Add((singlePeptide.Protein.Accession, huh.id, modInProtein.Key));
                        }
                    }

                    if (psm.FullSequence != null)
                    {
                        foreach (var nice in singlePeptide.allModsOneIsNterminus)
                        {
                            int locInProtein;
                            if (nice.Key == 1)
                                locInProtein = singlePeptide.OneBasedStartResidueInProtein;
                            else if (nice.Key == singlePeptide.Length + 2)
                                locInProtein = singlePeptide.OneBasedEndResidueInProtein;
                            else
                                locInProtein = singlePeptide.OneBasedStartResidueInProtein + nice.Key - 2;

                            modsSeenAndLocalized.Add((singlePeptide.Protein.Accession, nice.Value.id, locInProtein));
                        }
                    }
                }

                allModsOnProteins[j] = modsOnProteins.GroupBy(b => b.Item2).ToDictionary(b => b.Key, b => b.Count());
                allmodsSeenAndLocalized[j] = modsSeenAndLocalized.GroupBy(b => b.Item2).ToDictionary(b => b.Key, b => b.Count());
            }
            myAnalysisResults.NonLocalizedModsSeen = nonLocMods;
            myAnalysisResults.ModsSeenAndLocalized = allmodsSeenAndLocalized;
            myAnalysisResults.AllModsOnProteins = allModsOnProteins;
            return myAnalysisResults;
        }

        #endregion Protected Methods
    }
}