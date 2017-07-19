using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ModificationAnalysisEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private readonly List<MassDiffAcceptor> searchModes;
        private readonly List<SingleScanManyPeptidesMatch>[] newPsms;

        #endregion Private Fields

        #region Public Constructors

        public ModificationAnalysisEngine(List<SingleScanManyPeptidesMatch>[] newPsms, List<MassDiffAcceptor> searchModes, List<string> nestedIds) : base(nestedIds)
        {
            this.searchModes = searchModes;
            this.newPsms = newPsms;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ModificationAnalysisResults myAnalysisResults = new ModificationAnalysisResults(this);

            Dictionary<string, int>[] allModsSeen = new Dictionary<string, int>[searchModes.Count];
            Dictionary<string, int>[] allModsOnPeptides = new Dictionary<string, int>[searchModes.Count];

            for (int j = 0; j < searchModes.Count; j++)
            {
                Status("Running modification analysis...", nestedIds);

                Dictionary<string, int> modsSeen = new Dictionary<string, int>();
                Dictionary<string, int> modsOnPeptides = new Dictionary<string, int>();

                foreach (var highConfidencePSM in newPsms[j].Where(b => b.FdrInfo.QValue <= 0.01 && !b.MostProbableProteinInfo.IsDecoy).GroupBy(b => b.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence).Select(b => b.FirstOrDefault()))
                {
                    var singlePeptide = highConfidencePSM.MostProbableProteinInfo.PeptidesWithSetModifications.First();
                    var modsIdentified = singlePeptide.allModsOneIsNterminus;
                    foreach (var modSeen in modsIdentified)
                    {
                        if (modsSeen.ContainsKey(modSeen.Value.id))
                            modsSeen[modSeen.Value.id]++;
                        else
                            modsSeen.Add(modSeen.Value.id, 1);
                    }
                    var modsInProtein = singlePeptide.Protein.OneBasedPossibleLocalizedModifications.Where(b => b.Key >= singlePeptide.OneBasedStartResidueInProtein && b.Key <= singlePeptide.OneBasedEndResidueInProtein).SelectMany(b => b.Value);
                    foreach (var modInProtein in modsInProtein)
                    {
                        if (modsOnPeptides.ContainsKey(modInProtein.id))
                            modsOnPeptides[modInProtein.id]++;
                        else
                            modsOnPeptides.Add(modInProtein.id, 1);
                    }
                }
                allModsSeen[j] = modsSeen;
                allModsOnPeptides[j] = modsOnPeptides;
            }
            myAnalysisResults.AllModsSeen = allModsSeen;
            myAnalysisResults.AllModsOnPeptides = allModsOnPeptides;
            return myAnalysisResults;
        }

        #endregion Protected Methods

    }
}