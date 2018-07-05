using Chemistry;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.ModificationAnalysis
{
    public class ModificationAnalysisEngine : MetaMorpheusEngine
    {
        private readonly List<PeptideSpectralMatch> NewPsms;

        public ModificationAnalysisEngine(List<PeptideSpectralMatch> newPsms, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            NewPsms = newPsms;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Running modification analysis...");
            ModificationAnalysisResults myAnalysisResults = new ModificationAnalysisResults(this);

            var confidentTargetPsms = NewPsms.Where(b => b.FdrInfo.QValue <= 0.01 && !b.IsDecoy).ToList();

            // For the database ones, only need un-ambiguous protein and location in protein
            var forObserved = confidentTargetPsms
                .Where(b => b.ProteinAccesion != null && b.OneBasedEndResidueInProtein != null && b.OneBasedStartResidueInProtein != null);

            // For the unambiguously localized ones, need FullSequence and un-ambiguous protein and location in protein
            var forUnambiguouslyLocalized = confidentTargetPsms
                .Where(b => b.FullSequence != null && b.ProteinAccesion != null && b.OneBasedEndResidueInProtein != null && b.OneBasedStartResidueInProtein != null);

            // For the localized but ambiguous ones, need FullSequence
            var forAmbiguousButLocalized = confidentTargetPsms
                .Where(b => b.FullSequence != null && !(b.ProteinAccesion != null && b.OneBasedEndResidueInProtein != null && b.OneBasedStartResidueInProtein != null))
                .GroupBy(b => b.FullSequence);

            // For unlocalized but identified modifications, skip ones with full sequences!
            var forUnlocalized = confidentTargetPsms
                .Where(b => b.BaseSequence != null && b.FullSequence == null && b.ModsIdentified != null)
                .GroupBy(b => (b.BaseSequence, string.Join(" ", b.ModsIdentified.Values.OrderBy(c => c))));

            // For chemical formulas of modifications, skip ones with full sequences and identified mods!
            var forChemicalFormulas = confidentTargetPsms
                .Where(b => b.BaseSequence != null && b.FullSequence == null && b.ModsIdentified == null && b.ModsChemicalFormula != null)
                .GroupBy(b => (b.BaseSequence, b.ModsChemicalFormula));

            // We do not want to double-count modifications. Hence the HashSet!!!
            HashSet<(string, string, int)> modsOnProteins = new HashSet<(string, string, int)>();
            foreach (var psm in forObserved)
            {
                var singlePeptide = psm.CompactPeptides.First().Value.Item2.First();
                foreach (var modInProtein in singlePeptide.Protein.OneBasedPossibleLocalizedModifications.Where(b => b.Key >= singlePeptide.OneBasedStartResidueInProtein && b.Key <= singlePeptide.OneBasedEndResidueInProtein))
                {
                    foreach (var huh in modInProtein.Value)
                    {
                        modsOnProteins.Add((singlePeptide.Protein.Accession, huh.id, modInProtein.Key));
                    }
                }
            }

            // We do not want to double-count modifications. Hence the HashSet!!!
            HashSet<(string, string, int)> modsSeenAndLocalized = new HashSet<(string, string, int)>();
            foreach (var psm in forUnambiguouslyLocalized)
            {
                var singlePeptide = psm.CompactPeptides.First().Value.Item2.First();
                foreach (var nice in singlePeptide.AllModsOneIsNterminus)
                {
                    int locInProtein;
                    if (nice.Key == 1)
                    {
                        locInProtein = singlePeptide.OneBasedStartResidueInProtein;
                    }
                    else if (nice.Key == singlePeptide.Length + 2)
                    {
                        locInProtein = singlePeptide.OneBasedEndResidueInProtein;
                    }
                    else
                    {
                        locInProtein = singlePeptide.OneBasedStartResidueInProtein + nice.Key - 2;
                    }
                    modsSeenAndLocalized.Add((singlePeptide.Protein.Accession, nice.Value.id, locInProtein));
                }
            }

            // Might have some double counting...
            Dictionary<string, int> ambiguousButLocalizedModsSeen = new Dictionary<string, int>();
            foreach (var representativePsm in forAmbiguousButLocalized.Select(b => b.First()))
            {
                foreach (var modCountKvp in representativePsm.ModsIdentified)
                {
                    if (ambiguousButLocalizedModsSeen.ContainsKey(modCountKvp.Key))
                    {
                        ambiguousButLocalizedModsSeen[modCountKvp.Key] += modCountKvp.Value;
                    }
                    else
                    {
                        ambiguousButLocalizedModsSeen.Add(modCountKvp.Key, modCountKvp.Value);
                    }
                }
            }

            // Might have some double counting...
            Dictionary<string, int> unlocalizedMods = new Dictionary<string, int>();
            foreach (var representativePsm in forUnlocalized.Select(b => b.First()))
            {
                foreach (var modCountKvp in representativePsm.ModsIdentified)
                {
                    if (unlocalizedMods.ContainsKey(modCountKvp.Key))
                    {
                        unlocalizedMods[modCountKvp.Key] += modCountKvp.Value;
                    }
                    else
                    {
                        unlocalizedMods.Add(modCountKvp.Key, modCountKvp.Value);
                    }
                }
            }

            // Might have some double counting...
            Dictionary<ChemicalFormula, int> unlocalizedFormulas = new Dictionary<ChemicalFormula, int>();
            foreach (var representativePsm in forChemicalFormulas.Select(b => b.First()))
            {
                if (unlocalizedFormulas.ContainsKey(representativePsm.ModsChemicalFormula))
                {
                    unlocalizedFormulas[representativePsm.ModsChemicalFormula] += 1;
                }
                else
                {
                    unlocalizedFormulas.Add(representativePsm.ModsChemicalFormula, 1);
                }
            }

            myAnalysisResults.AllModsOnProteins = modsOnProteins.GroupBy(b => b.Item2).ToDictionary(b => b.Key, b => b.Count());
            myAnalysisResults.ModsSeenAndLocalized = modsSeenAndLocalized.GroupBy(b => b.Item2).ToDictionary(b => b.Key, b => b.Count());
            myAnalysisResults.AmbiguousButLocalizedModsSeen = ambiguousButLocalizedModsSeen;
            myAnalysisResults.UnlocalizedMods = unlocalizedMods;
            myAnalysisResults.UnlocalizedFormulas = unlocalizedFormulas;

            return myAnalysisResults;
        }
    }
}