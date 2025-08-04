using EngineLayer.GlycoSearch;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.ModSearch
{
    public class ModSearchEngine : GlycoSearchEngine
    {

        // For general modification settings
        private readonly int TopN; // DDA top Peak number
        private readonly int MaxModNumber;
        List<(string, string)> ListOfInterestedMods; // List of interested modifications


        public ModSearchEngine(List<SpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans,
            List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition,
            CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            string oglycanDatabase, string nglycanDatabase, List<(string, string)> ListOfInterestedMods, int modSearchTopNum,

        int maxModNum, bool oxoniumIonFilter, List<string> nestedIds)
            : base(globalCsms, listOfSortedms2Scans, peptideIndex, fragmentIndex, secondFragmentIndex, currentPartition, commonParameters, fileSpecificParameters, oglycanDatabase, nglycanDatabase, GlycoSearchType.ModSearch, modSearchTopNum, maxModNum, oxoniumIonFilter, nestedIds)
        {
            this.TopN = modSearchTopNum;
            this.MaxModNumber = maxModNum;
            this.ListOfInterestedMods = ListOfInterestedMods;

            //Load glycan databases and build the modBox 
            ModBox.GlobalModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => ListOfInterestedMods.Contains((b.ModificationType, b.IdWithMotif))).ToArray();
            ModBoxes = ModBox.BuildModBoxes(MaxModNumber, false).OrderBy(p => p.Mass).ToArray();
        }


        protected override HashSet<string> GetMotifs(int id = -1)
        {
            if (id == -1)
            {
                return ModBox.GlobalModifications.Select(p => p.Target.ToString()).ToHashSet();
            }

            return  ModBoxes[id].ModIds
                .Select(p => ModBox.GlobalModifications[p].Target.ToString())
                .ToHashSet();
        }

    }
}
