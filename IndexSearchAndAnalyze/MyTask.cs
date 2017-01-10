using MetaMorpheus;
using Spectra;
using System.Collections.Generic;
using System.Collections.ObjectModel;

namespace IndexSearchAndAnalyze
{
    public abstract class MyTask
    {
        public MyTask(int selectedIndex)
        {
            switch (selectedIndex)
            {
                case 0:
                    taskType = MyTaskEnum.Calibrate;
                    break;

                case 1:
                    taskType = MyTaskEnum.Search;
                    break;

                case 2:
                    taskType = MyTaskEnum.GPTMD;
                    break;
            }
        }

        public MyTaskEnum taskType { get; internal set; }
        public bool IsMySelected { get; set; }
        public string output_folder { get; private set; }
        public int maxMissedCleavages { get; set; }
        public Protease protease { get; set; }
        public int maxModificationIsoforms { get; set; }
        public InitiatorMethionineBehavior initiatorMethionineBehavior { get; set; }
        public Tolerance productMassTolerance { get; set; }
        public bool bIons { get; set; }
        public bool yIons { get; set; }

        public List<ModListForSearch> listOfModListsForSearch { get; set; }

        public abstract MyTaskResults DoTask(ObservableCollection<RawData> completeRawFileListCollection, ObservableCollection<XMLdb> completeXmlDbList, AllTasksParams po);

        protected static void GenerateModsFromStrings(List<string> listOfXMLdbs, List<MorpheusModification> modsKnown, out Dictionary<string, List<MorpheusModification>> modsToLocalize, out HashSet<string> modsInXMLtoTrim)
        {
            modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
            var modsInXML = ProteomeDatabaseReader.ReadXMLmodifications(listOfXMLdbs);
            modsInXMLtoTrim = new HashSet<string>(modsInXML);
            foreach (var knownMod in modsKnown)
                if (modsInXML.Contains(knownMod.NameInXML))
                {
                    if (modsToLocalize.ContainsKey(knownMod.NameInXML))
                        modsToLocalize[knownMod.NameInXML].Add(knownMod);
                    else
                        modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification>() { knownMod });
                    modsInXMLtoTrim.Remove(knownMod.NameInXML);
                }
        }

        public void setOutputFolder(string thisOutputPath)
        {
            this.output_folder = thisOutputPath;
        }
    }
}