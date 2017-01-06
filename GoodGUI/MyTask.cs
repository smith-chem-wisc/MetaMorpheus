using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using MetaMorpheus;
using static GoodGUI.MainWindow;
using IndexSearchAndAnalyze;
using Spectra;
using System.IO;
using MassSpectrometry;
using IO.MzML;
using IO.Thermo;
using System.Linq;
using mzCal;
using Proteomics;
using FragmentGeneration;
using System.Globalization;

namespace GoodGUI
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
        public bool IsMySelected { get; internal set; }
        internal abstract void DoTask(ObservableCollection<RawData> completeRawFileListCollection, ObservableCollection<XMLdb> completeXmlDbList, ParamsObject po);

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
    }
}