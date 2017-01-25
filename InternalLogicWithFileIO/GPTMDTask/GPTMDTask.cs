using InternalLogicEngineLayer;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Xml;

namespace InternalLogicTaskLayer
{
    public class GptmdTask : MyTaskEngine
    {

        #region Public Fields

        public List<ModListForGPTMDTask> listOfModListsForGPTMD;

        #endregion Public Fields

        #region Public Constructors

        public GptmdTask(ObservableCollection<ModList> modList)
        {
            // Set default values here:
            MaxMissedCleavages = 2;
            Protease = ProteaseDictionary.Instance["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            BIons = true;
            YIons = true;
            listOfModListsForGPTMD = new List<ModListForGPTMDTask>();
            foreach (var uu in modList)
                listOfModListsForGPTMD.Add(new ModListForGPTMDTask(uu));
            listOfModListsForGPTMD[0].Fixed = true;
            listOfModListsForGPTMD[1].Variable = true;
            listOfModListsForGPTMD[2].Localize = true;
            listOfModListsForGPTMD[3].Gptmd = true;
            listOfModListsForGPTMD[4].Gptmd = true;
            TaskType = MyTask.Gptmd;
            TolInDaltons = 0.01;
            IsotopeErrors = false;
            MaxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance ProductMassTolerance { get; set; }
        public double TolInDaltons { get; set; }
        public bool IsotopeErrors { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string SpecificTaskInfo
        {
            get
            {
                var sb = new StringBuilder();
                sb.AppendLine("isotopeErrors: " + IsotopeErrors);
                sb.AppendLine("Fixed mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Fixed).Select(b => b.FileName)));
                sb.AppendLine("Variable mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Variable).Select(b => b.FileName)));
                sb.AppendLine("Localized mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Localize).Select(b => b.FileName)));
                sb.AppendLine("GPTMD mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Gptmd).Select(b => b.FileName)));
                sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
                sb.Append("TolInDaltons: " + TolInDaltons);
                return sb.ToString();
            }
        }

        #endregion Protected Properties

        #region Public Methods

        public static void WriteXmlDatabase(Dictionary<string, HashSet<Tuple<int, string>>> Mods, List<Protein> proteinList, string outputFileName)
        {
            var xmlWriterSettings = new XmlWriterSettings
            {
                Indent = true,
                IndentChars = "  "
            };

            using (XmlWriter writer = XmlWriter.Create(outputFileName, xmlWriterSettings))
            {
                writer.WriteStartDocument();
                writer.WriteStartElement("uniprot");

                foreach (Protein protein in proteinList)
                {
                    writer.WriteStartElement("entry");
                    writer.WriteStartElement("accession");
                    writer.WriteString(protein.Accession);
                    writer.WriteEndElement();
                    writer.WriteStartElement("name");
                    writer.WriteString(protein.Name);
                    writer.WriteEndElement();

                    writer.WriteStartElement("protein");
                    writer.WriteStartElement("recommendedName");
                    writer.WriteStartElement("fullName");
                    writer.WriteString(protein.FullName);
                    writer.WriteEndElement();
                    writer.WriteEndElement();
                    writer.WriteEndElement();

                    for (int i = 0; i < protein.BigPeptideTypes.Count(); i++)
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", protein.BigPeptideTypes[i]);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("begin");
                        writer.WriteAttributeString("position", protein.OneBasedBeginPositions[i].ToString(CultureInfo.InvariantCulture));
                        writer.WriteEndElement();
                        writer.WriteStartElement("end");
                        writer.WriteAttributeString("position", protein.OneBasedEndPositions[i].ToString(CultureInfo.InvariantCulture));
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    IEnumerable<Tuple<int, string>> SortedMods = protein.OneBasedPossibleLocalizedModifications.SelectMany(
                        b => b.Value.Select(c => new Tuple<int, string>(b.Key, c.NameInXml)
                        ));
                    IEnumerable<Tuple<int, string>> FinalSortedMods;
                    if (Mods.ContainsKey(protein.Accession))
                        FinalSortedMods = SortedMods.Union(Mods[protein.Accession]).OrderBy(b => b.Item1);
                    else
                        FinalSortedMods = SortedMods.OrderBy(b => b.Item1);
                    foreach (var ye in FinalSortedMods)
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", "modified residue");
                        writer.WriteAttributeString("description", ye.Item2);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("position");
                        writer.WriteAttributeString("position", ye.Item1.ToString(CultureInfo.InvariantCulture));
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    writer.WriteStartElement("sequence");
                    writer.WriteAttributeString("length", protein.Length.ToString(CultureInfo.InvariantCulture));
                    writer.WriteString(protein.BaseSequence);
                    writer.WriteEndElement();

                    writer.WriteEndElement();
                }

                writer.WriteEndElement();
                writer.WriteEndDocument();
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            MyTaskResults myGPTMDresults = new MyGPTMDTaskResults(this);
            myGPTMDresults.newDatabases = new List<XmlForTask>();

            var currentRawFileList = rawDataFilenameList;

            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForGPTMD.Where(b => b.Variable).SelectMany(b => b.Mods).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForGPTMD.Where(b => b.Fixed).SelectMany(b => b.Mods).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForGPTMD.Where(b => b.Localize).SelectMany(b => b.Mods).ToList();
            List<MorpheusModification> gptmdModifications = listOfModListsForGPTMD.Where(b => b.Gptmd).SelectMany(b => b.Mods).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            MatchXMLmodsToKnownMods(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            IEnumerable<Tuple<double, double>> combos = LoadCombos().ToList();

            // Do not remove the zero!!! It's needed here
            SearchMode searchMode = new DotSearchMode("", gptmdModifications.Select(b => b.MonoisotopicMassShift).Concat(combos.Select(b => b.Item1 + b.Item2)).Concat(new List<double> { 0 }).OrderBy(b => b), new Tolerance(ToleranceUnit.Absolute, TolInDaltons));
            var searchModes = new List<SearchMode> { searchMode };

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[1];
            allPsms[0] = new List<ParentSpectrumMatch>();

            List<ProductType> lp = new List<ProductType>();
            if (BIons)
                lp.Add(ProductType.B);
            if (YIons)
                lp.Add(ProductType.Y);

            Status("Loading proteins...");
            var proteinList = xmlDbFilenameList.SelectMany(b => GetProteins(true, identifiedModsInXML, b)).ToList();
            AnalysisEngine analysisEngine;
            AnalysisResults analysisResults = null;
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                Status("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, MaxNumPeaksPerScan);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, MaxNumPeaksPerScan);
                Status("Opening spectra file...");
                myMsDataFile.Open();

                var listOfSortedms2Scans = GetMs2Scans(myMsDataFile).OrderBy(b => b.PrecursorMass).ToArray();

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, myMsDataFile.Name, lp);

                var searchResults = (ClassicSearchResults)searchEngine.Run();

                allPsms[0].AddRange(searchResults.OuterPsms[0]);

                analysisEngine = new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), null, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, TolInDaltons);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            if (currentRawFileList.Count > 1)
            {
                analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s), null, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, TolInDaltons);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            var gptmdEngine = new GptmdEngine(analysisResults.AllResultingIdentifications[0], IsotopeErrors, gptmdModifications, combos, TolInDaltons);
            var gptmdResults = (GptmdResults)gptmdEngine.Run();

            //output(gptmdResults.ToString());

            string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", xmlDbFilenameList.Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

            WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy).ToList(), outputXMLdbFullName);

            SucessfullyFinishedWritingFile(outputXMLdbFullName);

            // TODO: Fix so not always outputting a contaminant
            myGPTMDresults.newDatabases.Add(new XmlForTask(outputXMLdbFullName, false));

            return myGPTMDresults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<Tuple<double, double>> LoadCombos()
        {
            using (StreamReader r = new StreamReader(@"combos.txt"))
            {
                while (r.Peek() >= 0)
                {
                    var line = r.ReadLine().Split(' ');
                    yield return new Tuple<double, double>(double.Parse(line[0]), double.Parse(line[1]));
                }
            }
        }

        #endregion Private Methods

    }
}