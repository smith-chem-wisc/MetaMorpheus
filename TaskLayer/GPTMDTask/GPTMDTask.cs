using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using EngineLayer.Gptmd;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Xml;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class GptmdTask : MetaMorpheusTask
    {

        #region Private Fields

        private const double binTolInDaltons = 0.003;

        #endregion Private Fields

        #region Public Constructors

        public GptmdTask()
        {
            // Set default values here:
            MaxMissedCleavages = 2;
            Protease = ProteaseDictionary.Instance["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            PrecursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            BIons = true;
            YIons = true;

            ListOfModListsFixed = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("f.txt")) };
            ListOfModListsVariable = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("v.txt")) };
            ListOfModListsLocalize = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("ptmlist.txt")) };

            ListOfModListsGptmd = new List<ModList> {
                AllModLists.First(b => b.FileName.EndsWith("m.txt")),
                AllModLists.First(b => b.FileName.EndsWith("glyco.txt")),
            };

            TaskType = MyTask.Gptmd;
            IsotopeErrors = false;
            MaxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModList> ListOfModListsFixed { get; set; }
        public List<ModList> ListOfModListsVariable { get; set; }
        public List<ModList> ListOfModListsLocalize { get; set; }
        public List<ModList> ListOfModListsGptmd { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }
        public bool IsotopeErrors { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string SpecificTaskInfo
        {
            get
            {
                var sb = new StringBuilder();
                sb.AppendLine("isotopeErrors: " + IsotopeErrors);
                sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsFixed.Select(b => b.FileName)));
                sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsVariable.Select(b => b.FileName)));
                sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsLocalize.Select(b => b.FileName)));
                sb.AppendLine("GPTMD mod lists: " + string.Join(",", ListOfModListsGptmd.Select(b => b.FileName)));
                sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
                sb.Append("PrecursorMassTolerance: " + PrecursorMassTolerance);
                return sb.ToString();
            }
        }

        #endregion Protected Properties

        #region Public Methods

        public static void WriteXmlDatabase(Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> Mods, List<Protein> proteinList, string outputFileName)
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
                        writer.WriteAttributeString("position", protein.OneBasedBeginPositions[i].ToString());
                        writer.WriteEndElement();
                        writer.WriteStartElement("end");
                        writer.WriteAttributeString("position", protein.OneBasedEndPositions[i].ToString());
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }
                    foreach (var ye in protein.OneBasedPossibleLocalizedModifications.OrderBy(b => b.Key))
                    {
                        foreach (var nice in ye.Value)
                        {
                            writer.WriteStartElement("feature");
                            writer.WriteAttributeString("type", "modified residue");
                            writer.WriteAttributeString("description", nice.id);
                            //writer.WriteStartElement("db");
                            //writer.WriteString(ye.Item3);
                            //writer.WriteEndElement();
                            writer.WriteStartElement("location");
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", ye.Key.ToString(CultureInfo.InvariantCulture));
                            writer.WriteEndElement();
                            writer.WriteEndElement();
                            writer.WriteEndElement();

                        }
                    }
                    if (Mods.ContainsKey(protein.Accession))
                        foreach (var ye in Mods[protein.Accession].OrderBy(b => b.Item1))
                        {
                            writer.WriteStartElement("feature");
                            writer.WriteAttributeString("type", "modified residue");
                            writer.WriteAttributeString("description", ye.Item2.id);
                            //writer.WriteStartElement("db");
                            //writer.WriteString(ye.Item3);
                            //writer.WriteEndElement();
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
            myGPTMDresults.newDatabases = new List<DbForTask>();

            var currentRawFileList = rawDataFilenameList;

            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...");
            List<ModificationWithMass> variableModifications = ListOfModListsVariable.SelectMany(b => b.Mods).Where(b => b is ModificationWithMass).Select(b => b as ModificationWithMass).ToList();
            List<ModificationWithMass> fixedModifications = ListOfModListsFixed.SelectMany(b => b.Mods).Where(b => b is ModificationWithMass).Select(b => b as ModificationWithMass).ToList();
            List<ModificationWithMass> localizeableModifications = ListOfModListsLocalize.SelectMany(b => b.Mods).Where(b => b is ModificationWithMass).Select(b => b as ModificationWithMass).ToList();
            List<ModificationWithMass> gptmdModifications = ListOfModListsGptmd.SelectMany(b => b.Mods).Where(b => b is ModificationWithMass).Select(b => b as ModificationWithMass).ToList();

            IEnumerable<Tuple<double, double>> combos = LoadCombos().ToList();

            // Do not remove the zero!!! It's needed here
            SearchMode searchMode = new DotSearchMode("", gptmdModifications.SelectMany(b => b.massesObserved).Concat(combos.Select(b => b.Item1 + b.Item2)).Concat(new List<double> { 0 }).OrderBy(b => b), PrecursorMassTolerance);
            var searchModes = new List<SearchMode> { searchMode };

            List<PsmParent>[] allPsms = new List<PsmParent>[1];
            allPsms[0] = new List<PsmParent>();

            List<ProductType> lp = new List<ProductType>();
            if (BIons)
                lp.Add(ProductType.B);
            if (YIons)
                lp.Add(ProductType.Y);

            Status("Loading proteins...");
            var allKnownModifications = GetDict(localizeableModifications);
            var proteinList = dbFilenameList.SelectMany(b => ProteinDbLoader.LoadProteinDb(b.FileName, true, allKnownModifications, b.IsContaminant)).ToList();
            AnalysisEngine analysisEngine;
            AnalysisResults analysisResults = null;
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                Status("Loading spectra file...");
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
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

                analysisEngine = new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), null, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, binTolInDaltons);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            if (currentRawFileList.Count > 1)
            {
                analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s), null, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, binTolInDaltons);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            var gptmdEngine = new GptmdEngine(analysisResults.AllResultingIdentifications[0], IsotopeErrors, gptmdModifications, combos, PrecursorMassTolerance);
            var gptmdResults = (GptmdResults)gptmdEngine.Run();

            //output(gptmdResults.ToString());

            string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

            WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy).ToList(), outputXMLdbFullName);

            SucessfullyFinishedWritingFile(outputXMLdbFullName);

            // TODO: Fix so not always outputting a contaminant
            myGPTMDresults.newDatabases.Add(new DbForTask(outputXMLdbFullName, false));

            return myGPTMDresults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<Tuple<double, double>> LoadCombos()
        {
            using (StreamReader r = new StreamReader(Path.Combine("Data", @"combos.txt")))
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