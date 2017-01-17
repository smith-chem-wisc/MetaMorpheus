using InternalLogicEngineLayer;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text;
using System.Xml;

namespace InternalLogicTaskLayer
{
    public class GPTMDTask : MyTaskEngine
    {

        #region Public Fields

        public List<ModListForGPTMDTask> listOfModListsForGPTMD;
        public Tolerance precursorMassTolerance;

        #endregion Public Fields

        #region Public Constructors

        public GPTMDTask(ObservableCollection<ModList> modList)
        {
            // Set default values here:
            maxMissedCleavages = 2;
            protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            maxModificationIsoforms = 4096;
            initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            bIons = true;
            yIons = true;
            listOfModListsForGPTMD = new List<ModListForGPTMDTask>();
            foreach (var uu in modList)
                listOfModListsForGPTMD.Add(new ModListForGPTMDTask(uu));
            listOfModListsForGPTMD[0].Fixed = true;
            listOfModListsForGPTMD[1].Variable = true;
            listOfModListsForGPTMD[2].Localize = true;
            listOfModListsForGPTMD[3].GPTMD = true;
            precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            taskType = MyTaskEnum.GPTMD;
            tol = 0.003;
            isotopeErrors = false;
            maxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance productMassTolerance { get; set; }
        public double tol { get; set; }
        public bool isotopeErrors { get; set; }

        #endregion Public Properties

        #region Protected Methods

        protected override string GetSpecificTaskInfo()
        {
            var sb = new StringBuilder();
            sb.AppendLine("isotopeErrors: " + isotopeErrors);
            sb.AppendLine("Fixed mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Fixed).Select(b => b.FileName)));
            sb.AppendLine("Variable mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Variable).Select(b => b.FileName)));
            sb.AppendLine("Localized mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Localize).Select(b => b.FileName)));
            sb.AppendLine("GPTMD mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.GPTMD).Select(b => b.FileName)));
            sb.AppendLine("precursorMassTolerance: " + precursorMassTolerance);
            sb.AppendLine("productMassTolerance: " + productMassTolerance);
            sb.Append("tol: " + tol);
            return sb.ToString();
        }

        protected override MyResults RunSpecific()
        {
            string outputXMLdbFullName = Path.Combine(output_folder, string.Join("-", xmlDbFilenameList.Select(b => Path.GetFileNameWithoutExtension(b))) + "GPTMD.xml");

            MyTaskResults myGPTMDresults = new MyGPTMDTaskResults(this);
            myGPTMDresults.newDatabases = new List<string>();

            var currentRawFileList = rawDataFilenameList;

            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForGPTMD.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForGPTMD.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForGPTMD.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> gptmdModifications = listOfModListsForGPTMD.Where(b => b.GPTMD).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            MatchXMLmodsToKnownMods(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            IEnumerable<Tuple<double, double>> combos = LoadCombos();

            SearchMode searchMode = new DotSearchMode("", gptmdModifications.Select(b => b.MonoisotopicMassShift).Concat(combos.Select(b => b.Item1 + b.Item2)).OrderBy(b => b), precursorMassTolerance);
            var searchModes = new List<SearchMode> { searchMode };

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[1];
            allPsms[0] = new List<ParentSpectrumMatch>();

            status("Loading proteins...");
            var proteinList = xmlDbFilenameList.SelectMany(b => getProteins(true, identifiedModsInXML, b)).ToList();
            AnalysisEngine analysisEngine;
            AnalysisResults analysisResults = null;
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                status("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, maxNumPeaksPerScan);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, maxNumPeaksPerScan);
                status("Opening spectra file...");
                myMsDataFile.Open();

                var listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMs2Scan(b)).OrderBy(b => b.precursorMass).ToArray();

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, spectraFileIndex, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maxMissedCleavages, maxModificationIsoforms);

                var searchResults = (ClassicSearchResults)searchEngine.Run();

                allPsms[0].AddRange(searchResults.outerPsms[0]);

                analysisEngine = new AnalysisEngine(searchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, "aggregate" + s), null, false, maxMissedCleavages, maxModificationIsoforms);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            if (currentRawFileList.Count > 1)
            {
                analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, null, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, "aggregate" + s), null, false, maxMissedCleavages, maxModificationIsoforms);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            var gptmdEngine = new GPTMDEngine(analysisResults.allResultingIdentifications[0], isotopeErrors, gptmdModifications, combos, tol);
            var gptmdResults = (GPTMDResults)gptmdEngine.Run();

            //output(gptmdResults.ToString());

            WriteGPTMDdatabse(gptmdResults.mods, proteinList.Where(b => !b.isDecoy).ToList(), outputXMLdbFullName);

            myGPTMDresults.newDatabases.Add(outputXMLdbFullName);

            return myGPTMDresults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<Tuple<double, double>> LoadCombos()
        {
            yield return new Tuple<double, double>(15.994915, 15.994915);
        }

        private void WriteGPTMDdatabse(Dictionary<string, HashSet<Tuple<int, string>>> Mods, List<Protein> proteinList, string outputFileName)
        {
            var xmlWriterSettings = new XmlWriterSettings
            {
                Indent = true,
                IndentChars = "  "
            };

            status("Writing XML...");
            using (XmlWriter writer = XmlWriter.Create(outputFileName, xmlWriterSettings))
            {
                writer.WriteStartDocument();
                writer.WriteStartElement("uniprot");

                foreach (Protein protein in proteinList)
                {
                    writer.WriteStartElement("entry");
                    writer.WriteAttributeString("dataset", protein.dataset_abbreviation);
                    writer.WriteStartElement("accession");
                    writer.WriteString(protein.Accession);
                    writer.WriteEndElement();
                    writer.WriteStartElement("name");
                    writer.WriteString(protein.name);
                    writer.WriteEndElement();

                    writer.WriteStartElement("protein");
                    writer.WriteStartElement("recommendedName");
                    writer.WriteStartElement("fullName");
                    writer.WriteString(protein.fullName);
                    writer.WriteEndElement();
                    writer.WriteEndElement();
                    writer.WriteEndElement();

                    for (int i = 0; i < protein.bigPeptideTypes.Count(); i++)
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", protein.bigPeptideTypes[i]);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("begin");
                        writer.WriteAttributeString("position", protein.oneBasedBeginPositions[i].ToString());
                        writer.WriteEndElement();
                        writer.WriteStartElement("end");
                        writer.WriteAttributeString("position", protein.oneBasedEndPositions[i].ToString());
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    IEnumerable<Tuple<int, string>> SortedMods = protein.OneBasedPossibleLocalizedModifications.SelectMany(
                        b => b.Value.Select(c => new Tuple<int, string>(b.Key, c.NameInXML)
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
                        writer.WriteAttributeString("position", ye.Item1.ToString());
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    writer.WriteStartElement("sequence");
                    writer.WriteAttributeString("length", protein.Length.ToString());
                    writer.WriteString(protein.BaseSequence);
                    writer.WriteEndElement();

                    writer.WriteEndElement();
                }

                writer.WriteEndElement();
                writer.WriteEndDocument();
            }
            SucessfullyFinishedWritingFile(outputFileName);
        }

        #endregion Private Methods

    }
}