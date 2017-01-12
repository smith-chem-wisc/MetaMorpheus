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
            this.taskType = MyTaskEnum.GPTMD;
            tol = 0.003;
        }

        public bool isotopeErrors { get; set; }
        public List<ModListForGPTMDTask> listOfModListsForGPTMD { get; set; }
        public string outputFileName { get; set; }
        public Tolerance precursorMassTolerance { get; set; }
        public double tol { get; set; }

        protected override void ValidateParams()
        {
            if (listOfModListsForGPTMD == null)
                throw new EngineValidationException("listOfModListsForGPTMD should not be null");
            if (listOfModListsForGPTMD.Where(b => b.GPTMD).Count() == 0)
                throw new EngineValidationException("Need to marks some modification files for use in GPTMD");
        }

        internal override string GetSpecificTaskInfo()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("isotopeErrors: " + isotopeErrors.ToString());
            sb.AppendLine("Fixed mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Fixed).Select(b => b.FileName)));
            sb.AppendLine("Variable mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Variable).Select(b => b.FileName)));
            sb.AppendLine("Localized mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.Localize).Select(b => b.FileName)));
            sb.AppendLine("GPTMD mod lists: " + string.Join(",", listOfModListsForGPTMD.Where(b => b.GPTMD).Select(b => b.FileName)));
            sb.AppendLine("outputFileName: " + outputFileName);
            sb.AppendLine("precursorMassTolerance: " + precursorMassTolerance);
            sb.Append("tol: " + tol);
            return sb.ToString();
        }

        protected override MyResults RunSpecific()
        {
            outputFileName = Path.Combine(Path.GetDirectoryName(xmlDbFilenameList.First()), string.Join("-", xmlDbFilenameList.Select(b => Path.GetFileNameWithoutExtension(b))) + "GPTMD.xml");

            MyTaskResults myGPTMDresults = new MyGPTMDTaskResults(this);
            myGPTMDresults.newDatabases = new List<string>();

            var currentRawFileList = rawDataFilenameList;

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();

            status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForGPTMD.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForGPTMD.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForGPTMD.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> gptmdModifications = listOfModListsForGPTMD.Where(b => b.GPTMD).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            GenerateModsFromStrings(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            IEnumerable<Tuple<double, double>> combos = LoadCombos();

            SearchMode searchMode = new DotSearchMode("", gptmdModifications, combos, precursorMassTolerance);
            List<SearchMode> searchModes = new List<SearchMode>() { searchMode };

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
                    myMsDataFile = new Mzml(origDataFile, 400);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, 400);
                status("Opening spectra file...");
                myMsDataFile.Open();
                //output("Finished opening spectra file " + Path.GetFileName(origDataFile));

                ClassicSearchEngine searchEngine = new ClassicSearchEngine(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModes);

                ClassicSearchResults searchResults = (ClassicSearchResults)searchEngine.Run();

                //output(searchResults.ToString());

                allPsms[0].AddRange(searchResults.outerPsms[0]);

                analysisEngine = new AnalysisEngine(searchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => WriteToTabDelimitedTextFileWithDecoys(h, output_folder, "aggregate" + s), false);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            if (currentRawFileList.Count > 1)
            {
                analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, null, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => WriteToTabDelimitedTextFileWithDecoys(h, output_folder, "aggregate" + s), false);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                //output(analysisResults.ToString());
            }

            GPTMDEngine gptmdEngine = new GPTMDEngine(analysisResults.allResultingIdentifications, analysisResults.dict, variableModifications, localizeableModifications, isotopeErrors, gptmdModifications, combos, tol);
            GPTMDResults gptmdResults = (GPTMDResults)gptmdEngine.Run();

            //output(gptmdResults.ToString());

            WriteGPTMDdatabse(gptmdResults.mods, proteinList);

            myGPTMDresults.newDatabases.Add(outputFileName);
            return myGPTMDresults;
        }

        private IEnumerable<Tuple<double, double>> LoadCombos()
        {
            yield return new Tuple<double, double>(15.994915, 15.994915);
        }

        private void WriteGPTMDdatabse(Dictionary<string, HashSet<Tuple<int, string>>> Mods, List<Protein> proteinList)
        {
            XmlWriterSettings xmlWriterSettings = new XmlWriterSettings()
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
        }
    }
}