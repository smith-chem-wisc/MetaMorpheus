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
using System.Xml;

namespace InternalLogicTaskLayer
{
    public class GPTMDTask : MyTaskEngine
    {
        public Tolerance precursorMassTolerance { get; set; }
        public List<ModListForGPTMD> listOfModListsForGPTMD { get; set; }

        public double tol { get; set; }

        public bool isotopeErrors { get; set; }
        public string outputFileName { get; set; }

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
            listOfModListsForGPTMD = new List<ModListForGPTMD>();
            foreach (var uu in modList)
                listOfModListsForGPTMD.Add(new ModListForGPTMD(uu));
            listOfModListsForGPTMD[0].Fixed = true;
            listOfModListsForGPTMD[1].Variable = true;
            listOfModListsForGPTMD[2].Localize = true;
            listOfModListsForGPTMD[3].Use = true;
            precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            this.taskType = MyTaskEnum.GPTMD;
        }

        private static bool ModFits(MorpheusModification attemptToLocalize, char v1, char prevAA, int peptideIndex, int peptideLength, int proteinIndex, int proteinLength)
        {
            if (!attemptToLocalize.AminoAcid.Equals('\0') && !attemptToLocalize.AminoAcid.Equals(v1))
                return false;
            if (!attemptToLocalize.PrevAminoAcid.Equals('\0') && !attemptToLocalize.PrevAminoAcid.Equals(prevAA))
                return false;
            if (attemptToLocalize.Type == ModificationType.ProteinNTerminus &&
                ((proteinIndex > 2) || (proteinIndex == 2 && prevAA != 'M')))
                return false;
            if (attemptToLocalize.Type == ModificationType.PeptideNTerminus && peptideIndex > 1)
                return false;
            if (attemptToLocalize.Type == ModificationType.PeptideCTerminus && peptideIndex < peptideLength)
                return false;
            if (attemptToLocalize.Type == ModificationType.ProteinCTerminus && proteinIndex < proteinLength)
                return false;
            return true;
        }

        private static IEnumerable<MorpheusModification> GetMod(double massDiff, bool isotopeErrors, IEnumerable<MorpheusModification> allMods, IEnumerable<Tuple<double, double>> combos, double tol)
        {
            foreach (var Mod in allMods)
            {
                if (Mod.MonoisotopicMassShift > massDiff - tol && Mod.MonoisotopicMassShift < massDiff + tol)
                    yield return Mod;
                if (isotopeErrors && Mod.MonoisotopicMassShift > massDiff - tol - 1.003 && Mod.MonoisotopicMassShift < massDiff + tol - 1.003)
                    yield return Mod;
                if (!double.IsNaN(Mod.AlternativeMassShift) && Mod.AlternativeMassShift > massDiff - tol && Mod.AlternativeMassShift < massDiff + tol)
                    yield return Mod;
                if (!double.IsNaN(Mod.AlternativeMassShift) && isotopeErrors && Mod.AlternativeMassShift > massDiff - tol - 1.003 && Mod.AlternativeMassShift < massDiff + tol - 1.003)
                    yield return Mod;
            }

            foreach (var combo in combos)
            {
                var m1 = combo.Item1;
                var m2 = combo.Item2;
                var combined = m1 + m2;
                if (combined > massDiff - tol && combined < massDiff + tol)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                }
                if (isotopeErrors && combined > massDiff - tol - 1.003 && combined < massDiff + tol - 1.003)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                }
            }
        }


        protected override MyResults RunSpecific()
        {
            outputFileName = Path.Combine(Path.GetDirectoryName(xmlDbFilenameList.First()), string.Join("-", xmlDbFilenameList.Select(b => Path.GetFileNameWithoutExtension(b))) + "GPTMD.xml");

            MyTaskResults myGPTMDresults = new MyTaskResults(this);
            myGPTMDresults.newDatabases = new List<string>();

            var currentRawFileList = rawDataFilenameList;

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();

            status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForGPTMD.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForGPTMD.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForGPTMD.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> gptmdModifications = listOfModListsForGPTMD.Where(b => b.Use).SelectMany(b => b.getMods()).ToList();
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
                output("Finished opening spectra file " + Path.GetFileName(origDataFile));

                ClassicSearchEngine searchEngine = new ClassicSearchEngine(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModes);

                ClassicSearchResults searchResults = (ClassicSearchResults)searchEngine.Run();

                output(searchResults.ToString());

                allPsms[0].AddRange(searchResults.outerPsms[0]);

                analysisEngine = new AnalysisEngine(searchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => Writing.WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => Writing.WriteToTabDelimitedTextFileWithDecoys(h, output_folder, "aggregate" + s), false);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                output(analysisResults.ToString());
            }

            if (currentRawFileList.Count > 1)
            {
                analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, null, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => Writing.WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => Writing.WriteToTabDelimitedTextFileWithDecoys(h, output_folder, "aggregate" + s), false);
                analysisResults = (AnalysisResults)analysisEngine.Run();
                output(analysisResults.ToString());
            }

            Dictionary<string, HashSet<Tuple<int, string>>> Mods = new Dictionary<string, HashSet<Tuple<int, string>>>();

            int modsAdded = 0;
            foreach (var ye in analysisResults.allResultingIdentifications[0].Where(b => b.QValue <= 0.01 && !b.isDecoy))
            {
                var baseSequence = ye.thisPSM.BaseSequence;
                double massDiff = ye.thisPSM.scanPrecursorMass - ye.thisPSM.PeptideMonoisotopicMass;
                foreach (MorpheusModification mod in GetMod(massDiff, isotopeErrors, gptmdModifications, combos, tol))
                {
                    foreach (var peptide in analysisResults.dict[ye.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications)])
                    {
                        int proteinLength = peptide.protein.Length;
                        var proteinAcession = peptide.protein.Accession;
                        for (int i = 0; i < baseSequence.Length; i++)
                        {
                            int indexInProtein = peptide.OneBasedStartResidueInProtein + i;

                            if (ModFits(mod, baseSequence[i], i > 0 ? baseSequence[i - 1] : peptide.PreviousAminoAcid, i + 1, baseSequence.Length, indexInProtein, proteinLength))
                            {
                                if (!Mods.ContainsKey(proteinAcession))
                                    Mods[proteinAcession] = new HashSet<Tuple<int, string>>();
                                var theTuple = new Tuple<int, string>(indexInProtein, mod.NameInXML);
                                if (!Mods[proteinAcession].Contains(theTuple))
                                {
                                    Mods[proteinAcession].Add(theTuple);
                                    modsAdded++;
                                }
                            }
                        }
                    }
                }
            }

            output("Modifications added = " + modsAdded);

            XmlWriterSettings xmlWriterSettings = new XmlWriterSettings()
            {
                Indent = true,
                IndentChars = "  "
            };

            output("Writing XML...");
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

            return myGPTMDresults;
        }

        private IEnumerable<Tuple<double, double>> LoadCombos()
        {
            yield return new Tuple<double, double>(15.994915, 15.994915);
        }

        public override void ValidateParams()
        {
            if (listOfModListsForGPTMD == null)
                throw new EngineValidationException("listOfModListsForGPTMD should not be null");
            if (listOfModListsForGPTMD.Where(b => b.Use).Count() == 0)
                throw new EngineValidationException("Need to marks some modification files for use in GPTMD");
        }
    }
}