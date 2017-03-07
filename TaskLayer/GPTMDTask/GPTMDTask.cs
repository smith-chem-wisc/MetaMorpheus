using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using EngineLayer.Gptmd;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
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
            CIons = false;
            ZdotIons = false;

            ListOfModListsFixed = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("f.txt")) };
            ListOfModListsVariable = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("v.txt")) };
            ListOfModListsLocalize = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("ptmlist.txt")) };

            ListOfModListsGptmd = new List<ModList> {
                AllModLists.First(b => b.FileName.EndsWith("m.txt")),
                AllModLists.First(b => b.FileName.EndsWith("glyco.txt")),
                AllModLists.First(b => b.FileName.EndsWith("metals.txt")),
                AllModLists.First(b => b.FileName.EndsWith("pt.txt"))
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

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawFileList, string taskId)
        {
            MyTaskResults myGPTMDresults = new MyGptmdTaskResults(this);
            myGPTMDresults.newDatabases = new List<DbForTask>();

            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = ListOfModListsVariable.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> fixedModifications = ListOfModListsFixed.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> localizeableModifications = ListOfModListsLocalize.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> gptmdModifications = ListOfModListsGptmd.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();

            IEnumerable<Tuple<double, double>> combos = LoadCombos().ToList();

            // Do not remove the zero!!! It's needed here
            SearchMode searchMode = new DotSearchMode("", gptmdModifications.SelectMany(b => b.massesObserved).Concat(combos.Select(b => b.Item1 + b.Item2)).Concat(new List<double> { 0 }).GroupBy(b => Math.Round(b, 6)).Select(b => b.FirstOrDefault()).OrderBy(b => b), PrecursorMassTolerance);
            var searchModes = new List<SearchMode> { searchMode };

            List<PsmParent>[] allPsms = new List<PsmParent>[1];
            allPsms[0] = new List<PsmParent>();

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            List<ProductType> lp = new List<ProductType>();
            if (BIons)
                lp.Add(ProductType.B);
            if (YIons)
                lp.Add(ProductType.Y);
            if (CIons)
                lp.Add(ProductType.C);
            if (ZdotIons)
                lp.Add(ProductType.Zdot);

            Status("Loading proteins...", new List<string> { taskId });
            Dictionary<string, Modification> um = null;
            var proteinList = currentXmlDbFilenameList.SelectMany(b => LoadProteinDb(b.FileName, true, localizeableModifications, b.IsContaminant, out um)).ToList();

            AnalysisResults analysisResults = null;
            var numRawFiles = currentRawFileList.Count;
            for (int spectraFileIndex = 0; spectraFileIndex < numRawFiles; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                Status("Loading spectra file...", new List<string> { taskId });
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                else
                    myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                Status("Opening spectra file...", new List<string> { taskId });

                var searchResults = (ClassicSearchResults)new ClassicSearchEngine(MyEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray(), myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, origDataFile, lp, new List<string> { taskId }).Run();
                myGPTMDresults.AddResultText(searchResults);

                allPsms[0].AddRange(searchResults.OuterPsms[0]);

                analysisResults = (AnalysisResults)new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), null, false, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, binTolInDaltons, initiatorMethionineBehavior, false, 0, 0).Run();
                analysisResults = (AnalysisResults)new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, origDataFile }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, origDataFile }), null, false, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId }).Run();
                myGPTMDresults.AddResultText(analysisResults);
            }

            if (numRawFiles > 1)
            {
                analysisResults = (AnalysisResults)new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s), null, false, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, binTolInDaltons, initiatorMethionineBehavior, false, 0, 0).Run();
                analysisResults = (AnalysisResults)new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s, new List<string> { taskId }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s, new List<string> { taskId }), null, false, false, MaxMissedCleavages, MaxModificationIsoforms, true, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId }).Run();
                myGPTMDresults.AddResultText(analysisResults);
            }

            var gptmdResults = (GptmdResults)new GptmdEngine(analysisResults.AllResultingIdentifications[0], IsotopeErrors, gptmdModifications, combos, PrecursorMassTolerance).Run();
            myGPTMDresults.AddResultText(gptmdResults);

            if (currentXmlDbFilenameList.Any(b => !b.IsContaminant))
            {
                string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", currentXmlDbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

                ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });

                myGPTMDresults.newDatabases.Add(new DbForTask(outputXMLdbFullName, false));
            }
            if (currentXmlDbFilenameList.Any(b => b.IsContaminant))
            {
                string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", currentXmlDbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

                ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });

                myGPTMDresults.newDatabases.Add(new DbForTask(outputXMLdbFullNameContaminants, true));
            }
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