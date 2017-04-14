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

        public GptmdTask() : base(MyTask.Gptmd)
        {
            // Set default values here:
            MaxMissedCleavages = 2;
            MinPeptideLength = null;
            MaxPeptideLength = null;
            Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            PrecursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 2);
            BIons = true;
            YIons = true;
            CIons = false;
            ZdotIons = false;

            LocalizeAll = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();
            ListOfModsGptmd = GlobalTaskLevelSettings.AllModsKnown.Where(b => b.modificationType.Equals("metals")).Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            IsotopeErrors = false;

            DatabaseReferencesToKeep = new string[] { "" };
        }

        #endregion Public Constructors

        #region Public Properties

        public static List<string> AllModLists { get; private set; }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }

        public int? MinPeptideLength { get; set; }

        public int? MaxPeptideLength { get; set; }

        public int MaxModificationIsoforms { get; set; }

        public Protease Protease { get; set; }

        public bool BIons { get; set; }

        public bool YIons { get; set; }

        public bool ZdotIons { get; set; }

        public bool CIons { get; set; }
        public List<Tuple<string, string>> ListOfModsFixed { get; set; }
        public List<Tuple<string, string>> ListOfModsVariable { get; set; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; set; }
        public List<Tuple<string, string>> ListOfModsGptmd { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }
        public string[] DatabaseReferencesToKeep { get; set; }
        public bool IsotopeErrors { get; set; }
        public bool LocalizeAll { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(TaskType.ToString());
            sb.AppendLine("The initiator methionine behavior is set to "
                + InitiatorMethionineBehavior
                + " and the maximum number of allowed missed cleavages is "
                + MaxMissedCleavages);
            sb.AppendLine("MinPeptideLength: " + MinPeptideLength);
            sb.AppendLine("MaxPeptideLength: " + MaxPeptideLength);
            sb.AppendLine("maxModificationIsoforms: " + MaxModificationIsoforms);
            sb.AppendLine("protease: " + Protease);
            sb.AppendLine("bIons: " + BIons);
            sb.AppendLine("yIons: " + YIons);
            sb.AppendLine("cIons: " + CIons);
            sb.AppendLine("zdotIons: " + ZdotIons);
            sb.AppendLine("isotopeErrors: " + IsotopeErrors);
            //sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsFixed));
            //sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsVariable));
            //sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsLocalize));
            //sb.AppendLine("GPTMD mod lists: " + string.Join(",", ListOfModListsGptmd));
            sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
            sb.Append("PrecursorMassTolerance: " + PrecursorMassTolerance);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawFileList, string taskId)
        {
            myTaskResults = new MyTaskResults(this)
            {
                newDatabases = new List<DbForTask>()
            };
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", new List<string> { taskId });

            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            List<ModificationWithMass> localizeableModifications;
            if (LocalizeAll)
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            List<ModificationWithMass> gptmdModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsGptmd.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                if (!modsDictionary.ContainsKey(mod))
                    modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            IEnumerable<Tuple<double, double>> combos = LoadCombos(gptmdModifications).ToList();

            // Daltons around zero = 3.5
            SearchMode searchMode;
            if (IsotopeErrors)
                searchMode = new DaltonsAroundZeroWithInnerSearchMode(new DotSearchMode("", gptmdModifications.Select(b => b.monoisotopicMass).Concat(gptmdModifications.Select(b => b.monoisotopicMass + 1.003)).Concat(GetObservedMasses(variableModifications.Concat(fixedModifications), gptmdModifications)).Concat(combos.Select(b => b.Item1 + b.Item2)).GroupBy(b => Math.Round(b, 6)).Select(b => b.FirstOrDefault()).OrderBy(b => b), PrecursorMassTolerance), 3.5);
            else
                searchMode = new DaltonsAroundZeroWithInnerSearchMode(new DotSearchMode("", gptmdModifications.Select(b => b.monoisotopicMass).Concat(GetObservedMasses(variableModifications.Concat(fixedModifications), gptmdModifications)).Concat(combos.Select(b => b.Item1 + b.Item2)).GroupBy(b => Math.Round(b, 6)).Select(b => b.FirstOrDefault()).OrderBy(b => b), PrecursorMassTolerance), 3.5);

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
            var proteinList = currentXmlDbFilenameList.SelectMany(b => LoadProteinDb(b.FileName, true, localizeableModifications, b.IsContaminant, DatabaseReferencesToKeep, out um)).ToList();

            AnalysisResults analysisResults = null;
            var numRawFiles = currentRawFileList.Count;

            Status("Running G-PTM-D...", new List<string> { taskId });

            for (int spectraFileIndex = 0; spectraFileIndex < numRawFiles; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];

                NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Searches", origDataFile });
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Searches", origDataFile });

                Status("Loading spectra file...", new List<string> { taskId, "Individual Searches", origDataFile });
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                else
                    myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                Status("Opening spectra file...", new List<string> { taskId, "Individual Searches", origDataFile });

                var searchResults = (ClassicSearchResults)new ClassicSearchEngine(MetaMorpheusEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray(), myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, origDataFile, lp, new List<string> { taskId, "Individual Searches", origDataFile }, false).Run();

                allPsms[0].AddRange(searchResults.OuterPsms[0]);

                analysisResults = (AnalysisResults)new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s, ss), null, false, false, false, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, true, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }, false, 0, 0, modsDictionary).Run();

                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Searches", origDataFile });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Searches", origDataFile }));
            }
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Searches" }));

            if (numRawFiles > 1)
            {
                analysisResults = (AnalysisResults)new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + "_" + s, new List<string> { taskId }), (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, "aggregate" + "_" + s, ss), null, false, false, false, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, true, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId }, false, 0, 0, modsDictionary).Run();
            }

            var gptmdResults = (GptmdResults)new GptmdEngine(analysisResults.AllResultingIdentifications[0], IsotopeErrors, gptmdModifications, combos, PrecursorMassTolerance).Run();

            if (currentXmlDbFilenameList.Any(b => !b.IsContaminant))
            {
                string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", currentXmlDbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

                ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });

                myTaskResults.newDatabases.Add(new DbForTask(outputXMLdbFullName, false));
            }
            if (currentXmlDbFilenameList.Any(b => b.IsContaminant))
            {
                string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", currentXmlDbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

                ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });

                myTaskResults.newDatabases.Add(new DbForTask(outputXMLdbFullNameContaminants, true));
            }
            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<double> GetObservedMasses(IEnumerable<ModificationWithMass> enumerable, List<ModificationWithMass> gptmdModifications)
        {
            foreach (var modOnPeptide in enumerable)
            {
                foreach (var modToLocalize in gptmdModifications)
                {
                    if (modOnPeptide.motif.Motif.Equals(modToLocalize.motif.Motif))
                    {
                        yield return modToLocalize.monoisotopicMass - modOnPeptide.monoisotopicMass;
                    }
                }
            }
        }

        private IEnumerable<Tuple<double, double>> LoadCombos(List<ModificationWithMass> allowedCombos)
        {
            using (StreamReader r = new StreamReader(Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Data", @"combos.txt")))
            {
                while (r.Peek() >= 0)
                {
                    var line = r.ReadLine().Split(' ');
                    var mass1 = double.Parse(line[0]);
                    var mass2 = double.Parse(line[1]);
                    if (allowedCombos.Any(b => b.monoisotopicMass - mass1 < 1e-3) &&
                        allowedCombos.Any(b => b.monoisotopicMass - mass2 < 1e-3))
                        yield return new Tuple<double, double>(mass1, mass2);
                }
            }
        }

        #endregion Private Methods

    }
}