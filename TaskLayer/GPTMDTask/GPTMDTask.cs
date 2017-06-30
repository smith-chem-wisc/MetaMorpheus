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
            ListOfModsLocalize = new List<Tuple<string, string>>();
            ListOfModsGptmd = GlobalTaskLevelSettings.AllModsKnown.Where(b =>
            b.modificationType.Equals("Glycan") ||
            b.modificationType.Equals("Mod") ||
            b.modificationType.Equals("PeptideTermMod") ||
            b.modificationType.Equals("Metal") ||
            b.modificationType.Equals("ProteinTermMod")).Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();
            ConserveMemory = false;
        }

        #endregion Public Constructors

        #region Public Properties

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }

        public int? MinPeptideLength { get; set; }

        public int? MaxPeptideLength { get; set; }
        public bool ConserveMemory { get; set; }

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

            MassDiffAcceptor searchMode = new DotMassDiffAcceptor("", gptmdModifications.Select(b => b.monoisotopicMass).Concat(GetObservedMasses(variableModifications.Concat(fixedModifications), gptmdModifications)).Concat(combos.Select(b => b.Item1 + b.Item2)).Concat(new List<double> { 0 }).GroupBy(b => Math.Round(b, 6)).Select(b => b.FirstOrDefault()).OrderBy(b => b), PrecursorMassTolerance);

            var searchModes = new List<MassDiffAcceptor> { searchMode };

            List<PsmParent>[] allPsms = new List<PsmParent>[1];
            allPsms[0] = new List<PsmParent>();

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

            FdrAnalysisResults analysisResults = null;
            var numRawFiles = currentRawFileList.Count;

            Status("Running G-PTM-D...", new List<string> { taskId });

            for (int spectraFileIndex = 0; spectraFileIndex < numRawFiles; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];

                NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Spectra Files", origDataFile });
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                else
                    myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);

                Status("Getting ms2 scans...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                bool findAllPrecursors = true;
                bool useProvidedPrecursorInfo = true;
                var intensityRatio = 4;
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, origDataFile).OrderBy(b => b.PrecursorMass).ToArray();
                var searchResults = (SearchResults)new ClassicSearchEngine(arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Spectra Files", origDataFile }, ConserveMemory).Run();

                allPsms[0].AddRange(searchResults.Psms[0]);

                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
            }
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            List<PsmParent>[] allPsmsTest = new List<PsmParent>[1];
            allPsmsTest[0] = allPsms[0].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, double>(b.FullFilePath, b.ScanNumber, b.PeptideMonoisotopicMass)).Select(b => b.First()).ToList();
            // Group and order psms

            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngineTest = new SequencesToActualProteinPeptidesEngine(allPsmsTest, modsDictionary, proteinList, searchModes, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, InitiatorMethionineBehavior, fixedModifications, variableModifications, MaxModificationIsoforms);
            var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngineTest.Run();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

            analysisResults = (FdrAnalysisResults)new FdrAnalysisEngine(allPsmsTest, compactPeptideToProteinPeptideMatchingTest,
                        searchModes,
                        false, false, false,
                        new List<string> { taskId }).Run();

            var gptmdResults = (GptmdResults)new GptmdEngine(allPsmsTest[0], gptmdModifications, combos, PrecursorMassTolerance).Run();

            if (currentXmlDbFilenameList.Any(b => !b.IsContaminant))
            {
                string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", currentXmlDbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

                var newModsActuallyWritten = ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });

                myTaskResults.newDatabases.Add(new DbForTask(outputXMLdbFullName, false));
                myTaskResults.AddNiceText("Modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                myTaskResults.AddNiceText("Mods types and counts:");
                myTaskResults.AddNiceText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            }
            if (currentXmlDbFilenameList.Any(b => b.IsContaminant))
            {
                string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", currentXmlDbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "GPTMD.xml");

                var newModsActuallyWritten = ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });

                myTaskResults.newDatabases.Add(new DbForTask(outputXMLdbFullNameContaminants, true));
                myTaskResults.AddNiceText("Contaminant modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                myTaskResults.AddNiceText("Mods types and counts:");
                myTaskResults.AddNiceText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
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

        private IEnumerable<Tuple<double, double>> LoadCombos(List<ModificationWithMass> modificationsThatCanBeCombined)
        {
            using (StreamReader r = new StreamReader(Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Data", @"combos.txt")))
            {
                while (r.Peek() >= 0)
                {
                    var line = r.ReadLine().Split(' ');
                    var mass1 = double.Parse(line[0]);
                    var mass2 = double.Parse(line[1]);
                    if (modificationsThatCanBeCombined.Any(b => Math.Abs(b.monoisotopicMass - mass1) < 1e-3) &&
                        modificationsThatCanBeCombined.Any(b => Math.Abs(b.monoisotopicMass - mass2) < 1e-3))
                        yield return new Tuple<double, double>(mass1, mass2);
                }
            }
        }

        #endregion Private Methods

    }
}