using System;
using System.Collections.ObjectModel;
using System.Text;
using EngineLayer;
using MzLibUtil;
using Nett;
using Proteomics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;
using System.Collections.Generic;

namespace MetaMorpheusGUI
{
    internal static class NeoLoadTomls
    {
        private static string novelAddition = @"NeoTomlFiles\";
        private static string defaultAddition = @"EngineLayer\Neo\TomlFiles\";

        public static void LoadTomls(NeoSearchTask ye5, ObservableCollection<PreRunTask> staticTasksObservableCollection, string draggedFilePath)
        {
            string novelFolderPath = "";
            string[] splitPath = draggedFilePath.Split('\\').ToArray();
            for (int i = 0; i < splitPath.Length - 1; i++)
                novelFolderPath += splitPath[i] + '\\';
            novelFolderPath += novelAddition;

            string defaultFolderPath = "";
            splitPath = Directory.GetCurrentDirectory().Split('\\').ToArray();
            for (int i = 0; i < splitPath.Length - 3; i++)
                defaultFolderPath += splitPath[i] + '\\';
            defaultFolderPath += defaultAddition;

            #region write TOML

            var tomlFileName = Path.Combine(defaultFolderPath, ye5.GetType().Name + "config.toml");
            Toml.WriteFile(ye5, tomlFileName, MetaMorpheusTask.tomlConfig);

            #endregion write TOML

            if (ye5.NeoParameters.Calibrate)
            {
                string caliFileName = "CalibrationTaskconfig.toml";
                caliFileName = File.Exists(novelFolderPath + caliFileName) ? novelFolderPath + caliFileName : defaultFolderPath + caliFileName;
                UpdateTomls(tomlFileName, caliFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<CalibrationTask>(caliFileName, MetaMorpheusTask.tomlConfig);//FIXME, An item with the same key has already been added, dictionary in Toml.ReadFile
                AddTaskToCollection(yeo, staticTasksObservableCollection);//multiple protease issue
            }

            if (ye5.NeoParameters.GPTMD)
            {
                string gptmdFileName = "GptmdTaskconfig.toml";
                gptmdFileName = File.Exists(novelFolderPath + gptmdFileName) ? novelFolderPath + gptmdFileName : defaultFolderPath + gptmdFileName;
                UpdateTomls(tomlFileName, gptmdFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<GptmdTask>(gptmdFileName, MetaMorpheusTask.tomlConfig);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            if (ye5.NeoParameters.TargetSearch)
            {
                string targetFileName = "SearchTaskTargetconfig.toml";
                targetFileName = File.Exists(novelFolderPath + targetFileName) ? novelFolderPath + targetFileName : defaultFolderPath + targetFileName;
                UpdateTomls(tomlFileName, targetFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            if (ye5.NeoParameters.DecoySearch)
            {
                string decoyFileName = "SearchTaskDecoyconfig.toml";
                decoyFileName = File.Exists(novelFolderPath + decoyFileName) ? novelFolderPath + decoyFileName : defaultFolderPath + decoyFileName;
                UpdateTomls(tomlFileName, decoyFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<SearchTask>(decoyFileName, MetaMorpheusTask.tomlConfig);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            var yeo5_1 = ye5.Clone();
            yeo5_1.AggregateTargetDecoyFiles = true;
            yeo5_1.AggregateNormalSplicedFiles = false;
            yeo5_1.GenerateSplicedPeptides = false;

            AddTaskToCollection(yeo5_1, staticTasksObservableCollection);

            if (ye5.NeoParameters.SearchNTerminus)
            {
                string targetFileName = "SearchTaskNconfig.toml";
                targetFileName = File.Exists(novelFolderPath + targetFileName) ? novelFolderPath + targetFileName : defaultFolderPath + targetFileName;
                UpdateTomls(tomlFileName, targetFileName, ye5.CommonParameters, TerminusType.N, false);
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            if (ye5.NeoParameters.SearchCTerminus)
            {
                string targetFileName = "SearchTaskCconfig.toml";
                targetFileName = File.Exists(novelFolderPath + targetFileName) ? novelFolderPath + targetFileName : defaultFolderPath + targetFileName;
                UpdateTomls(tomlFileName, targetFileName, ye5.CommonParameters, TerminusType.C, false);
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            var yeo5_2 = ye5.Clone();
            yeo5_2.AggregateTargetDecoyFiles = false;
            yeo5_2.AggregateNormalSplicedFiles = false;
            yeo5_2.GenerateSplicedPeptides = true;

            AddTaskToCollection(yeo5_2, staticTasksObservableCollection);

            string cisFileName = "SearchTaskCisconfig.toml";
            cisFileName = File.Exists(novelFolderPath + cisFileName) ? novelFolderPath + cisFileName : defaultFolderPath + cisFileName;
            UpdateTomls(tomlFileName, cisFileName, ye5.CommonParameters, TerminusType.None, true);
            var yeocis = Toml.ReadFile<SearchTask>(cisFileName, MetaMorpheusTask.tomlConfig);
            AddTaskToCollection(yeocis, staticTasksObservableCollection);

            string transFileName = "SearchTaskTransconfig.toml";
            transFileName = File.Exists(novelFolderPath + transFileName) ? novelFolderPath + transFileName : defaultFolderPath + transFileName;
            UpdateTomls(tomlFileName, transFileName, ye5.CommonParameters, TerminusType.None, true);
            var yeotrans = Toml.ReadFile<SearchTask>(transFileName, MetaMorpheusTask.tomlConfig);
            AddTaskToCollection(yeotrans, staticTasksObservableCollection);

            var yeo5_3 = ye5.Clone();
            yeo5_3.AggregateTargetDecoyFiles = false;
            yeo5_3.AggregateNormalSplicedFiles = true;
            yeo5_3.GenerateSplicedPeptides = false;

            #region DeleteTomlFile

            File.Delete(tomlFileName);

            #endregion DeleteTomlFile

            AddTaskToCollection(yeo5_3, staticTasksObservableCollection);
        }

        private static void AddTaskToCollection(MetaMorpheusTask yeo, ObservableCollection<PreRunTask> staticTasksObservableCollection)
        {
            PreRunTask teo = new PreRunTask(yeo);
            staticTasksObservableCollection.Add(teo);
            staticTasksObservableCollection.Last().DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(teo) + 1) + "-" + yeo.CommonParameters.TaskDescriptor;
        }

        private static void UpdateTomls(string tomlFileName, string fileName, ICommonParameters ye5, TerminusType terminusType, bool spliceSearch)
        {
            string[] oldTomlLines = File.ReadAllLines(@fileName);
            List<string> newTomlLines = new List<string>();
            foreach (string line in oldTomlLines)
            {
                if (line.Contains("LocalizeAll") && terminusType.Equals(TerminusType.None))
                    newTomlLines.Add(GetCorrectValue("LocalizeAll", tomlFileName, line));
                else if (line.Contains("ListOfModsFixed"))
                    newTomlLines.Add(GetCorrectValue("ListOfModsFixed", tomlFileName, line));
                else if (line.Contains("ListOfModsVariable") && terminusType.Equals(TerminusType.None))
                    newTomlLines.Add(GetCorrectValue("ListOfModsVariable", tomlFileName, line));
                else if (line.Contains("BIons"))
                {
                    if (terminusType.Equals(TerminusType.N) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("BIons", tomlFileName, line));
                    else
                        newTomlLines.Add("BIons = false");
                }
                else if (line.Contains("YIons"))
                {
                    if (terminusType.Equals(TerminusType.C) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("YIons", tomlFileName, line));
                    else
                        newTomlLines.Add("YIons = false");
                }
                else if (line.Contains("ZdotIons"))
                {
                    if (terminusType.Equals(TerminusType.C) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("ZdotIons", tomlFileName, line));
                    else
                        newTomlLines.Add("ZdotIons = false");
                }
                else if (line.Contains("CIons"))
                {
                    if (terminusType.Equals(TerminusType.N) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("CIons", tomlFileName, line));
                    else
                        newTomlLines.Add("CIons = false");
                }
                else if (line.Contains("ProductMassTolerance"))
                    newTomlLines.Add(GetCorrectValue("ProductMassTolerance", tomlFileName, line));
                else if (line.Contains("PrecursorMassTolerance"))
                    newTomlLines.Add(GetCorrectValue("PrecursorMassTolerance", tomlFileName, line));
                else if (line.Contains("MaxMissedCleavages"))
                    newTomlLines.Add(GetCorrectValue("MaxMissedCleavages", tomlFileName, line));
                else if (line.Contains("InitiatorMethionineBehavior"))
                    newTomlLines.Add(GetCorrectValue("InitiatorMethionineBehavior", tomlFileName, line));
                else if (line.Contains("MinPeptideLength") && !!terminusType.Equals(TerminusType.None))
                    newTomlLines.Add(GetCorrectValue("MinPeptideLength", tomlFileName, line));
                else if (line.Contains("MaxPeptideLength"))
                    newTomlLines.Add(GetCorrectValue("MaxPeptideLength", tomlFileName, line));
                else if (line.Contains("MaxModificationIsoforms"))
                    newTomlLines.Add(GetCorrectValue("MaxModificationIsoforms", tomlFileName, line));
                else if (line.Contains("MaxModsForPeptide"))
                    newTomlLines.Add(GetCorrectValue("MaxModsForPeptide", tomlFileName, line));
                else if (line.Contains("SemiProteaseDigestion"))
                    newTomlLines.Add(GetCorrectValue("SemiProteaseDigestion", tomlFileName, line));
                else if (line.Contains("TerminusTypeSemiProtease"))
                    newTomlLines.Add(GetCorrectValue("TerminusTypeSemiProtease", tomlFileName, line));
                else if (line.Contains("Protease") && terminusType.Equals(TerminusType.None) && !spliceSearch) //this must be last, else other names including protease will be overwritten and crash.
                    newTomlLines.Add(GetCorrectValue("Protease", tomlFileName, line));
                else
                    newTomlLines.Add(line);
            }
            using (StreamWriter file = new StreamWriter(fileName))
                foreach (string line in newTomlLines)
                    file.WriteLine(line);
        }

        private static string GetCorrectValue(string parameter, string tomlFileName, string oldLine)
        {
            string[] newTomlLines = File.ReadAllLines(@tomlFileName);
            foreach (string line in newTomlLines)
                if (line.Contains(parameter))
                    return line;
            return oldLine;
        }
    }
}
