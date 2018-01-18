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

namespace MetaMorpheusGUI
{
    internal static class NeoLoadTomls
    {
        public static string novelAddition = @"NeoTomlFiles\";
        public static string defaultAddition = @"EngineLayer\Neo\TomlFiles\";

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

            if (ye5.NeoParameters.Calibrate)
            {
                string caliFileName = "CalibrationTaskconfig.toml";
                caliFileName = File.Exists(novelFolderPath + caliFileName) ? novelFolderPath + caliFileName : defaultFolderPath + caliFileName;
                var yeo = Toml.ReadFile<CalibrationTask>(caliFileName, MetaMorpheusTask.tomlConfig);
                yeo.CommonParameters = ye5.CommonParameters.Clone(yeo.CommonParameters.TaskDescriptor);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            if (ye5.NeoParameters.GPTMD)
            {
                string gptmdFileName = "GptmdTaskconfig.toml";
                gptmdFileName = File.Exists(novelFolderPath + gptmdFileName) ? novelFolderPath + gptmdFileName : defaultFolderPath + gptmdFileName;
                var yeo = Toml.ReadFile<GptmdTask>(gptmdFileName, MetaMorpheusTask.tomlConfig);
                yeo.CommonParameters = ye5.CommonParameters.Clone(yeo.CommonParameters.TaskDescriptor);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            if (ye5.NeoParameters.TargetSearch)
            {
                string targetFileName = "SearchTaskTargetconfig.toml";
                targetFileName = File.Exists(novelFolderPath + targetFileName) ? novelFolderPath + targetFileName : defaultFolderPath + targetFileName;
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                yeo.CommonParameters = ye5.CommonParameters.Clone(yeo.CommonParameters.TaskDescriptor);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            if (ye5.NeoParameters.DecoySearch)
            {
                string targetFileName = "SearchTaskDecoyconfig.toml";
                targetFileName = File.Exists(novelFolderPath + targetFileName) ? novelFolderPath + targetFileName : defaultFolderPath + targetFileName;
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                yeo.CommonParameters = ye5.CommonParameters.Clone(yeo.CommonParameters.TaskDescriptor);
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
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                yeo.CommonParameters = ye5.CommonParameters.Clone(yeo.CommonParameters.TaskDescriptor);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            if (ye5.NeoParameters.SearchCTerminus)
            {
                string targetFileName = "SearchTaskCconfig.toml";
                targetFileName = File.Exists(novelFolderPath + targetFileName) ? novelFolderPath + targetFileName : defaultFolderPath + targetFileName;
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                yeo.CommonParameters = ye5.CommonParameters.Clone(yeo.CommonParameters.TaskDescriptor);
                AddTaskToCollection(yeo, staticTasksObservableCollection);
            }

            var yeo5_2 = ye5.Clone();
            yeo5_2.AggregateTargetDecoyFiles = false;
            yeo5_2.AggregateNormalSplicedFiles = false;
            yeo5_2.GenerateSplicedPeptides = true;

            AddTaskToCollection(yeo5_2, staticTasksObservableCollection);

            string cisFileName = "SearchTaskCisconfig.toml";
            cisFileName = File.Exists(novelFolderPath + cisFileName) ? novelFolderPath + cisFileName : defaultFolderPath + cisFileName;
            var yeocis = Toml.ReadFile<SearchTask>(cisFileName, MetaMorpheusTask.tomlConfig);
            yeocis.CommonParameters = ye5.CommonParameters.Clone(yeocis.CommonParameters.TaskDescriptor);
            AddTaskToCollection(yeocis, staticTasksObservableCollection);

            string transFileName = "SearchTaskTransconfig.toml";
            transFileName = File.Exists(novelFolderPath + transFileName) ? novelFolderPath + transFileName : defaultFolderPath + transFileName;
            var yeotrans = Toml.ReadFile<SearchTask>(transFileName, MetaMorpheusTask.tomlConfig);
            yeotrans.CommonParameters = ye5.CommonParameters.Clone(yeotrans.CommonParameters.TaskDescriptor);
            AddTaskToCollection(yeotrans, staticTasksObservableCollection);

            var yeo5_3 = ye5.Clone();
            yeo5_3.AggregateTargetDecoyFiles = false;
            yeo5_3.AggregateNormalSplicedFiles = true;
            yeo5_3.GenerateSplicedPeptides = false;

            AddTaskToCollection(yeo5_3, staticTasksObservableCollection);
        }

        private static void AddTaskToCollection(MetaMorpheusTask yeo, ObservableCollection<PreRunTask> staticTasksObservableCollection)
        {
            PreRunTask teo = new PreRunTask(yeo);
            staticTasksObservableCollection.Add(teo);
            staticTasksObservableCollection.Last().DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(teo) + 1) + "-" + yeo.CommonParameters.TaskDescriptor;
        }
    }
}
