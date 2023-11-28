using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Input;
using Easy.Common.Extensions;
using EngineLayer;
using Nett;
using Org.BouncyCastle.Asn1.X509.Qualified;
using OxyPlot;
using pepXML.Generated;
using TaskLayer;
using ThermoFisher.CommonCore.Data;
using TopDownProteomics.Chemistry;

namespace GuiFunctions
{
    public class TaskSettingViewModel : BaseViewModel
    {
        #region Private Properties

        private string selectedSettings;
        private string typedSettingsName;
        private Dictionary<string, MetaMorpheusTask> allSettingsDict;
        private ObservableCollection<string> allSettings => new ObservableCollection<string>(allSettingsDict.Keys);

        #endregion

        #region Public Properties

        public ObservableCollection<string> AllSettings => allSettings;

        public string TypedSettingsName
        {
            get => typedSettingsName;
            set
            {
                typedSettingsName = value;
                OnPropertyChanged(nameof(TypedSettingsName));
            }
        }

        public string SelectedSettings
        {
            get => selectedSettings;
            set
            {
                selectedSettings = value;
                if (UpdateFieldsInGuiWithNewTask is not null)
                    UpdateFieldsInGuiWithNewTask.Invoke(AllSettingsDict[SelectedSettings]);
                OnPropertyChanged(nameof(SelectedSettings));
            }
        }

        public Dictionary<string, MetaMorpheusTask> AllSettingsDict
        {
            get => allSettingsDict;
            set
            {
                allSettingsDict = value;
                OnPropertyChanged(nameof(AllSettingsDict));
            }
        }

      
        public ICommand SaveSettingsCommand { get; set; }
        public ICommand DeleteSettingsCommand { get; set; }
        public ICommand SaveAsDefaultSettingsCommand { get; set; }

        /// <summary>
        /// The task that is currently working on
        /// </summary>
        public MetaMorpheusTask TheTask => AllSettingsDict[SelectedSettings];

        /// <summary>
        /// Task from gui
        /// </summary>
        public Func<MetaMorpheusTask> GetTaskFromGui { get; set; }

        /// <summary>
        /// update the gui from the task
        /// </summary>
        public Action<MetaMorpheusTask> UpdateFieldsInGuiWithNewTask { get; set; }

        #endregion

        #region Constructor

        public TaskSettingViewModel(MetaMorpheusTask task, Action<MetaMorpheusTask> updateFieldsInGuiEventHandler, Func<MetaMorpheusTask> getTaskFromGui)
        {
            // if directory does not exist, create new directory and default settings
            if (!Directory.Exists(TomlFileFolderSerializer.PathToCheck))
            {
                Directory.CreateDirectory(Path.Combine(TomlFileFolderSerializer.PathToCheck));
            }

            // check if default is saved in previous location, if so, move to new location and rename
            if (Directory.Exists(TomlFileFolderSerializer.PathToCheck))
            {
                string oldDefaultPathToCheck = task switch
                {
                    SearchTask search => "SearchTaskDefault.toml",
                    CalibrationTask calibration => "CalibrationTaskDefault.toml",
                    GptmdTask gptmd => "GptmdTaskDefault.toml",
                    XLSearchTask xLSearch => "XLSearchTaskDefault.toml",
                    GlycoSearchTask glycoSearch => "GlycoSearchTaskDefault.toml",
                    SpectralAveragingTask spectralAverage => "SpectralAverageTaskDefault.toml",
                    _ => "",
                };

                String finalPath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", oldDefaultPathToCheck);
                if (File.Exists(finalPath))
                {
                    // load in
                    MetaMorpheusTask previouslySavedDefaultTask = task switch
                    {
                        SearchTask search => Toml.ReadFile<SearchTask>(finalPath, MetaMorpheusTask.tomlConfig),
                        CalibrationTask calibration => Toml.ReadFile<CalibrationTask>(finalPath, MetaMorpheusTask.tomlConfig),
                        GptmdTask gptmd => Toml.ReadFile<GptmdTask>(finalPath, MetaMorpheusTask.tomlConfig),
                        XLSearchTask xLSearch => Toml.ReadFile<XLSearchTask>(finalPath, MetaMorpheusTask.tomlConfig),
                        GlycoSearchTask glycoSearch => Toml.ReadFile<GlycoSearchTask>(finalPath, MetaMorpheusTask.tomlConfig),
                        SpectralAveragingTask spectralAverage => Toml.ReadFile<SpectralAveragingTask>(finalPath, MetaMorpheusTask.tomlConfig),
                        _ => Toml.ReadFile<MetaMorpheusTask>(finalPath, MetaMorpheusTask.tomlConfig),
                    };


                    String s = oldDefaultPathToCheck.Substring(0, oldDefaultPathToCheck.Length - 5) + "(Default)";
                    TomlFileFolderSerializer.Save(oldDefaultPathToCheck.Substring(0, oldDefaultPathToCheck.Length - 5) + "(Default)", previouslySavedDefaultTask);
                    TomlFileFolderSerializer.Delete(previouslySavedDefaultTask.GetType(), oldDefaultPathToCheck);
                    File.Delete(finalPath);
                }
            }

            // initial construction
            var type = task.GetType();
            UpdateFieldsInGuiWithNewTask = updateFieldsInGuiEventHandler;
            GetTaskFromGui = getTaskFromGui;
            SaveSettingsCommand = new RelayCommand(SaveSettings);
            DeleteSettingsCommand = new RelayCommand(DeleteSettings);
            SaveAsDefaultSettingsCommand = new RelayCommand(SaveAsDefaultSettings);

            // load in all settings from directory of specific object type
            MethodInfo method = typeof(TomlFileFolderSerializer).GetMethod("LoadAllOfTypeT");
            MethodInfo genericMethod = method.MakeGenericMethod(type);
            System.Collections.IDictionary result = ((System.Collections.IDictionary)genericMethod.Invoke(null, new object[] { task }));

            // parse directory to get object type and MMTask from toml files
            var keys = result.Keys.Cast<string>().ToList();
            var values = result.Values.Cast<MetaMorpheusTask>().ToList();
            Dictionary<string, MetaMorpheusTask> tempDict = new();
            for (int i = 0; i < keys.Count; i++)
            {
                tempDict.Add(keys[i], values[i]);
            }
            AllSettingsDict = tempDict;
            
            if (!AllSettingsDict.Any())
            {
                string name = "DefaultSetting(Default)";
                TomlFileFolderSerializer.Save("DefaultSetting(Default)", task);
                AllSettingsDict.Add(name, task);
                SelectedSettings = AllSettingsDict.First(p => p.Key.Contains("(Default)")).Key;
            }

            if (AllSettingsDict.Any(p => p.Key.Contains("(Default)")))
            {
                SelectedSettings = AllSettingsDict.First(p => p.Key.Contains("(Default)")).Key;
            }
        }

        
        #endregion



        #region Command Methods

        /// <summary>
        /// Modify a saved task
        /// </summary>
        public void SaveSettings()
        {
            if (SelectedSettings == null)
            {
                return;
            }

            if (SelectedSettings.Contains("(Default)"))
            {
                if (!GlobalVariables.MetaMorpheusVersion.Contains("DEBUG"))
                {
                    MessageBox.Show("Default Setting cannot be modified. Please use \"Save as\" button",
                        "Modifying Default Settings Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                }

                return;
            }

            var taskFromGUI = GetTaskFromGui.Invoke();
            TomlFileFolderSerializer.Delete(TheTask.GetType(), SelectedSettings);
            AllSettingsDict.Remove(SelectedSettings);
            TomlFileFolderSerializer.Save(SelectedSettings, taskFromGUI);
            AllSettingsDict.Add(SelectedSettings, taskFromGUI);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            SelectedSettings = this.SelectedSettings;
        }

        /// <summary>
        /// Save a new task with a window popped up
        /// </summary>
        public void SaveSettingsFromWindow()
        {

            if (TypedSettingsName is null || TypedSettingsName.IsNullOrEmptyOrWhiteSpace())
            {
                if (!GlobalVariables.MetaMorpheusVersion.Contains("DEBUG"))
                {
                    MessageBox.Show("The name cannot be empty", "Empty Settings Error", MessageBoxButton.OK,
                        MessageBoxImage.Warning);
                    TypedSettingsName = "";
                }

                return;
            }

            string settingsName = TypedSettingsName.ToString();

            if (AllSettingsDict.ContainsKey(settingsName))
            {
                if (!GlobalVariables.MetaMorpheusVersion.Contains("DEBUG"))
                {
                    MessageBox.Show("The name already exists", "Repeated Name Error", MessageBoxButton.OK,
                        MessageBoxImage.Warning);
                    TypedSettingsName = "";
                }

                return;
            }

            // get current gui representation
            var taskFromGUI = GetTaskFromGui.Invoke();
            AllSettingsDict.Add(settingsName, taskFromGUI);

            // revert modified settings before saving
            string settingsThatWereModifiedButNotMeantToBeSavedPath = TomlFileFolderSerializer.GetFilePath(taskFromGUI.GetType(), SelectedSettings);
            var type = taskFromGUI.GetType();
            MethodInfo method = typeof(TomlFileFolderSerializer).GetMethod("Deserialize");
            MethodInfo genericMethod = method.MakeGenericMethod(type);
            var reloadedModifiedSettings = genericMethod.Invoke(null, new object[] { settingsThatWereModifiedButNotMeantToBeSavedPath });
            AllSettingsDict[SelectedSettings] = reloadedModifiedSettings as MetaMorpheusTask;

            // save modified settings with new name
            SelectedSettings = settingsName;
            TomlFileFolderSerializer.Save(settingsName, taskFromGUI);

            // update gui representation
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            TypedSettingsName = "";
            //SelectDefaultAuto();
        }

        /// <summary>
        /// Delete the selected task
        /// </summary>
        public void DeleteSettings()
        {
            if(SelectedSettings == null)
            {
                return;
            }

            if (SelectedSettings.Contains("(Default)"))
            {
                if(!GlobalVariables.MetaMorpheusVersion.Contains("DEBUG")){
                    MessageBox.Show("Default Setting cannot be deleted", "Deleting Default Settings Error",
                        MessageBoxButton.OK, MessageBoxImage.Warning);
                }

                return;
            }

            TomlFileFolderSerializer.Delete(TheTask.GetType(), SelectedSettings);
            AllSettingsDict.Remove(SelectedSettings);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            //SelectDefaultAuto();
        }

        /// <summary>
        /// Set a saved set as a default setting
        /// </summary>
        public void SaveAsDefaultSettings()
        {
            if (SelectedSettings == null)
            {
                return;
            }

            if (SelectedSettings.Contains("(Default)"))
            {
                if (!GlobalVariables.MetaMorpheusVersion.Contains("DEBUG"))
                {
                    MessageBox.Show("It's already a default setting", "Repeatted Setting Default Settings Error",
                        MessageBoxButton.OK, MessageBoxImage.Warning);
                }

                return;
            }

            var temp = TheTask;
            string currentName = SelectedSettings;

            // if there is currently a default, remove "(Default)" and save with new name
            if (AllSettingsDict.Any(p => p.Key.Contains("(Default)")))
            {
                var oldDefault = AllSettingsDict.First(p => p.Key.Contains("(Default)"));
                string newOldDefaultName = oldDefault.Key.Remove(oldDefault.Key.Length - 9);
                TomlFileFolderSerializer.Delete(oldDefault.Value.GetType(), oldDefault.Key);
                AllSettingsDict.Remove(oldDefault.Key);
                TomlFileFolderSerializer.Save(newOldDefaultName, oldDefault.Value);
                AllSettingsDict.Add(newOldDefaultName, oldDefault.Value);
            }

            // save currently selected settings as new default
            TomlFileFolderSerializer.Delete(TheTask.GetType(), SelectedSettings); 
            AllSettingsDict.Remove(SelectedSettings);
            TomlFileFolderSerializer.Save(currentName + "(Default)", temp);
            AllSettingsDict.Add(currentName + "(Default)", temp);

            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            //SelectDefaultAuto();
        }

        

    #endregion

        #region Helpers
        
    }

    #endregion

}

