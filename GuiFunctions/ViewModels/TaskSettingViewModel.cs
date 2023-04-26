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
using Org.BouncyCastle.Asn1.X509.Qualified;
using OxyPlot;
using TaskLayer;

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

        public ObservableCollection<string> AllSettings
        {
            get => allSettings;
        }

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
        /// 
        /// </summary>
        public MetaMorpheusTask TheTask { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public Func<MetaMorpheusTask> GetTaskFromGui { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public Action<MetaMorpheusTask> UpdateFieldsInGuiWithNewTask { get; set; }

        #endregion

        #region Constructor

        public TaskSettingViewModel(MetaMorpheusTask task, Action<MetaMorpheusTask> updateFieldsInGuiEventHandler, Func<MetaMorpheusTask> getTaskFromGui)
        {
            
            TheTask = task;
            UpdateFieldsInGuiWithNewTask = updateFieldsInGuiEventHandler;
            GetTaskFromGui = getTaskFromGui;

            var type = task.GetType();
            MethodInfo method = typeof(TomlFileFolderSerializer).GetMethod("LoadAllOfTypeT");
            MethodInfo genericMethod = method.MakeGenericMethod(type);
            System.Collections.IDictionary result = ((System.Collections.IDictionary)genericMethod.Invoke(null, new object[] {TheTask}));

            var keys = result.Keys.Cast<string>().ToList();
            var values = result.Values.Cast<MetaMorpheusTask>().ToList();
            Dictionary<string, MetaMorpheusTask> tempDict = new();
            for(int i = 0; i < keys.Count; i++)
            {
                tempDict.Add(keys[i], values[i]);
            }
            AllSettingsDict = tempDict;

            
            SaveSettingsCommand = new RelayCommand(SaveSettings);
            DeleteSettingsCommand = new RelayCommand(DeleteSettings);
            SaveAsDefaultSettingsCommand = new RelayCommand(SaveAsDefaultSettings);
            SelectDefaultAuto();
        }

        /// <summary>
        /// ONLY USED FOR TESTING
        /// </summary>
        /// <param name="task"></param>
        internal TaskSettingViewModel(MetaMorpheusTask task)
        {
            TheTask = task;

            var type = task.GetType();
            MethodInfo method = typeof(TomlFileFolderSerializer).GetMethod("LoadAllOfTypeT");
            MethodInfo genericMethod = method.MakeGenericMethod(type);
            System.Collections.IDictionary result = ((System.Collections.IDictionary)genericMethod.Invoke(null, new object[] { TheTask }));

            var keys = result.Keys.Cast<string>().ToList();
            var values = result.Values.Cast<MetaMorpheusTask>().ToList();
            Dictionary<string, MetaMorpheusTask> tempDict = new();
            for (int i = 0; i < keys.Count; i++)
            {
                tempDict.Add(keys[i], values[i]);
            }
            AllSettingsDict = tempDict;

            SaveSettingsCommand = new RelayCommand(SaveSettings);
            DeleteSettingsCommand = new RelayCommand(DeleteSettings);
            SaveAsDefaultSettingsCommand = new RelayCommand(SaveAsDefaultSettings);
            SelectDefaultAuto();
        }

        /// <summary>
        ///  Used only for design time viewing of data 
        /// </summary>
        public TaskSettingViewModel()
        {
            MetaMorpheusTask task = new SearchTask();
            TheTask = task;
            var temp = typeof(SearchTask);
            var type = task.GetType();
            MethodInfo method = typeof(TomlFileFolderSerializer).GetMethod("LoadAllOfTypeT");
            MethodInfo genericMethod = method.MakeGenericMethod(type);
            System.Collections.IDictionary result = ((System.Collections.IDictionary)genericMethod.Invoke(null, new object[] { TheTask }));

            var keys = result.Keys.Cast<string>().ToList();
            var values = result.Values.Cast<MetaMorpheusTask>().ToList();
            Dictionary<string, MetaMorpheusTask> tempDict = new();
            for (int i = 0; i < keys.Count; i++)
            {
                tempDict.Add(keys[i], values[i]);
            }
            AllSettingsDict = tempDict;

            SaveSettingsCommand = new RelayCommand(SaveSettings);
            DeleteSettingsCommand = new RelayCommand(DeleteSettings);
            SaveAsDefaultSettingsCommand = new RelayCommand(SaveAsDefaultSettings);
            SelectDefaultAuto();
        }

        
        #endregion



        #region Command Methods

        public void SaveSettings()
        {
            if (SelectedSettings == null)
            {
                return;
            }

            if (SelectedSettings.Contains("(Default)"))
            {
                MessageBox.Show("Default Setting cannot be modified. Please use \"Save as Default\" button", "Modifying Default Settings Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }
            
            TomlFileFolderSerializer.Delete(TheTask.GetType(), SelectedSettings);
            AllSettingsDict.Remove(SelectedSettings);

            TheTask = GetTaskFromGui.Invoke();

            TomlFileFolderSerializer.Save(SelectedSettings, TheTask);
            AllSettingsDict.Add(SelectedSettings, TheTask);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            SelectDefaultAuto();
        }

    public void SaveSettingsFromWindow()
    {

            if (TypedSettingsName is null || TypedSettingsName.IsNullOrEmptyOrWhiteSpace())
            {
                MessageBox.Show("The name cannot be empty", "Empty Settings Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                TypedSettingsName = "";
                return;
            }

            string settingsName = TypedSettingsName.ToString();

            if (AllSettingsDict.ContainsKey(settingsName))
            {
                MessageBox.Show("The name already exists", "Repeated Name Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                TypedSettingsName = "";
                return;
            }

            // get current gui representation
            TheTask = GetTaskFromGui.Invoke();
            AllSettingsDict.Add(settingsName, TheTask);

            // revert modified settings before saving
            string settingsThatWereModifiedButNotMeantToBeSavedPath = TomlFileFolderSerializer.GetFilePath(TheTask.GetType(), SelectedSettings);
            var type = TheTask.GetType();
            MethodInfo method = typeof(TomlFileFolderSerializer).GetMethod("Deserialize");
            MethodInfo genericMethod = method.MakeGenericMethod(type);
            var reloadedModifiedSettings = genericMethod.Invoke(null, new object[] { settingsThatWereModifiedButNotMeantToBeSavedPath });
            AllSettingsDict[SelectedSettings] = reloadedModifiedSettings as MetaMorpheusTask;

            // save modified settings with new name
            SelectedSettings = settingsName;
            TomlFileFolderSerializer.Save(settingsName, TheTask);

            // update gui representation
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            TypedSettingsName = "";
            SelectDefaultAuto();
        }

        public void DeleteSettings()
        {
            if(SelectedSettings == null)
            {
                return;
            }

            if (SelectedSettings.Contains("(Default)"))
            {
                MessageBox.Show("Default Setting cannot be deleted", "Deleting Default Settings Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            TomlFileFolderSerializer.Delete(TheTask.GetType(), SelectedSettings);
            AllSettingsDict.Remove(SelectedSettings);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            SelectDefaultAuto();
        }

        public void SaveAsDefaultSettings()
        {
            if (SelectedSettings == null)
            {
                return;
            }

            if (SelectedSettings.Contains("(Default)"))
            {
                MessageBox.Show("It's already a default setting", "Repeatted Setting Default Settings Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            // if there is currently a default
            if (AllSettingsDict.Any(p => p.Key.Contains("(Default)")))
            {
                var oldDefault = AllSettingsDict.First(p => p.Key.Contains("(Default)"));
                string oldDefaultName = oldDefault.Key;
                string changedOldDefaultName = oldDefaultName.Remove(oldDefaultName.Length - 9);
                var oldDefaultValue = oldDefault.Value;
                TomlFileFolderSerializer.Delete(oldDefaultValue.GetType(), oldDefaultName);
                AllSettingsDict.Remove(oldDefaultName);
                TomlFileFolderSerializer.Save(changedOldDefaultName, oldDefaultValue);
                AllSettingsDict.Add(changedOldDefaultName, oldDefaultValue);               
            }

            var temp = TheTask;
            TomlFileFolderSerializer.Delete(TheTask.GetType(), SelectedSettings); 
            AllSettingsDict.Remove(SelectedSettings);
            TomlFileFolderSerializer.Save(SelectedSettings + "(Default)", temp);
            AllSettingsDict.Add(SelectedSettings + "(Default)", temp);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            SelectDefaultAuto();
        }

        

    #endregion

        #region Helpers
        private void SelectDefaultAuto()
        {
            var defaultSetting = AllSettingsDict.First(p => p.Key.Contains("(Default)"));
            SelectedSettings = defaultSetting.Key;
        }
   

        #endregion

    }
}
