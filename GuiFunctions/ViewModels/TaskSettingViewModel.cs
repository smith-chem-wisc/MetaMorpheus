using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Input;
using Easy.Common.Extensions;
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

        public MetaMorpheusTask TheTask { get; set; }


        #endregion

        #region Constructor

        public TaskSettingViewModel(MetaMorpheusTask task)
        {
            TheTask = task;
            var temp = typeof(SearchTask);
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
        }

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
        }


            #endregion

            #region Command Methods

        public void SaveSettings()
        {
            if (SelectedSettings == null)
            {
                return;
            }
            AllSettingsDict.Add(SelectedSettings, TheTask);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            TomlFileFolderSerializer.Save(SelectedSettings, TheTask);
        }

        public void SaveSettingsFromWindow()
        {

            if (TypedSettingsName is null || TypedSettingsName.IsNullOrEmptyOrWhiteSpace())
            {
                MessageBox.Show("The name cannot be empty", "Save Settings Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            string settingsName = TypedSettingsName.ToString();

            if (AllSettingsDict.ContainsKey(settingsName))
            {
                MessageBox.Show("The name already exists", "Repeated Name Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            AllSettingsDict.Add(settingsName, TheTask);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            TomlFileFolderSerializer.Save(settingsName, TheTask);
            
        }

        public void DeleteSettings()
        {
            if(SelectedSettings == null)
            {
                return;
            }
            AllSettingsDict.Remove(SelectedSettings);
            OnPropertyChanged(nameof(AllSettings));
            OnPropertyChanged(nameof(AllSettingsDict));
            TomlFileFolderSerializer.Delete(TheTask.GetType(), SelectedSettings);
        }

        #endregion

        #region Helpers

   

        #endregion

    }
}
