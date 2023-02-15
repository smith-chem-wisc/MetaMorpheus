using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using EngineLayer;
using GuiFunctions;
using SpectralAveraging;
using Nett;
using TaskLayer;
using System.IO;
using MassSpectrometry;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SpectralAveragingTaskWindow.xaml
    /// </summary>
    public partial class SpectralAveragingTaskWindow : Window
    {
        public SpectralAveragingTaskWindow() : this(null)
        {
        }

        internal SpectralAveragingTask TheTask { get; private set; }
        public SpectralAveragingTaskWindow(SpectralAveragingTask task)
        {
            InitializeComponent();
            TheTask = task ?? new SpectralAveragingTask(new SpectralAveragingParameters());
            DataContext = new SpectralAveragingParametersViewModel(TheTask.Parameters);
        }

        private void CancelButton_OnClick(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_OnClick(object sender, RoutedEventArgs e)
        {
            // TODO: Add max threads per file as a param that can be set on the GUI for task
            // setting dissociation type ensures no filters occur when loading files in MM Task
            CommonParameters commParamsToSave = new CommonParameters(
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "AveragingTask", 
                dissociationType: DissociationType.LowCID,
                maxThreadsToUsePerFile: 2);

            TheTask.CommonParameters = commParamsToSave;
            DialogResult = true;
        }

        private void SetDefaultbutton_OnClick(object sender, RoutedEventArgs e)
        {
            SaveButton_OnClick(sender, e);
            Toml.WriteFile(TheTask, Path.Combine(GlobalVariables.DataDir, "DefaultParameters", "SpectralAverageTaskDefault.toml"), MetaMorpheusTask.tomlConfig);
        }

        private void SpectralAveragingTaskWindow_OnKeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                SaveButton_OnClick(sender, e);
            }
            else if (e.Key == Key.Escape)
            {
                CancelButton_OnClick(sender, e);
            }
        }
    }
}
