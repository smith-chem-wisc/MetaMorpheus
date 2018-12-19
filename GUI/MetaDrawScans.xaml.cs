using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
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
using System.Windows.Shapes;
using ViewModels;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawScans.xaml
    /// </summary>
    public partial class MetaDrawScans : Window
    {
        public MetaDrawScans()
        {
            InitializeComponent();
            DataContext = new ItemsControlSampleViewModel();
        }

        public class ItemsControlSampleViewModel
        {
            public ObservableCollection<ItemsControlSampleData> Data { get; set; }

            public Command AddNewRowCommand { get; set; }

            public Command<ItemsControlSampleData> Command1 { get; set; }

            public ItemsControlSampleViewModel()
            {
                var sampledata = Enumerable.Range(0, 1)
                                           .Select(x => new ItemsControlSampleData()
                                           {
                                               PsmAnnotationViewModel = new PsmAnnotationViewModel(),
                                               Label1Text = "Label1 " + x.ToString(),
                                           });

                Data = new ObservableCollection<ItemsControlSampleData>(sampledata);
                AddNewRowCommand = new Command(AddNewRow);
                Command1 = new Command<ItemsControlSampleData>(ExecuteCommand1);

            }

            private void AddNewRow()
            {
                Data.Add(new ItemsControlSampleData() { PsmAnnotationViewModel = new PsmAnnotationViewModel(), Label1Text = "Label 1 - New Row"});
            }

            private void ExecuteCommand1(ItemsControlSampleData data)
            {
                MessageBox.Show("Command1 - " + data.Label1Text);
            }
        }

        public class ItemsControlSampleData
        {
            public PsmAnnotationViewModel PsmAnnotationViewModel { get; set; }

            public string Label1Text { get; set; }

        }

        public class Command : ICommand
        {
            public Action Action { get; set; }

            public string DisplayName { get; set; }

            public void Execute(object parameter)
            {
                if (Action != null)
                    Action();
            }

            public bool CanExecute(object parameter)
            {
                return IsEnabled;
            }

            private bool _isEnabled = true;
            public bool IsEnabled
            {
                get { return _isEnabled; }
                set
                {
                    _isEnabled = value;
                    if (CanExecuteChanged != null)
                        CanExecuteChanged(this, EventArgs.Empty);
                }
            }

            public event EventHandler CanExecuteChanged;

            public Command(Action action)
            {
                Action = action;
            }
        }

        public class Command<T> : ICommand
        {
            public Action<T> Action { get; set; }

            public void Execute(object parameter)
            {
                if (Action != null && parameter is T)
                    Action((T)parameter);
            }

            public bool CanExecute(object parameter)
            {
                return IsEnabled;
            }

            private bool _isEnabled = true;
            public bool IsEnabled
            {
                get { return _isEnabled; }
                set
                {
                    _isEnabled = value;
                    if (CanExecuteChanged != null)
                        CanExecuteChanged(this, EventArgs.Empty);
                }
            }

            public event EventHandler CanExecuteChanged;

            public Command(Action<T> action)
            {
                Action = action;
            }
        }
    }
}
