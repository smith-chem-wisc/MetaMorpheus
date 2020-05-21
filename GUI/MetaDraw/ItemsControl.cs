using System.ComponentModel;
using ViewModels;
using System.Collections.ObjectModel;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    public class ItemsControlSampleViewModel : INotifyPropertyChanged
    {
        private int _myColumnCount = 1;
        public event PropertyChangedEventHandler PropertyChanged = (s, e) => { };
        public int MyColumnCount
        {
            get { return _myColumnCount; }
            set
            {
                _myColumnCount = value;
                this.PropertyChanged(this, new PropertyChangedEventArgs("MyColumnCount"));
            }
        }

        public ObservableCollection<ItemsControlSampleData> Data { get; set; }

        public ItemsControlSampleViewModel()
        {
            var sampledata = new ItemsControlSampleData()
            {
                PsmAnnotationViewModel = new PsmAnnotationViewModel(),
                SpectrumLabel = "Spectra info here",
                TheCanvas = new Canvas()
            };

            Data = new ObservableCollection<ItemsControlSampleData>();
            Data.Add(sampledata);
        }

        public void AddNewRow(PsmAnnotationViewModel psmAnnotationViewModel, string annotation, Canvas canvas)
        {
            Data.Add(new ItemsControlSampleData { PsmAnnotationViewModel = psmAnnotationViewModel, SpectrumLabel = annotation, TheCanvas = canvas });
        }
    }

    public class ItemsControlSampleData
    {
        public PsmAnnotationViewModel PsmAnnotationViewModel { get; set; }

        public string SpectrumLabel { get; set; }

        public Canvas TheCanvas { get; set; }

    }
}