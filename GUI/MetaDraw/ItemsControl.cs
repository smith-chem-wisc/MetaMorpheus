using System.ComponentModel;
using ViewModels;
using System.Collections.ObjectModel;

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
                SpectrumLabel = "Spectra info here"
            };

            Data = new ObservableCollection<ItemsControlSampleData>();
            Data.Add(sampledata);
        }

        public void AddNewRow(PsmAnnotationViewModel psmAnnotationViewModel, string annotation)
        {
            Data.Add(new ItemsControlSampleData() { PsmAnnotationViewModel = psmAnnotationViewModel, SpectrumLabel = annotation });
        }
    }

    public class ItemsControlSampleData
    {
        public PsmAnnotationViewModel PsmAnnotationViewModel { get; set; }

        public string SpectrumLabel { get; set; }

    }
}