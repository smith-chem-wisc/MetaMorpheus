using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.ObjectModel;
using System.ComponentModel;

namespace MetaMorpheusGUI
{
    internal class Parameter : INotifyPropertyChanged
    {
        private object _Value;
        private bool _Status;

        public Parameter()
        {
        }

        public Parameter(string name, string valueType)
        {
            ParamName = name;
            ValueType = valueType;

            ProtList = new ObservableCollection<Protease>();
            InitList = new ObservableCollection<string>();
            ProductMassToleranceList = new ObservableCollection<string>();

            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
            {
                ProtList.Add(protease);
            }
            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
            {
                InitList.Add(initiatior_methionine_behavior);
            }
            ProductMassToleranceList.Add("Absolute");
            ProductMassToleranceList.Add("ppm");
        }

        public event PropertyChangedEventHandler PropertyChanged;

        public string ParamName { get; set; }
        public string ValueType { get; set; }

        public object Value
        {
            get { return _Value; }
            set
            {
                _Value = value;
                OnPropertyChanged("Value");
            }
        }

        public bool Different { get; set; }

        public bool HasChanged
        {
            get { return _Status; }
            set
            {
                _Status = value;
                if (value == _Status) return;
                _Status = value;
            }
        }

        public ObservableCollection<Protease> ProtList { get; private set; }
        public ObservableCollection<string> InitList { get; private set; }
        public ObservableCollection<string> ProductMassToleranceList { get; private set; }

        protected virtual void OnPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(propertyName));
            }
            HasChanged = true;
        }
    }
}