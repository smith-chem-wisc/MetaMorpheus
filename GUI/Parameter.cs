using EngineLayer;
using System;
using System.Collections.ObjectModel;
using System.ComponentModel;

namespace MetaMorpheusGUI
{
    internal class Parameter : INotifyPropertyChanged
    {
        #region Private Fields

        private object _value;
        private bool status;

        #endregion Private Fields

        #region Public Constructors

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

            foreach (Protease protease in GlobalVariables.ProteaseDictionary.Values)
                ProtList.Add(protease);
            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                InitList.Add(initiatior_methionine_behavior);
            ProductMassToleranceList.Add("Absolute");
            ProductMassToleranceList.Add("ppm");
        }

        #endregion Public Constructors

        #region Public Events

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion Public Events

        #region Public Properties

        public string ParamName { get; set; }

        public string ValueType { get; set; }

        public object Value
        {
            get { return _value; }
            set
            {
                _value = value;
                OnPropertyChanged("Value");
            }
        }

        public bool Different { get; set; }

        public bool HasChanged
        {
            get { return status; }
            set
            {
                status = value;
                if (value == status) return;
                status = value;
            }
        }

        public ObservableCollection<Protease> ProtList { get; private set; }

        public ObservableCollection<string> InitList { get; private set; }

        public ObservableCollection<string> ProductMassToleranceList { get; private set; }

        #endregion Public Properties

        #region Protected Methods

        protected virtual void OnPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
                handler(this, new PropertyChangedEventArgs(propertyName));

            this.HasChanged = true;
        }

        #endregion Protected Methods
    }
}