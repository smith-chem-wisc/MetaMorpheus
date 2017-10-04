using EngineLayer;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    class Parameter : INotifyPropertyChanged
    {
        public string ParamName { get; set; }

        public string ValueType { get; set; }

        object _value;
        public object Value
        {
            get { return _value; }
            set
            {
                _value = value;
               OnPropertyChanged("Value");
            }
        }
       public event PropertyChangedEventHandler PropertyChanged;

        public bool Different { get; set; }

        private bool status;

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

            foreach (Protease protease in GlobalEngineLevelSettings.ProteaseDictionary.Values)
                ProtList.Add(protease);
            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                InitList.Add(initiatior_methionine_behavior);
            ProductMassToleranceList.Add("Absolute");
            ProductMassToleranceList.Add("ppm");

        }

         protected virtual void OnPropertyChanged(string propertyName)
         {
             PropertyChangedEventHandler handler = PropertyChanged;
             if (handler != null)
                 handler(this, new PropertyChangedEventArgs(propertyName));

             this.HasChanged = true;
         }

    }


}
