using EngineLayer;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    class Parameter
    {
        public string ParamName { get; set; }

        public string ValueType { get; set; }

        public object Value { get; set; }

        public bool Different { get; set; }
        
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
            ProductMassToleranceList.Add("Ppm");

        }

    }


}
