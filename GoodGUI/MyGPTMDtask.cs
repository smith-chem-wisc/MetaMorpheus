using MetaMorpheus;
using System.Collections.ObjectModel;
using System.Linq;
using System;

namespace GoodGUI
{
    internal class MyGPTMDtask : MyTask
    {
        public MyGPTMDtask() : base(2)
        {
        }
        public string precursorMassToleranceTextBox { get; internal set; }
        public int precursorMassToleranceComboBox { get; internal set; }
        public string missedCleavagesTextBox { get; internal set; }
        public Protease protease { get; internal set; }
        public string maxModificationIsoformsTextBox { get; internal set; }
        public int initiatorMethionineBehaviorComboBox { get; internal set; }
        public string productMassToleranceTextBox { get; internal set; }
        public int productMassToleranceComboBox { get; internal set; }
        public bool? bCheckBox { get; internal set; }
        public bool? yCheckBox { get; internal set; }
        public bool? checkBoxDecoy { get; internal set; }
        public string acceptedPrecursorMassErrorsTextBox { get; internal set; }
        public bool? checkBoxMonoisotopic { get; internal set; }

        internal override void DoTask(ObservableCollection<RawData> completeRawFileListCollection, ObservableCollection<XMLdb> completeXmlDbList, ParamsObject po)
        {
            throw new NotImplementedException();
        }
    }
}