using MetaMorpheus;
using System;
using System.Collections.ObjectModel;

namespace GoodGUI
{
    public class MyGPTMDtask : MyTask
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

        internal override void DoTask(ObservableCollection<RawData> completeRawFileListCollection, ObservableCollection<XMLdb> completeXmlDbList, AllTasksParams po)
        {
            throw new NotImplementedException();
        }
    }
}