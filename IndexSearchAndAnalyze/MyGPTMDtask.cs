using MetaMorpheus;
using System;
using System.Collections.ObjectModel;

namespace IndexSearchAndAnalyze
{
    public class MyGPTMDtask : MyTask
    {
        public MyGPTMDtask() : base(2)
        {
        }

        public override void DoTask(ObservableCollection<RawData> completeRawFileListCollection, ObservableCollection<XMLdb> completeXmlDbList, AllTasksParams po)
        {
            throw new NotImplementedException();
        }
    }
}