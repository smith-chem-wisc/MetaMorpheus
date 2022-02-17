using EngineLayer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class SingleCellProteomicsTask : MetaMorpheusTask
    {
        public SingleCellProteomicsTask() : base(MyTask.SingleCellProteomics)
        {
            CommonParameters = new CommonParameters();
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            throw new NotImplementedException();
        }
    }
}
