using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibSpectralAveraging;

namespace TaskLayer
{
    public class SpectralAveragingTask : MetaMorpheusTask
    {
        public SpectralAveragingTask(MyTask taskType) : base(taskType)
        {
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            throw new NotImplementedException();
        }
    }
}
