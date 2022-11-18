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
        public MzLibSpectralAveragingOptions Options { get; set; }

        public SpectralAveragingTask(MzLibSpectralAveragingOptions options) : base(MyTask.Average)
        {
            Options = options;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            throw new NotImplementedException();
        }
    }
}
