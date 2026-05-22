using System.Collections.Generic;
using EngineLayer;
using EngineLayer.DatabaseLoading;

namespace TaskLayer
{
    /// <summary>
    /// Identifies N- and C-terminally truncated proteoforms by re-searching MS2 scans against a
    /// fragment index built from the proteoforms found by an upstream top-down <see cref="SearchTask"/>.
    /// See 01_Architecture.md for the full three-pass design (index + dual single-series scoring,
    /// terminus-directed chopping, pooled FDR/PEP).
    ///
    /// Phase 0 scaffolding stub: <see cref="RunSpecific"/> returns an empty result. The algorithm
    /// is implemented in Phases 1-3.
    /// </summary>
    public class TruncationSearchTask : MetaMorpheusTask
    {
        public TruncationSearchTask() : base(MyTask.Truncation)
        {
            CommonParameters = new CommonParameters();
            TruncationSearchParameters = new TruncationSearchParameters();
        }

        public TruncationSearchParameters TruncationSearchParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList,
            List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            // TODO Phases 1-3: load the deduped Pass 1 proteoform list (from TaskChainContext, or the
            // AllProteoforms.psmtsv disk fallback), build the Pass 2 fragment index, run dual
            // single-series scoring + terminus-directed chopping, then pooled FDR/PEP and write
            // AllTruncatedPSMs.psmtsv / AllTruncatedProteoforms.psmtsv.
            MyTaskResults = new MyTaskResults(this);
            return MyTaskResults;
        }
    }
}
