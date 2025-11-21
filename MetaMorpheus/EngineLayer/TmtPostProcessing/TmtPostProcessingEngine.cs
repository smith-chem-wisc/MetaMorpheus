using EngineLayer.FdrAnalysis;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.TmtPostProcessing
{
    public class TmtPostProcessingEngine : MetaMorpheusEngine
    {
        private readonly string AnalysisType;
        public TmtPostProcessingEngine(CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            TmtPostProcessingResults myAnalysisResults = new TmtPostProcessingResults(this, AnalysisType);

            Status("Running TMT Post-processing...");
            ProcessTmtData(myAnalysisResults);
            Status("Done.");
            //myAnalysisResults.PsmsWithin1PercentFdr = AllPsms.Count(b => b.FdrInfo.QValue <= 0.01 && !b.IsDecoy);

            return myAnalysisResults;
        }

        private void ProcessTmtData(TmtPostProcessingResults myAnalysisResults)
        {
            throw new NotImplementedException();
        }
    }
}
