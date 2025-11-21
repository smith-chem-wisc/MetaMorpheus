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
        private List<SpectralMatch> AllPsms;
        private TmtProcessingParams Tmtparams;
        private readonly string AnalysisType;
        public TmtPostProcessingEngine(List<SpectralMatch> psms, TmtProcessingParams tmtParams, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIdsstring, string analysisType = "PSM") : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllPsms = psms.OrderByDescending(p => p).ToList();
            AnalysisType = analysisType;
            Tmtparams = tmtParams;
            if (fileSpecificParameters == null) throw new ArgumentNullException("file specific parameters cannot be null");
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

        private void ProcessTmtData(TmtPostProcessingResults myAnalysisResults) => ProcessTmtData(AllPsms, Tmtparams, myAnalysisResults);
        internal static void ProcessTmtData(List<SpectralMatch> allPsms, TmtProcessingParams tmtParams, TmtPostProcessingResults myAnalysisResults)
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }
            
            throw new NotImplementedException();
        }
    }
}
