using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.TmtPostProcessing
{
    public class TmtPostProcessingResults : MetaMorpheusEngineResults
    {
        public TmtPostProcessingResults(MetaMorpheusEngine s, object analysisType) : base(s)
        {
        }
        public override string ToString()
        {
            var sb = new StringBuilder();
            //sb.AppendLine(base.ToString());
            //sb.AppendLine($"{AnalysisType}s within 1% FDR: {PsmsWithin1PercentFdr.ToString()}");
            //sb.AppendLine($"Delta Score Used for FDR Analysis: {DeltaScoreImprovement.ToString()}");
            //sb.AppendLine(BinarySearchTreeMetrics);
            return sb.ToString();
        }
    }
}
