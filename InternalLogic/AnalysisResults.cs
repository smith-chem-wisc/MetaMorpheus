using OldInternalLogic;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class AnalysisResults : MyResults
    {
        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> dict { get; private set; }

        public List<NewPsmWithFDR>[] allResultingIdentifications { get; private set; }

        public AnalysisResults(AnalysisEngine s, List<NewPsmWithFDR>[] allResultingIdentifications, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> dict) : base(s)
        {
            this.allResultingIdentifications = allResultingIdentifications;
            this.dict = dict;
        }

        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("\t\tAll PSMS within 1% FDR: " + string.Join(", ", allResultingIdentifications.Select(b => b.Where(c => c.QValue <= 0.01).Count())));
            return sb.ToString();
        }
    }
}