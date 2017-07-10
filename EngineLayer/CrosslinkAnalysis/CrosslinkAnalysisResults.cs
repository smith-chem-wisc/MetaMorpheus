using System;
using System.Collections.Generic;
using EngineLayer.CrosslinkSearch;

namespace EngineLayer.CrosslinkAnalysis
{
    public class CrosslinkAnalysisResults : MetaMorpheusEngineResults
    {

        #region Internal Fields

        internal Dictionary<string, int> allModsOnPeptides;
        internal Dictionary<string, int> allModsSeen;

        #endregion Internal Fields

        #region Public Constructors

        public CrosslinkAnalysisResults(CrosslinkAnalysisEngine s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public List<Tuple<PsmCross, PsmCross>> AllResultingIdentificationFdrPairs { get; set; }

        #endregion Public Properties

    }
}
