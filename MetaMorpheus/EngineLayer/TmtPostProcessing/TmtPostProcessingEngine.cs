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
        /// <summary>
        /// Protein parsimony must be completed prior to TMT post-processing and correct protein IDs must be present on PSMs
        /// </summary>
        /// <param name="psms"></param>
        /// <param name="tmtParams"></param>
        /// <param name="commonParameters"></param>
        /// <param name="fileSpecificParameters"></param>
        /// <param name="nestedIdsstring"></param>
        /// <param name="analysisType"></param>
        /// <exception cref="ArgumentNullException"></exception>
        public TmtPostProcessingEngine(List<SpectralMatch> psms, TmtProcessingParams tmtParams, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIdsstring, string analysisType = "PSM") : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllPsms = psms.Where(p=>p != null && p.IsobaricMassTagReporterIonIntensities.Any()).ToList();
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
            ConstructIntensityTable();
            EliminateLowIntensityPsms();
            NormalizeLoading();
            NormalizeToReference();
            SelectPsmForProteinQuant();
            throw new NotImplementedException();
        }
        private static void ConstructIntensityTable()
        {
                       // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }
        }
        private static void EliminateLowIntensityPsms()
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }

        }
        private static void NormalizeLoading()
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }

        }
        private static void NormalizeToReference()
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }
        }
        private static void SelectPsmForProteinQuant() 
        { 
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }
        }
        /// <summary>
        /// PAW Pipeline sums PSMs for peptide quantification
        /// </summary>
        private static void SumPsmsForPeptideQuant()
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }
        }
        /// <summary>
        /// PAW Pipeline sums peptides for protein quantification
        /// </summary>
        private static void SumPeptidesForProteinQuant()
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }
        }
    }
}
