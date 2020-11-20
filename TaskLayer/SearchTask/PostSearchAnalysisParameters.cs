using EngineLayer;
using FlashLFQ;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;

namespace TaskLayer
{
    public class PostSearchAnalysisParameters
    {
        public MyTaskResults SearchTaskResults { get; set; }
        public string SearchTaskId { get; set; }
        public SearchParameters SearchParameters { get; set; }
        public List<Protein> ProteinList { get; set; }
        public List<Modification> VariableModifications { get; set; }
        public List<Modification> FixedModifications { get; set; }
        public HashSet<DigestionParams> ListOfDigestionParams { get; set; }
        public List<PeptideSpectralMatch> AllPsms { get; set; }
        public FlashLfqResults FlashLfqResults { get; set; }
        public int NumNotches { get; set; }
        public string OutputFolder { get; set; }
        public string IndividualResultsOutputFolder { get; set; }
        public FileSpecificParameters[] FileSettingsList { get; set; }
        public Dictionary<string, int[]> NumMs2SpectraPerFile { get; set; }
        public MyFileManager MyFileManager { get; set; }
        public List<DbForTask> DatabaseFilenameList { get; set; }
        public List<string> CurrentRawFileList { get; set; }
    }
}