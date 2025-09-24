using EngineLayer.GlycoSearch;
using FlashLFQ;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using Omics.Modifications;

namespace TaskLayer
{
    public class PostGlycoSearchAnalysisParameters
    {
        public MyTaskResults GlycoSearchTaskResults { get; set; }
        public string SearchTaskId { get; set; }
        public HashSet<DigestionParams> ListOfDigestionParams { get; set; }
        public GlycoSearchParameters GlycoSearchParameters { get; set; }
        public List<Protein> ProteinList { get; set; }
        public List<Modification> VariableModifications { get; set; }
        public List<Modification> FixedModifications { get; set; }
        public List<GlycoSpectralMatch> AllPsms { get; set; }
        public string OutputFolder { get; set; }
        public FileSpecificParameters[] FileSettingsList { get; set; }
        public List<DbForTask> DatabaseFilenameList { get; set; }
        public List<string> CurrentRawFileList { get; set; }
        public FlashLfqResults FlashLfqResults { get; set; }
        public string IndividualResultsOutputFolder { get; set; }
        public List<Protein> BioPolymerList { get; set; }
    }
}
