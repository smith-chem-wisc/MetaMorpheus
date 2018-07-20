using EngineLayer;
using FlashLFQ;
using MassSpectrometry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;

namespace TaskLayer
{
    public class PostSearchAnalysisParameters
    {
        public MyTaskResults SearchTaskResults { get; set; }
        public string SearchTaskId { get; set; }
        public CommonParameters CommonParameters { get; set; }
        public SearchParameters SearchParameters { get; set; }
        public List<Protein> ProteinList { get; set; }
        public List<ProductType> IonTypes { get; set; }
        public List<ModificationWithMass> VariableModifications { get; set; }
        public HashSet<DigestionParams> ListOfDigestionParams { get; set; }
        public List<PeptideSpectralMatch> AllPsms { get; set; }
        public FlashLFQResults FlashLfqResults { get; set; }
        public List<ModificationWithMass> FixedModifications { get; set; }
        public int NumNotches { get; set; }
        public string OutputFolder { get; set; }
        public FileSpecificParameters[] FileSettingsList { get; set; }
        public Dictionary<string, int[]> NumMs2SpectraPerFile { get; set; }
        public MyFileManager MyFileManager { get; set; }
        public List<DbForTask> DatabaseFilenameList { get; set; }
        public List<string> CurrentRawFileList { get; set; }
    }
}