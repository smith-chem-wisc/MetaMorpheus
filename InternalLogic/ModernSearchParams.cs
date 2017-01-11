using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System.Collections.Generic;

namespace InternalLogic
{
    public class ModernSearchParams : MyParams
    {
        public List<MorpheusModification> fixedModifications { get; private set; }
        public List<int>[] fragmentIndex { get; private set; }
        public double fragmentTolerance { get; private set; }
        public float[] keys { get; private set; }
        public List<MorpheusModification> localizeableModifications { get; private set; }
        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile { get; private set; }
        public List<CompactPeptide> peptideIndex { get; private set; }
        public Protease protease { get; private set; }
        public List<Protein> proteinList { get; private set; }
        public List<SearchMode> searchModes { get; private set; }
        public int spectraFileIndex { get; private set; }
        public List<MorpheusModification> variableModifications { get; private set; }

        public ModernSearchParams(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int spectraFileIndex, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, List<Protein> proteinList, double fragmentTolerance, Protease protease, List<SearchMode> searchModes, AllTasksParams a2) : base(a2)
        {
            this.myMsDataFile = myMsDataFile;
            this.spectraFileIndex = spectraFileIndex;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.proteinList = proteinList;
            this.fragmentTolerance = fragmentTolerance;
            this.protease = protease;
            this.searchModes = searchModes;
        }

        internal override void Validate()
        {
        }
    }
}