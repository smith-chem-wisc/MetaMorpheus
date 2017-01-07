using MassSpectrometry;
using MetaMorpheus;
using Spectra;
using System;
using System.Collections.Generic;

namespace IndexSearchAndAnalyze
{
    public class ClassicSearchParams : MyParams
    {
        public List<Protein> proteinList { get; private set; }
        public Protease protease { get; private set; }
        public List<MorpheusModification> fixedModifications { get; private set; }
        public List<MorpheusModification> localizeableModifications { get; private set; }
        public List<MorpheusModification> variableModifications { get; private set; }
        public Tolerance productMassTolerance { get; private set; }
        public SearchMode searchMode { get; private set; }
        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile { get; private set; }
        public int spectraFileIndex { get; private set; }

        public ClassicSearchParams(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int spectraFileIndex, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, List<Protein> proteinList, Tolerance fragmentTolerance, Protease protease, SearchMode searchMode, Action<string> a1, Action<int> a2) : base(a1, a2)
        {
            this.myMsDataFile = myMsDataFile;
            this.spectraFileIndex = spectraFileIndex;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.proteinList = proteinList;
            this.productMassTolerance = fragmentTolerance;
            this.protease = protease;
            this.searchMode = searchMode;
        }

        internal override void Validate()
        {
        }
    }
}