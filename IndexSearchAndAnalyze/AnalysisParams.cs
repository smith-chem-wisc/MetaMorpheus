using MassSpectrometry;
using MetaMorpheus;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;

namespace IndexSearchAndAnalyze
{
    public class AnalysisParams : MyParams
    {
        public List<NewPsm>[] newPsms { get; private set; }
        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching { get; private set; }
        public List<Protein> proteinList { get; private set; }
        public List<MorpheusModification> variableModifications { get; private set; }
        public List<MorpheusModification> fixedModifications { get; private set; }
        public List<MorpheusModification> localizeableModifications { get; private set; }
        public Protease protease { get; private set; }
        public List<SearchMode> searchModes { get; private set; }
        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile { get; private set; }
        public double fragmentTolerance { get; private set; }
        public UsefulProteomicsDatabases.Generated.unimod unimodDeserialized { get; internal set; }
        public Dictionary<int, ChemicalFormulaModification> uniprotDeseralized { get; private set; }
        public Action<MyNewTreeStructure, string> action1 { get; private set; }
        public Action<List<NewPsmWithFDR>, string> action2 { get; private set; }

        public AnalysisParams(List<NewPsm>[] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, double fragmentTolerance, UsefulProteomicsDatabases.Generated.unimod unimodDeserialized, Dictionary<int, ChemicalFormulaModification> uniprotDeseralized, Action<MyNewTreeStructure, string> action1, Action<List<NewPsmWithFDR>, string> action2)
        {
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.protease = protease;
            this.searchModes = searchModes;
            this.myMsDataFile = myMsDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.unimodDeserialized = unimodDeserialized;
            this.uniprotDeseralized = uniprotDeseralized;
            this.action1 = action1;
            this.action2 = action2;
        }
    }
}
