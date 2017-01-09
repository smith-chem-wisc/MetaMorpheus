using MassSpectrometry;
using MetaMorpheus;
using Proteomics;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using UsefulProteomicsDatabases.Generated;

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
        public Tolerance fragmentTolerance { get; private set; }
        public unimod unimodDeserialized { get; internal set; }
        public Dictionary<int, ChemicalFormulaModification> uniprotDeseralized { get; private set; }
        public Action<BinTreeStructure, string> action1 { get; private set; }
        public Action<List<NewPsmWithFDR>, string> action2 { get; private set; }

        public AnalysisParams(List<NewPsm>[] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFDR>, string> action2, AllTasksParams a2) : base(a2)
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
            this.action1 = action1;
            this.action2 = action2;
            this.unimodDeserialized = AllTasksParams.unimodDeserialized;
            this.uniprotDeseralized = AllTasksParams.uniprotDeseralized;
        }

        // For a single search mode
        public AnalysisParams(List<NewPsm> newPsms1, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching1, List<Protein> proteinList1, List<MorpheusModification> variableModifications1, List<MorpheusModification> fixedModifications1, List<MorpheusModification> localizeableModifications1, Protease protease1, SearchMode searchMode, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile1, Tolerance fragmentTolerance, Action<BinTreeStructure, string> p1, Action<List<NewPsmWithFDR>, string> p2, AllTasksParams a2) : base(a2)
        {
            this.newPsms = new List<NewPsm>[1] { newPsms1 };
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching1;
            this.proteinList = proteinList1;
            this.variableModifications = variableModifications1;
            this.fixedModifications = fixedModifications1;
            this.localizeableModifications = localizeableModifications1;
            this.protease = protease1;
            this.searchModes = new List<SearchMode> { searchMode };
            this.myMsDataFile = myMsDataFile1;
            this.fragmentTolerance = fragmentTolerance;
            this.action1 = p1;
            this.action2 = p2;
        }

        internal override void Validate()
        {
            if (newPsms == null)
                throw new ValidationException("newPsms is null");
            if (proteinList == null)
                throw new ValidationException("proteinList is null");
        }
    }
}