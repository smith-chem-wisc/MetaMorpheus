using MassSpectrometry;
using OldInternalLogic;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using UsefulProteomicsDatabases.Generated;

namespace InternalLogic
{
    public class AnalysisParams : MyParams
    {
        public ParentSpectrumMatch[][] newPsms { get; private set; }
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
        public bool doParsimony { get; internal set; }

        public AnalysisParams(ParentSpectrumMatch[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFDR>, string> action2, AllTasksParams a2, bool doParsimony) : base(a2)
        {
            this.doParsimony = doParsimony;
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

        internal override void Validate()
        {
            if (newPsms == null)
                throw new ValidationException("newPsms is null");
            if (proteinList == null)
                throw new ValidationException("proteinList is null");
        }
    }
}