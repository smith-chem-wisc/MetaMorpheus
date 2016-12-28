using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MetaMorpheus;
using Spectra;
using MassSpectrometry;

namespace IndexSearchAndAnalyze
{
    public class AnalysisParams : MyParams
    {
        public NewPsm[][] newPsms { get; private set; }
        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching { get; private set; }
        public List<Protein> proteinList { get; private set; }
        public List<MorpheusModification> variableModifications { get; private set; }
        public List<MorpheusModification> fixedModifications { get; private set; }
        public List<MorpheusModification> localizeableModifications { get; private set; }
        public Protease protease { get; private set; }
        public List<SearchMode> searchModes { get; private set; }
        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile { get; private set; }
        public double fragmentTolerance { get; private set; }

        public AnalysisParams(NewPsm[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, double fragmentTolerance)
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
        }

        internal void onfinished1(MyNewTreeStructure hm)
        {
            throw new NotImplementedException();
        }

        internal void onfinished2(List<NewPsmWithFDR> orderedPsmsWithFDR)
        {
            throw new NotImplementedException();
        }
    }
}
