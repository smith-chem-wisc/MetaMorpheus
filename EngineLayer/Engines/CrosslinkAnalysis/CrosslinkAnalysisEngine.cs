using EngineLayer.CrosslinkSearch;
using Proteomics;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkAnalysis
{
    public class CrosslinkAnalysisEngine : MetaMorpheusEngine
    {
        #region Protected Fields

        protected readonly TerminusType terminusType;

        #endregion Protected Fields

        #region Private Fields

        private readonly List<PsmCross> newPsms;
        private readonly List<Protein> proteinList;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly List<ModificationWithMass> fixedModifications;
        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;

        private readonly List<ProductType> lp;
        private readonly CrosslinkerTypeClass crosslinker;

        private readonly Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;
        private readonly string OutputFolder;
        private readonly ICommonParameters CommonParameters;

        #endregion Private Fields

        #region Public Constructors

        public CrosslinkAnalysisEngine(List<PsmCross> newPsms, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ProductType> lp, Dictionary<ModificationWithMass, ushort> modsDictionary, string OutputFolder, CrosslinkerTypeClass crosslinker, TerminusType terminusType, ICommonParameters CommonParameters, List<string> nestedIds) : base(nestedIds)
        {
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.lp = lp;
            this.modsDictionary = modsDictionary;
            this.OutputFolder = OutputFolder;
            this.crosslinker = crosslinker;
            this.terminusType = terminusType;
            this.CommonParameters = CommonParameters;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            MassDiffAcceptor XLsearchMode = new OpenSearchMode();

            CrosslinkAnalysisResults myAnalysisResults = new CrosslinkAnalysisResults(this);
            Status("Running analysis engine!");
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy

            #region Match Seqeunces to PeptideWithSetModifications

            //myAnalysisResults.AddText("Starting compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            Status("Adding observed peptides to dictionary...");
            foreach (var psmpair in newPsms)
            {
                if (psmpair != null)
                {
                    var cp = psmpair.compactPeptide;
                    if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                        compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());

                    if (psmpair.BetaPsmCross != null)
                    {
                        var cp1 = psmpair.BetaPsmCross.compactPeptide;
                        if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp1))
                            compactPeptideToProteinPeptideMatching.Add(cp1, new HashSet<PeptideWithSetModifications>());
                    }
                }
            }
            //myAnalysisResults.AddText("Ending compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            int totalProteins = proteinList.Count;
            int proteinsSeen = 0;
            int old_progress = 0;
            var obj = new object();
            Status("Adding possible sources to peptide dictionary...");
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                for (int i = fff.Item1; i < fff.Item2; i++)
                    foreach (var peptideWithSetModifications in proteinList[i].Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList())
                        if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications, terminusType), out HashSet<PeptideWithSetModifications> v))
                            v.Add(peptideWithSetModifications);

                lock (obj)
                {
                    foreach (var ye in local)
                    {
                        if (compactPeptideToProteinPeptideMatching.TryGetValue(ye.Key, out HashSet<PeptideWithSetModifications> v))
                            foreach (var huh in ye.Value)
                                v.Add(huh);
                    }
                    proteinsSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In adding possible" +
                            " sources to peptide dictionary loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });

            #endregion Match Seqeunces to PeptideWithSetModifications

            Status("Computing info about actual peptides with modifications...");
            for (int myScanWithMassIndex = 0; myScanWithMassIndex < newPsms.Count; myScanWithMassIndex++)
            {
                var huh = newPsms[myScanWithMassIndex];
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
                if (huh != null)
                {
                    var huh1 = newPsms[myScanWithMassIndex].BetaPsmCross;
                    if (huh1 != null)
                        huh1.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
                }
            }

            return myAnalysisResults;
        }

        #endregion Protected Methods
    }
}