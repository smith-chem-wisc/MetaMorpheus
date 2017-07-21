using EngineLayer.CrosslinkSearch;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkAnalysis
{
    public class CrosslinkAnalysisEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private const int max_mods_for_peptide = 3;
        private readonly int maximumMissedCleavages;
        private readonly int? minPeptideLength;
        private readonly int? maxPeptideLength;
        private readonly int maxModIsoforms;
        private readonly List<Tuple<PsmCross, PsmCross>> newPsms;
        private readonly List<Protein> proteinList;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly List<ModificationWithMass> fixedModifications;
        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;
        private readonly Protease protease;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private readonly Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans;
        private readonly Tolerance fragmentTolerance;
        private readonly List<ProductType> lp;
        private readonly InitiatorMethionineBehavior initiatorMethionineBehavior;

        //Draw Control
        private readonly CrosslinkerTypeClass crosslinker;

        private Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;
        private string OutputFolder;

        #endregion Private Fields

        //private bool draw = true;

        #region Public Constructors

        public CrosslinkAnalysisEngine(List<Tuple<PsmCross, PsmCross>> newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, Protease protease, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, Tolerance fragmentTolerance, int maximumMissedCleavages, int? minPeptideLength, int? maxPeptideLength, int maxModIsoforms, List<ProductType> lp, InitiatorMethionineBehavior initiatorMethionineBehavior, Dictionary<ModificationWithMass, ushort> modsDictionary, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, string OutputFolder, CrosslinkerTypeClass crosslinker, List<string> nestedIds) : base(nestedIds)
        {
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.protease = protease;

            this.myMsDataFile = myMSDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.minPeptideLength = minPeptideLength;
            this.maxPeptideLength = maxPeptideLength;
            this.maxModIsoforms = maxModIsoforms;
            this.lp = lp;
            this.initiatorMethionineBehavior = initiatorMethionineBehavior;
            this.modsDictionary = modsDictionary;
            this.arrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            this.OutputFolder = OutputFolder;
            this.crosslinker = crosslinker;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            MassDiffAcceptor XLsearchMode = new OpenSearchMode();

            CrosslinkAnalysisResults myAnalysisResults = new CrosslinkAnalysisResults(this);
            Status("Running analysis engine!", nestedIds);
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy

            #region Match Seqeunces to PeptideWithSetModifications

            //myAnalysisResults.AddText("Starting compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            Status("Adding observed peptides to dictionary...", nestedIds);
            foreach (var psmpair in newPsms)
            {
                if (psmpair != null)
                {
                    if (psmpair.Item1 != null)
                    {
                        var cp = psmpair.Item1.CompactPeptide;
                        if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                            compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                    }
                    if (psmpair.Item2 != null)
                    {
                        var cp = psmpair.Item2.CompactPeptide;
                        if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                            compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                    }
                }
            }
            //myAnalysisResults.AddText("Ending compactPeptideToProteinPeptideMatching count: " + compactPeptideToProteinPeptideMatching.Count);
            int totalProteins = proteinList.Count;
            int proteinsSeen = 0;
            int old_progress = 0;
            var obj = new object();
            Status("Adding possible sources to peptide dictionary...", nestedIds);
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());
                for (int i = fff.Item1; i < fff.Item2; i++)
                    foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications))
                    {
                        if (peptideWithPossibleModifications.Length <= 1)
                            continue;
                        foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModIsoforms, max_mods_for_peptide))
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications), out v))
                                v.Add(peptideWithSetModifications);
                        }
                    }
                lock (obj)
                {
                    foreach (var ye in local)
                    {
                        HashSet<PeptideWithSetModifications> v;
                        if (compactPeptideToProteinPeptideMatching.TryGetValue(ye.Key, out v))
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

            List<PsmCross> allResultingIdentifications = new List<PsmCross>();
            List<Tuple<PsmCross, PsmCross>> allResultingIdentificationsfdr = new List<Tuple<PsmCross, PsmCross>>();

            Status("Computing info about actual peptides with modifications...", nestedIds);
            for (int myScanWithMassIndex = 0; myScanWithMassIndex < newPsms.Count; myScanWithMassIndex++)
            {
                var huh = newPsms[myScanWithMassIndex].Item1;
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.ResolveProteinsAndMostProbablePeptide(compactPeptideToProteinPeptideMatching);
                var huh1 = newPsms[myScanWithMassIndex].Item2;
                if (huh1 != null && huh1.MostProbableProteinInfo == null)
                    huh1.ResolveProteinsAndMostProbablePeptide(compactPeptideToProteinPeptideMatching);
                newPsms[myScanWithMassIndex].Item1.XLTotalScore = newPsms[myScanWithMassIndex].Item1.XLBestScore + newPsms[myScanWithMassIndex].Item2.XLBestScore;
            }


            //Calculate Crosslink peptide FDR
            var CrosslinkOrderedPsmCrossWithPeptides = newPsms.OrderByDescending(b => b.Item1.XLTotalScore).ToList();
            var CrosslinkOrderedPsmsWithFDR = CrosslinkDoFalseDiscoveryRateAnalysis(CrosslinkOrderedPsmCrossWithPeptides, XLsearchMode);

            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        //Calculate the FDR of crosslinked peptide FP/(FP+TP)
        private static List<Tuple<PsmCross, PsmCross>> CrosslinkDoFalseDiscoveryRateAnalysis(List<Tuple<PsmCross, PsmCross>> items, MassDiffAcceptor sm)
        {
            var ids = new List<Tuple<PsmCross, PsmCross>>();
            foreach (var item in items)
            {
                ids.Add(new Tuple<PsmCross, PsmCross>(item.Item1, item.Item2));
            }

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches];

            for (int i = 0; i < ids.Count; i++)
            {
                var item1 = ids[i].Item1; var item2 = ids[i].Item2;

                var isDecoy1 = item1.MostProbableProteinInfo.IsDecoy; var isDecoy2 = item1.MostProbableProteinInfo.IsDecoy;
                int notch1 = item1.MostProbableProteinInfo.Notch; int notch2 = item1.MostProbableProteinInfo.Notch;
                if (isDecoy1 || isDecoy2)
                    cumulative_decoy++;
                else
                    cumulative_target++;


                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                item1.SetValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
                item2.SetValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
            }

            double min_q_value = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                PsmCross id = ids[i].Item1;
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;
            }

            return ids;
        }

        //Calculate n-score based on the equation from xlinkx
        private double XLCalculateNScore(int N, int n, int la, int lb, int ftotal, int ionType, double tolerance)
        {
            double x = 1 / 111.1 * ionType * tolerance * 2;
            double f = (double)lb / ((double)la + (double)lb) * ftotal;
            double e = Math.E;
            double p;
            double px = 0;
            double ifactorial = 1;
            for (int i = 0; i < n; i++)
            {
                if (i == 0) { ifactorial = 1; } else { ifactorial *= i; }
                px += Math.Pow(e, -x * f) * Math.Pow(x * f, i) / ifactorial;
            }
            p = 1 - px;
            return p * N;
        }

        #endregion Private Methods

    }
}