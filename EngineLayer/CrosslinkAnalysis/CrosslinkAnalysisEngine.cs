//using OxyPlot;
//using OxyPlot.Axes;
//using OxyPlot.Series;
//using OxyPlot.Annotations;
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
                    foreach (var peptideWithPossibleModifications in proteinList[i].Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications, false))
                    {
                        if (peptideWithPossibleModifications.Length <= 1)
                            continue;
                        foreach (var peptideWithSetModifications in peptideWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maxModIsoforms, max_mods_for_peptide))
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(peptideWithSetModifications, false), out v))
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

            #region Calculate single peptide FDR

            //List<PsmCross> PsmCrossForfdrList = new List<PsmCross>();
            //foreach (var PsmCrossForfdr in newPsms)
            //{
            //    //Filter matches that do not have improved XLLocalScore
            //    if (PsmCrossForfdr.Item1.XLBestScore > PsmCrossForfdr.Item1.Score)
            //    {
            //        PsmCrossForfdrList.Add(PsmCrossForfdr.Item1);
            //        PsmCrossForfdrList.Add(PsmCrossForfdr.Item2);
            //    }
            //}
            //Status("Sorting and grouping psms..", nestedIds);
            //var orderedPsmsWithPeptides = PsmCrossForfdrList.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.Pli.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, string>(b.FullFilePath, b.ScanNumber, b.Pli.FullSequence)).Select(b => b.First());

            //Status("Running FDR analysis...", nestedIds);
            //var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides, XLsearchMode);
            //allResultingIdentifications = orderedPsmsWithFDR.OrderBy(p => p.ScanNumber).ToList();
            //myAnalysisResults.AllResultingIdentifications = allResultingIdentifications;

            #endregion Calculate single peptide FDR

            //Calculate Crosslink peptide FDR
            var CrosslinkOrderedPsmCrossWithPeptides = newPsms.OrderByDescending(b => b.Item1.XLTotalScore).ToList();
            var CrosslinkOrderedPsmsWithFDR = CrosslinkDoFalseDiscoveryRateAnalysis(CrosslinkOrderedPsmCrossWithPeptides, XLsearchMode);

            //Filter out crosslink peptide from Decoy and Qvalue > 0.01
            //allResultingIdentificationsfdr = CrosslinkOrderedPsmsWithFDR.Where(p => p.Item1.Pli.IsDecoy != true && p.Item2.Pli.IsDecoy != true && p.Item1.FdrInfo.QValue <= 0.01).ToList();

            #region Draw spetra anotation

            //if (draw)
            //{
            //    foreach (var drawCompactPepList in CrosslinkOrderedPsmsWithFDR)
            //    {
            //        XLDrawMSMatchToPdf(arrayOfSortedMS2Scans?[drawCompactPepList.Item1.ScanIndex], drawCompactPepList);
            //    }
            //}

            #endregion Draw spetra anotation

            //myAnalysisResults.AllResultingIdentificationFdrPairs = allResultingIdentificationsfdr;
            return myAnalysisResults;
        }

        #endregion Protected Methods

        /* DoFalseDiscoveryRateAnalysis
        private static List<PsmCross> DoFalseDiscoveryRateAnalysis(IEnumerable<PsmCross> items, MassDiffAcceptor sm)
        {
            var ids = new List<PsmCross>();
            foreach (PsmCross item in items)
                ids.Add(item);

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches];

            for (int i = 0; i < ids.Count; i++)
            {
                var item = ids[i];
                var isDecoy = item.Pli.IsDecoy;
                int notch = item.Notch;
                if (isDecoy)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                if (isDecoy)
                    cumulative_decoy_per_notch[notch]++;
                else
                    cumulative_target_per_notch[notch]++;

                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                double temp_q_value_for_notch = (double)cumulative_decoy_per_notch[notch] / (cumulative_target_per_notch[notch] + cumulative_decoy_per_notch[notch]);
                item.SetValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch);
            }

            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[sm.NumNotches];
            for (int i = 0; i < sm.NumNotches; i++)
                min_q_value_notch[i] = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                PsmCross id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;

                int notch = id.Notch;
                if (id.FdrInfo.QValueNotch > min_q_value_notch[notch])
                    id.FdrInfo.QValueNotch = min_q_value_notch[notch];
                else if (id.FdrInfo.QValueNotch < min_q_value_notch[notch])
                    min_q_value_notch[notch] = id.FdrInfo.QValueNotch;
            }

            return ids;
        }
        */

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

                //if (isDecoy1 || isDecoy2)
                //    cumulative_decoy_per_notch[notch]++;
                //else
                //    cumulative_target_per_notch[notch]++;

                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                //double temp_q_value_for_notch = (double)cumulative_decoy_per_notch[notch] / (cumulative_target_per_notch[notch] + cumulative_decoy_per_notch[notch]);
                //item1.SetValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch);
                //item2.SetValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch);
                item1.SetValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
                item2.SetValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
            }

            double min_q_value = double.PositiveInfinity;
            //double[] min_q_value_notch = new double[sm.NumNotches];
            //for (int i = 0; i < sm.NumNotches; i++)
            //    min_q_value_notch[i] = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                PsmCross id = ids[i].Item1;
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;

                //int notch = id.thisPSM.Notch;
                //if (id.QValueNotch > min_q_value_notch[notch])
                //    id.QValueNotch = min_q_value_notch[notch];
                //else if (id.QValueNotch < min_q_value_notch[notch])
                //    min_q_value_notch[notch] = id.QValueNotch;
            }

            return ids;
        }

        /* Draw spectra annotation
        private void XLDrawMSMatchToPdf(Ms2ScanWithSpecificMass MsScanForDraw, Tuple<NewPsmWithFdr, NewPsmWithFdr> PsmCrosssForDraw)
        {
            var x = MsScanForDraw.TheScan.MassSpectrum.XArray;
            var y = MsScanForDraw.TheScan.MassSpectrum.YArray;

            string scanNum = MsScanForDraw.TheScan.OneBasedScanNumber.ToString();
            string sequence1 = PsmCrosssForDraw.Item1.thisPSM.Pli.FullSequence;
            string sequence2 = PsmCrosssForDraw.Item2.thisPSM.Pli.FullSequence;

            var matchedIonDic1 = PsmCrosssForDraw.Item1.thisPSM.matchedIonInfo;
            var matchedIonDic2 = PsmCrosssForDraw.Item2.thisPSM.matchedIonInfo;

            PlotModel model = new PlotModel { Title = "Spectrum anotation of Scan " + scanNum + " for Crosslinked Peptide", DefaultFontSize = 15 };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = x.Max() * 1.02 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity(counts)", Minimum = 0, Maximum = y.Max() * 1.2 });
            var textAnnoSeq1 = new TextAnnotation() { };
            textAnnoSeq1.FontSize = 9; textAnnoSeq1.TextColor = OxyColors.Red; textAnnoSeq1.StrokeThickness = 0; textAnnoSeq1.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.15); textAnnoSeq1.Text = sequence1;
            var textAnnoSeq2 = new TextAnnotation() { };
            textAnnoSeq2.FontSize = 9; textAnnoSeq2.TextColor = OxyColors.Blue; textAnnoSeq2.StrokeThickness = 0; textAnnoSeq2.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.1); textAnnoSeq2.Text = sequence2;
            model.Annotations.Add(textAnnoSeq1);
            model.Annotations.Add(textAnnoSeq2);

            LineSeries[] s0 = new LineSeries[x.Length];
            LineSeries[] s1 = new LineSeries[x.Length];
            LineSeries[] s2 = new LineSeries[x.Length];

            //Draw the ms/ms scan peaks
            for (int i = 0; i < x.Length; i++)
            {
                s0[i] = new LineSeries();
                s0[i].Color = OxyColors.DimGray;
                s0[i].StrokeThickness = 0.15;
                s0[i].Points.Add(new DataPoint(x[i], 0));
                s0[i].Points.Add(new DataPoint(x[i], y[i]));
                model.Series.Add(s0[i]);
            }
            //Draw the ms/ms scan matched peaks

            //OxyColor ionColor = OxyColors.White;

            //if (t == ProductType.B) { ionColor = OxyColors.OrangeRed; }
            //if (t == ProductType.Y) { ionColor = OxyColors.Blue; }  //DeepSkyBlue
            //if (t == ProductType.C) { ionColor = OxyColors.Crimson; }
            //if (t == ProductType.Zdot) { ionColor = OxyColors.Teal; }

            for (int i = 0; i < matchedIonDic1.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Red;
                if (matchedIonDic1.MatchedIonMz[i] > 0)
                {
                    s1[i] = new LineSeries();
                    s1[i].Color = ionColor;
                    s1[i].StrokeThickness = 0.2;
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i], 0));
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i], matchedIonDic1.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    textAnno1.FontSize = 6;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s1[i].Points[1];
                    textAnno1.Text = matchedIonDic1.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    textAnno2.FontSize = 6;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s1[i].Points[1].X, s1[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic1.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s1[i]);
                }
            }

            for (int i = 0; i < matchedIonDic2.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Blue;
                if (matchedIonDic2.MatchedIonMz[i] > 0)
                {
                    s2[i] = new LineSeries();
                    s2[i].Color = ionColor;
                    s2[i].StrokeThickness = 0.2;
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i], 0));
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i], matchedIonDic2.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    textAnno1.FontSize = 6;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s2[i].Points[1];
                    textAnno1.Text = matchedIonDic2.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    textAnno2.FontSize = 6;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s2[i].Points[1].X, s2[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic2.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s2[i]);
                }
            }

            using (var stream = File.Create(OutputFolder + "\\" + "Scan" + scanNum + ".pdf"))
            {
                PdfExporter pdf = new PdfExporter { Width = 700, Height = 300 };
                pdf.Export(model, stream);
            }
        }
        */

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