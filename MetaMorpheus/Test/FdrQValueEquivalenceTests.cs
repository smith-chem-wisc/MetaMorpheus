using Chemistry;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// CalculateQValue now counts each match's hypotheses per notch without LINQ and without the sort that
    /// BestMatchingBioPolymersWithSetMods performs on every enumeration, and DoFalseDiscoveryRateAnalysis sorts the
    /// matches by PEP once for its peptide- and PSM-level passes instead of twice. The FDR results must be exactly
    /// what the previous code produced. These tests run the current code and a verbatim copy of the previous code on
    /// identical, independently built sets of matches - notch-ambiguous ones, decoys, tied scores and PEPs - and
    /// compare every value both write, bit for bit.
    /// </summary>
    [TestFixture]
    internal static class FdrQValueEquivalenceTests
    {
        private const int PsmCountThresholdForInvertedQvalue = 1000;

        [Test]
        [TestCase(2500, false, TestName = "InvertedQValuesWithoutPep")]
        [TestCase(2500, true, TestName = "InvertedQValuesWithPep")]
        [TestCase(400, false, TestName = "TraditionalQValuesWithoutPep")]
        [TestCase(400, true, TestName = "TraditionalQValuesWithPep")]
        public static void FdrValuesAreUnchanged(int psmCount, bool withPep)
        {
            var commonParameters = new CommonParameters();
            var fileSpecificParameters = new List<(string fileName, CommonParameters fileSpecificParameters)>();

            List<SpectralMatch> previous = BuildMatches(psmCount, withPep, seed: 7);
            List<SpectralMatch> current = BuildMatches(psmCount, withPep, seed: 7);
            Assert.That(previous.Count(p => p.Notch == null), Is.GreaterThan(psmCount / 20), "notch-ambiguous matches must be exercised");

            PreviousDoFalseDiscoveryRateAnalysis(previous, false, fileSpecificParameters, null, null, commonParameters);
            FdrAnalysisEngine.DoFalseDiscoveryRateAnalysis(current, false, fileSpecificParameters, null, null, commonParameters);

            for (int i = 0; i < previous.Count; i++)
            {
                AssertSame(previous[i].PsmFdrInfo, current[i].PsmFdrInfo, $"PSM {i}");
                AssertSame(previous[i].PeptideFdrInfo, current[i].PeptideFdrInfo, $"peptide {i}");
                Assert.That(current[i].PsmCount, Is.EqualTo(previous[i].PsmCount));

                var expectedHypotheses = previous[i].BestMatchingBioPolymersWithSetMods.ToList();
                var actualHypotheses = current[i].BestMatchingBioPolymersWithSetMods.ToList();
                Assert.That(actualHypotheses.Count, Is.EqualTo(expectedHypotheses.Count));
                for (int h = 0; h < expectedHypotheses.Count; h++)
                {
                    var e = expectedHypotheses[h];
                    var a = actualHypotheses[h];
                    Assert.That(a.QValueNotch, Is.EqualTo(e.QValueNotch));
                    Assert.That(a.CumulativeTargetNotch, Is.EqualTo(e.CumulativeTargetNotch));
                    Assert.That(a.CumulativeDecoyNotch, Is.EqualTo(e.CumulativeDecoyNotch));
                    Assert.That(a.PeptideQValueNotch, Is.EqualTo(e.PeptideQValueNotch));
                    Assert.That(a.PeptideCumulativeTargetNotch, Is.EqualTo(e.PeptideCumulativeTargetNotch));
                    Assert.That(a.PeptideCumulativeDecoyNotch, Is.EqualTo(e.PeptideCumulativeDecoyNotch));
                }
            }
        }

        private static void AssertSame(FdrInfo expected, FdrInfo actual, string what)
        {
            Assert.That(Bits(actual.CumulativeTarget), Is.EqualTo(Bits(expected.CumulativeTarget)), what);
            Assert.That(Bits(actual.CumulativeDecoy), Is.EqualTo(Bits(expected.CumulativeDecoy)), what);
            Assert.That(Bits(actual.CumulativeTargetNotch), Is.EqualTo(Bits(expected.CumulativeTargetNotch)), what);
            Assert.That(Bits(actual.CumulativeDecoyNotch), Is.EqualTo(Bits(expected.CumulativeDecoyNotch)), what);
            Assert.That(Bits(actual.QValue), Is.EqualTo(Bits(expected.QValue)), what);
            Assert.That(Bits(actual.QValueNotch), Is.EqualTo(Bits(expected.QValueNotch)), what);
            Assert.That(Bits(actual.PEP), Is.EqualTo(Bits(expected.PEP)), what);
            Assert.That(Bits(actual.PEP_QValue), Is.EqualTo(Bits(expected.PEP_QValue)), what);
        }

        private static long Bits(double value) => BitConverter.DoubleToInt64Bits(value);

        /// <summary>
        /// A deterministic set of matches: the same seed builds the same matches, in the same order, every time
        /// </summary>
        private static List<SpectralMatch> BuildMatches(int count, bool withPep, int seed)
        {
            var random = new Random(seed);
            const string residues = "ACDEFGHILMNPQSTVWY";
            var peptides = new List<PeptideWithSetModifications>();
            for (int p = 0; p < 160; p++)
            {
                var sequence = string.Concat(Enumerable.Range(0, 6).Select(_ =>
                    string.Concat(Enumerable.Range(0, random.Next(6, 14)).Select(__ => residues[random.Next(residues.Length)])) + "K"));
                var protein = new Protein(sequence, "P" + p, isDecoy: p % 3 == 0);
                peptides.AddRange(protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()));
            }

            var commonParameters = new CommonParameters();
            var matches = new List<SpectralMatch>();
            for (int i = 0; i < count; i++)
            {
                var peptide = peptides[random.Next(peptides.Count)];
                int charge = random.Next(1, 5);
                double precursorMz = peptide.MonoisotopicMass.ToMz(charge) * (1 + (random.NextDouble() - 0.5) * 2e-5);
                var spectrum = new MzSpectrum(new double[] { 100, 200 }, new double[] { 1, 1 }, false);
                var scan = new MsDataScan(spectrum, random.Next(1, count / 2), 2, true, Polarity.Positive, random.NextDouble() * 100, null, "", MZAnalyzerType.Orbitrap, 2, null, null, null);
                var ms2 = new Ms2ScanWithSpecificMass(scan, precursorMz, charge, random.Next(2) == 0 ? "fileA" : "fileB", commonParameters);
                double score = random.Next(3, 40) + Math.Round(random.NextDouble(), 1);
                var psm = new PeptideSpectralMatch(peptide, random.Next(3), score, i, ms2, commonParameters, Ions(random));

                // a worse hypothesis sets the runner-up score; tied ones make the match ambiguous, often across notches
                if (random.Next(3) == 0)
                    psm.AddOrReplace(peptides[random.Next(peptides.Count)], score - random.Next(1, 4), random.Next(3), true, Ions(random));
                int ties = random.Next(5) == 0 ? random.Next(1, 4) : 0;
                for (int t = 0; t < ties; t++)
                    psm.AddOrReplace(peptides[random.Next(peptides.Count)], score, random.Next(3), true, Ions(random));

                psm.ResolveAllAmbiguities();
                psm.PsmFdrInfo = new FdrInfo();
                psm.PeptideFdrInfo = new FdrInfo();
                if (withPep)
                    psm.PsmFdrInfo.PEP = Math.Round(random.NextDouble(), 2) + 1e-6; // ties included, all above zero
                matches.Add(psm);
            }
            return matches;
        }

        private static List<MatchedFragmentIon> Ions(Random random)
        {
            return Enumerable.Range(0, random.Next(1, 6))
                .Select(n => new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 100 + n, n + 1, n + 1, 0), 101 + n, 10, 1))
                .ToList();
        }

        #region Verbatim copy of the previous FDR methods

        private static void PreviousDoFalseDiscoveryRateAnalysis(List<SpectralMatch> allPsms, bool doPep, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, string outputFolder, FdrAnalysisResults myAnalysisResults, CommonParameters commonParameters, bool iterativePepTraining = false)
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }

            // calculate FDR on a per-protease basis (targets and decoys for a specific protease)
            var psmsGroupedByProtease = allPsms.GroupBy(p => p.DigestionParams.DigestionAgent);

            foreach (var proteasePsms in psmsGroupedByProtease)
            {
                var psms = proteasePsms.OrderByDescending(p => p).ToList();
                var peptides = FilteredPsms.Filter(psms,
                    commonParameters,
                    includeDecoys: true,
                    includeContaminants: true,
                    includeAmbiguous: true,
                    includeHighQValuePsms: true,
                    filterAtPeptideLevel: true).FilteredPsmsList;

                if ((psms.Count > PsmCountThresholdForInvertedQvalue || FdrAnalysisEngine.QvalueThresholdOverride) & doPep)
                {
                    PreviousCalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: false);
                    if (peptides.Count > PsmCountThresholdForInvertedQvalue || FdrAnalysisEngine.QvalueThresholdOverride)
                    {
                        PreviousCalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: false);

                        //PEP will model will be developed using peptides and then applied to all PSMs. 
                        FdrAnalysisEngine.Compute_PEPValue(myAnalysisResults, psms, fileSpecificParameters, outputFolder, iterativePepTraining);

                        // peptides are ordered by PEP score from good to bad in order to select the best PSM for each peptide
                        peptides = psms
                            .OrderBy(p => p.FdrInfo.PEP)
                            .ThenByDescending(p => p)
                            .GroupBy(p => p.FullSequence)
                            .Select(p => p.FirstOrDefault()) // Get the best psm for each peptide based on PEP (default comparer is used to break ties)
                            .ToList();
                        PreviousCalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: true);

                        psms = psms.OrderBy(p => p.FdrInfo.PEP).ThenByDescending(p => p).ToList();
                        PreviousCalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: true);
                        
                    }
                    else //we have more than 100 psms but less than 100 peptides so
                    {
                        //this will be done using PSMs because we dont' have enough peptides
                        FdrAnalysisEngine.Compute_PEPValue(myAnalysisResults, psms, fileSpecificParameters, outputFolder, iterativePepTraining);
                        psms = psms.OrderBy(p => p.FdrInfo.PEP).ThenByDescending(p => p).ToList();
                        PreviousCalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: true);
                    }
                }
                else if(psms.Any(psm => psm.FdrInfo.PEP > 0)) 
                {
                    // If PEP's have been calculated, but doPEP = false, then we don't want to train another model,
                    // but we do want to calculate pep q-values
                    // really, in this case, we only need to run one or the other (i.e., only peptides or psms are passed in)
                    // but there's no mechanism to pass that info to the FDR analysis engine, so we'll do this for now
                    peptides = psms
                            .OrderBy(p => p.FdrInfo.PEP)
                            .ThenByDescending(p => p)
                            .GroupBy(p => p.FullSequence)
                            .Select(p => p.FirstOrDefault()) // Get the best psm for each peptide based on PEP (default comparer is used to break ties)
                            .ToList();
                    PreviousCalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: true);

                    psms = psms
                        .OrderBy(p => p.FdrInfo.PEP)
                        .ThenByDescending(p => p)
                        .ToList();
                    PreviousCalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: true);
                }

                //we do this section last so that target and decoy counts written in the psmtsv files are appropriate for the sort order which is by MM score
                peptides = psms
                    .OrderBy(p => p.FdrInfo.PEP)
                    .ThenByDescending(p => p)
                    .GroupBy(b => b.FullSequence)
                    .Select(b => b.FirstOrDefault()) // Get the best psm for each peptide based on PEP (default comparer is used to break ties)
                    .OrderByDescending(p => p) // Sort using the default comparer before calculating Q Values
                    .ToList();
                PreviousCalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: false);

                psms = psms.OrderByDescending(p => p).ToList();
                PreviousCalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: false);
                
                PreviousCountPsm(psms);
            }
        }

        private static void PreviousCalculateQValue(List<SpectralMatch> psms, bool peptideLevelCalculation, bool pepCalculation = false)
        {
            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;

            //set up arrays for local FDRs
            //set up arrays for local FDRs
            var allNotches = psms
                .SelectMany(p => p.BestMatchingBioPolymersWithSetMods.Select(h => h.Notch))
                .ToHashSet();

            var cumulativeTargetPerNotch = allNotches.ToDictionary(n => n, _ => 0.0);
            var cumulativeDecoyPerNotch = allNotches.ToDictionary(n => n, _ => 0.0);

            //Assign FDR values to PSMs
            List<IGrouping<int, SpectralMatchHypothesis>> ambigResultsGroupedByNotch = new(8);
            foreach (var psm in psms)
            {
                if (GlobalVariables.StopLoops) { break; }

                // Group all hypotheses by notch
                ambigResultsGroupedByNotch.Clear();
                ambigResultsGroupedByNotch.AddRange(psm.BestMatchingBioPolymersWithSetMods
                    .GroupBy(h => h.Notch));

                // Aggregate FDR Counts
                double totalHits = psm.BestMatchingBioPolymersWithSetMods.Count();
                double totalTargetHits = 0;
                double totalDecoyHits = 0;

                foreach (var group in ambigResultsGroupedByNotch)
                {
                    int notch = group.Key;
                    double groupTargetHits = group.Count(h => !h.IsDecoy);
                    double groupDecoyHits = group.Count(h => h.IsDecoy);

                    totalTargetHits += groupTargetHits;
                    totalDecoyHits += groupDecoyHits;

                    cumulativeTargetPerNotch[notch] += groupTargetHits / totalHits;
                    cumulativeDecoyPerNotch[notch] += groupDecoyHits / totalHits;
                }

                cumulativeTarget += totalTargetHits / totalHits;
                cumulativeDecoy += totalDecoyHits / totalHits;

                // Set FDR Counts
                psm.GetFdrInfo(peptideLevelCalculation).CumulativeDecoy = cumulativeDecoy;
                psm.GetFdrInfo(peptideLevelCalculation).CumulativeTarget = cumulativeTarget;

                if (ambigResultsGroupedByNotch.Count == 1) // Notch Unambiguous Psm
                {
                    int psmNotch = ambigResultsGroupedByNotch.First().Key;
                    psm.GetFdrInfo(peptideLevelCalculation).CumulativeDecoyNotch = cumulativeDecoyPerNotch[psmNotch];
                    psm.GetFdrInfo(peptideLevelCalculation).CumulativeTargetNotch = cumulativeTargetPerNotch[psmNotch];
                }
                else // Notch Ambiguous Psm 
                {
                    foreach (var hyp in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        int psmNotch = hyp.Notch;
                        if (peptideLevelCalculation)
                        {
                            hyp.PeptideCumulativeTargetNotch = cumulativeTargetPerNotch[psmNotch];
                            hyp.PeptideCumulativeDecoyNotch = cumulativeDecoyPerNotch[psmNotch];
                        }
                        else
                        {
                            hyp.CumulativeTargetNotch = cumulativeTargetPerNotch[psmNotch];
                            hyp.CumulativeDecoyNotch = cumulativeDecoyPerNotch[psmNotch];
                        }
                    }
                }
            }

            // Set FDR values
            if (pepCalculation)
            {
                PreviousPepQValueInverted(psms, peptideLevelAnalysis: peptideLevelCalculation);
            }
            else
            {
                //the QValueThreshodOverride condition here can be problematic in unit tests. 
                if (psms.Count < PsmCountThresholdForInvertedQvalue && !FdrAnalysisEngine.QvalueThresholdOverride)
                {
                   PreviousQValueTraditional(psms, peptideLevelAnalysis: peptideLevelCalculation);
                }
                else
                {
                    PreviousQValueInverted(psms, peptideLevelAnalysis: peptideLevelCalculation);
                }
            }
        }

        private static void PreviousQValueTraditional(List<SpectralMatch> psms, bool peptideLevelAnalysis)
        {
            double qValue = 0;
            Dictionary<int, double> qValueNotch = psms
                .SelectMany(p => p.BestMatchingBioPolymersWithSetMods.Select(h => h.Notch))
                .Distinct()
                .ToDictionary(p => p, _ => 0.0);

            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }
                qValue = Math.Max(qValue, psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget, 1));
                psms[i].GetFdrInfo(peptideLevelAnalysis).QValue = Math.Min(qValue, 1);

                // notch ambiguous
                if (psms[i].Notch is null)
                {
                    foreach (var notchGroup in psms[i].BestMatchingBioPolymersWithSetMods.GroupBy(p => p.Notch))
                    {
                        int notch = notchGroup.Key;
                        var first = notchGroup.First();

                        double cumDecoyNotch = peptideLevelAnalysis
                            ? first.PeptideCumulativeDecoyNotch!.Value
                            : first.CumulativeDecoyNotch!.Value;
                        double cumTargetNotch = peptideLevelAnalysis
                            ? first.PeptideCumulativeTargetNotch!.Value
                            : first.CumulativeTargetNotch!.Value;

                        qValueNotch[notch] = Math.Max(qValueNotch[notch], cumDecoyNotch / Math.Max(cumTargetNotch, 1));
                        foreach (var hyp in notchGroup)
                        {
                            if (peptideLevelAnalysis)
                            {
                                hyp.PeptideQValueNotch = Math.Min(qValueNotch[notch], 1);
                            }
                            else
                            {
                                hyp.QValueNotch = Math.Min(qValueNotch[notch], 1);
                            }
                        }
                    }
                }
                else // notch unambiguous - original implementation. 
                {
                    int notch = psms[i].Notch.Value;
                    qValueNotch[notch] = Math.Max(qValueNotch[notch], psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoyNotch / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTargetNotch, 1));
                    psms[i].GetFdrInfo(peptideLevelAnalysis).QValueNotch = Math.Min(qValueNotch[notch], 1);
                }
            }
        }

        private static void PreviousQValueInverted(List<SpectralMatch> psms, bool peptideLevelAnalysis)
        {
            Dictionary<int, double> qValueNotch = psms
                .SelectMany(p => p.BestMatchingBioPolymersWithSetMods.Select(h => h.Notch))
                .Distinct()
                .ToDictionary(p => p, _ => 2.0);
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Min(qValue, (psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget, 1));
                psms[i].GetFdrInfo(peptideLevelAnalysis).QValue = Math.Min(qValue, 1);

                // Notch ambiguous
                if (psms[i].Notch is null)
                {
                    foreach (var notchGroup in psms[i].BestMatchingBioPolymersWithSetMods.GroupBy(p => p.Notch))
                    {
                        int notch = notchGroup.Key;
                        var first = notchGroup.First();

                        double cumDecoyNotch = peptideLevelAnalysis
                            ? first.PeptideCumulativeDecoyNotch!.Value
                            : first.CumulativeDecoyNotch!.Value;
                        double cumTargetNotch = peptideLevelAnalysis
                            ? first.PeptideCumulativeTargetNotch!.Value
                            : first.CumulativeTargetNotch!.Value;

                        qValueNotch[notch] = Math.Min(qValueNotch[notch], (cumDecoyNotch + 1) / Math.Max(cumTargetNotch, 1));
                        foreach (var hyp in notchGroup)
                        {
                            if (peptideLevelAnalysis)
                            {
                                hyp.PeptideQValueNotch = Math.Min(qValueNotch[notch], 1);
                            }
                            else
                            {
                                hyp.QValueNotch = Math.Min(qValueNotch[notch], 1);
                            }
                        }
                    }
                }
                else
                {
                    int notch = psms[i].Notch.Value;
                    qValueNotch[notch] = Math.Min(qValueNotch[notch], (psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoyNotch + 1) / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTargetNotch, 1));
                    psms[i].GetFdrInfo(peptideLevelAnalysis).QValueNotch = Math.Min(qValueNotch[notch], 1);
                }
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        private static void PreviousPepQValueInverted(List<SpectralMatch> psms, bool peptideLevelAnalysis)
        {
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Min(qValue, (psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget);

                psms[i].GetFdrInfo(peptideLevelAnalysis).PEP_QValue = qValue;
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        private static void PreviousCountPsm(List<SpectralMatch> proteasePsms)
        {
            // exclude ambiguous psms and has a fdr cutoff = 0.01
            var allUnambiguousPsms = proteasePsms.Where(psm => psm.FullSequence != null).ToList();

            var unambiguousPsmsLessThanOnePercentFdr = allUnambiguousPsms.Where(psm =>
                    psm.FdrInfo.QValue <= 0.01
                    && psm.FdrInfo.QValueNotch <= 0.01)
                .GroupBy(p => p.FullSequence);

            Dictionary<string, int> sequenceToPsmCount = new Dictionary<string, int>();

            foreach (var sequenceGroup in unambiguousPsmsLessThanOnePercentFdr)
            {
                sequenceToPsmCount.TryAdd(sequenceGroup.First().FullSequence, sequenceGroup.Count());
            }

            foreach (SpectralMatch psm in allUnambiguousPsms)
            {
                if (sequenceToPsmCount.TryGetValue(psm.FullSequence, out int count))
                {
                    psm.PsmCount = count;
                }
            }
        }

        #endregion
    }
}
