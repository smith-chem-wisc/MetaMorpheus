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
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// PEP is trained on one match per peptide: each full sequence's best-scoring match
    /// (<see cref="SpectralMatchGroup.GetBestMatches"/>). The model then scores every match. Every peptide-level step
    /// after PEP stands each peptide on that same match: the classic q-value ranks the peptides by its score, the PEP
    /// q-value by its PEP, and the result files and peptide-level filters collapse each peptide to it. No step picks a
    /// peptide's match again by PEP, so a peptide's values never depend on PEP noise among its other matches.
    /// </summary>
    [TestFixture]
    internal static class PeptideStandInTests
    {
        private static readonly CommonParameters Common = new();
        private const string Residues = "ACDEFGHILNQSTVWY";

        /// <summary>
        /// A tryptic peptide of its own protein, one per index: the same index always gives the same sequence.
        /// </summary>
        private static PeptideWithSetModifications Peptide(int index, bool decoy)
        {
            var letters = new char[7];
            int n = index;
            for (int i = 0; i < letters.Length; i++)
            {
                letters[i] = Residues[n % Residues.Length];
                n /= Residues.Length;
            }
            string sequence = "G" + new string(letters) + "K";
            var protein = new Protein(sequence, (decoy ? "DECOY_" : "") + "P" + index, isDecoy: decoy);
            return protein.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .Single(p => p.BaseSequence == sequence);
        }

        private static SpectralMatch Psm(PeptideWithSetModifications peptide, double score, double pep, int scanNumber)
        {
            var spectrum = new MzSpectrum(new double[] { 100, 200 }, new double[] { 1, 1 }, false);
            var scan = new MsDataScan(spectrum, scanNumber, 2, true, Polarity.Positive, scanNumber / 10.0, null, "", MZAnalyzerType.Orbitrap, 2, null, null, null);
            var ms2 = new Ms2ScanWithSpecificMass(scan, peptide.MonoisotopicMass.ToMz(2), 2, "file", Common);
            var ions = new List<MatchedFragmentIon>
            {
                new(new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0), 101, 10, 1)
            };
            var psm = new PeptideSpectralMatch(peptide, 0, score, scanNumber, ms2, Common, ions);
            psm.ResolveAllAmbiguities();
            psm.PsmFdrInfo = new FdrInfo { PEP = pep };
            psm.PeptideFdrInfo = new FdrInfo();
            return psm;
        }

        /// <summary>
        /// Ten target peptides scored 40 to 31 with PEP 0.01, and five decoys scored 10 to 6, two with PEP 0.1 and three
        /// with PEP 0.9. Peptide X (index 100) has a best-scoring match (score 50, PEP 0.3) and a weak match with a lower
        /// PEP (score 1, PEP 0.001). By score, X's best match is above every decoy and its weak match below all five. By
        /// PEP, its weak match is above every decoy and its best match below two of them.
        /// </summary>
        private static (List<SpectralMatch> All, SpectralMatch XBest, SpectralMatch XLowestPep) PeptideXFixture(bool withXBest = true, bool withXLowestPep = true)
        {
            var all = new List<SpectralMatch>();
            int scan = 1;
            for (int t = 0; t < 10; t++)
                all.Add(Psm(Peptide(t, decoy: false), 40 - t, 0.01, scan++));
            for (int d = 0; d < 5; d++)
                all.Add(Psm(Peptide(50 + d, decoy: true), 10 - d, d < 2 ? 0.1 : 0.9, scan++));

            var x = Peptide(100, decoy: false);
            SpectralMatch xBest = withXBest ? Psm(x, 50, 0.3, scan++) : null;
            SpectralMatch xLowestPep = withXLowestPep ? Psm(x, 1, 0.001, scan++) : null;
            if (xBest != null) all.Add(xBest);
            if (xLowestPep != null) all.Add(xLowestPep);
            return (all, xBest, xLowestPep);
        }

        /// <summary>
        /// The FDR passes over PEPs already computed (doPep false): every peptide-level pass that follows PEP training.
        /// </summary>
        private static void RunFdr(List<SpectralMatch> psms) =>
            FdrAnalysisEngine.DoFalseDiscoveryRateAnalysis(psms, false, new List<(string, CommonParameters)>(), null, null, Common);

        /// <summary>
        /// The peptides as the result files and every peptide-level filter see them: one match per peptide.
        /// </summary>
        private static List<SpectralMatch> PeptideRows(List<SpectralMatch> psms, bool includeHighQValuePsms = true) =>
            FilteredPsms.Filter(psms, Common, includeDecoys: true, includeContaminants: true, includeAmbiguous: true,
                includeHighQValuePsms: includeHighQValuePsms, filterAtPeptideLevel: true).FilteredPsmsList;

        [Test]
        public static void APeptideStandsOnItsBestScoringMatch_ForBothQValuesAndItsRow()
        {
            var (all, xBest, xLowestPep) = PeptideXFixture();
            RunFdr(all);

            // What X gets when its best-scoring match is its only match, and when its weak match is.
            var (onlyBest, xBestAlone, _) = PeptideXFixture(withXLowestPep: false);
            RunFdr(onlyBest);
            var (onlyWeak, _, xWeakAlone) = PeptideXFixture(withXBest: false);
            RunFdr(onlyWeak);
            Assume.That(xBestAlone.PeptideFdrInfo.QValue, Is.Not.EqualTo(xWeakAlone.PeptideFdrInfo.QValue),
                "the fixture must tell the two matches apart by q-value");
            Assume.That(xBestAlone.PeptideFdrInfo.PEP_QValue, Is.Not.EqualTo(xWeakAlone.PeptideFdrInfo.PEP_QValue),
                "the fixture must tell the two matches apart by PEP q-value");

            var row = PeptideRows(all).Single(p => p.FullSequence == xBest.FullSequence);
            Assert.That(row, Is.SameAs(xBest), "the peptide's row is its best-scoring match");
            Assert.That(row.PeptideFdrInfo.QValue, Is.EqualTo(xBestAlone.PeptideFdrInfo.QValue));
            Assert.That(row.PeptideFdrInfo.QValueNotch, Is.EqualTo(xBestAlone.PeptideFdrInfo.QValueNotch));
            Assert.That(row.PeptideFdrInfo.CumulativeTarget, Is.EqualTo(xBestAlone.PeptideFdrInfo.CumulativeTarget));
            Assert.That(row.PeptideFdrInfo.CumulativeDecoy, Is.EqualTo(xBestAlone.PeptideFdrInfo.CumulativeDecoy));
            Assert.That(row.PeptideFdrInfo.PEP_QValue, Is.EqualTo(xBestAlone.PeptideFdrInfo.PEP_QValue),
                "the peptide's PEP q-value ranks it by its best-scoring match's PEP");

            // The weak match carries no peptide-level values of its own: no step stood the peptide on it.
            Assert.That(xLowestPep.PeptideFdrInfo.QValue, Is.EqualTo(2));
            Assert.That(xLowestPep.PeptideFdrInfo.PEP_QValue, Is.EqualTo(2));

            // X passes the default peptide-level filter at 1% on its best score.
            Assert.That(PeptideRows(all, includeHighQValuePsms: false), Does.Contain(xBest));
        }

        [Test]
        public static void NearZeroPepNoiseMovesNoPeptideLevelValue()
        {
            // Peptide Y has two matches, scored 30 and 20, with two decoys scored between them. Their PEPs are 0 and one
            // float32 step above it (5.96e-8), swapped between the two runs, as two PEP trainings differing only in
            // noise would give them.
            List<SpectralMatch> Build(double pepOfY30, double pepOfY20)
            {
                var all = new List<SpectralMatch>();
                int scan = 1;
                for (int t = 0; t < 10; t++)
                    all.Add(Psm(Peptide(t, decoy: false), 40 - t, 0.5, scan++));
                all.Add(Psm(Peptide(50, decoy: true), 25, 0.9, scan++));
                all.Add(Psm(Peptide(51, decoy: true), 24, 0.9, scan++));
                for (int t = 10; t < 20; t++)
                    all.Add(Psm(Peptide(t, decoy: false), 19 - (t - 10) * 0.5, 0.5, scan++));
                var y = Peptide(200, decoy: false);
                all.Add(Psm(y, 30, pepOfY30, scan++));
                all.Add(Psm(y, 20, pepOfY20, scan++));
                return all;
            }

            var a = Build(pepOfY30: 5.96e-8, pepOfY20: 0);
            var b = Build(pepOfY30: 0, pepOfY20: 5.96e-8);
            RunFdr(a);
            RunFdr(b);

            var rowsA = PeptideRows(a).ToDictionary(p => p.FullSequence);
            var rowsB = PeptideRows(b).ToDictionary(p => p.FullSequence);
            Assert.That(rowsB.Keys, Is.EquivalentTo(rowsA.Keys));
            foreach (var (sequence, rowA) in rowsA)
            {
                var rowB = rowsB[sequence];
                Assert.That(rowB.ScanNumber, Is.EqualTo(rowA.ScanNumber), sequence + ": the same match stands for the peptide");
                Assert.That(rowB.PeptideFdrInfo.QValue, Is.EqualTo(rowA.PeptideFdrInfo.QValue), sequence);
                Assert.That(rowB.PeptideFdrInfo.QValueNotch, Is.EqualTo(rowA.PeptideFdrInfo.QValueNotch), sequence);
                Assert.That(rowB.PeptideFdrInfo.PEP_QValue, Is.EqualTo(rowA.PeptideFdrInfo.PEP_QValue), sequence);
            }

            int TargetsAtOnePercent(List<SpectralMatch> psms) => PeptideRows(psms, includeHighQValuePsms: false).Count(p => !p.IsDecoy);
            Assert.That(TargetsAtOnePercent(b), Is.EqualTo(TargetsAtOnePercent(a)));
        }

        [Test]
        public static void EveryPeptideStandsOnTheMatchPepTrainedOn()
        {
            // Twelve target and four decoy peptides, each with three matches whose PEPs run opposite to their scores, so a
            // peptide's lowest-PEP match is never its best-scoring one.
            var all = new List<SpectralMatch>();
            int scan = 1;
            for (int i = 0; i < 16; i++)
            {
                bool decoy = i >= 12;
                var peptide = Peptide(300 + i, decoy);
                double top = decoy ? 12 - i * 0.1 : 40 - i;
                for (int k = 0; k < 3; k++)
                    all.Add(Psm(peptide, top - k, (decoy ? 0.5 : 0.01) + (2 - k) * 0.001, scan++));
            }
            RunFdr(all);

            var trainedOn = SpectralMatchGroup.GroupByBaseSequence(all)
                .SelectMany(group => group.GetBestMatches())
                .ToDictionary(p => p.FullSequence);
            var rows = PeptideRows(all);
            Assert.That(rows.Select(p => p.FullSequence), Is.EquivalentTo(trainedOn.Keys));
            foreach (var row in rows)
            {
                Assert.That(row, Is.SameAs(trainedOn[row.FullSequence]), row.FullSequence);
            }

            // Only those matches carry peptide-level values.
            foreach (var psm in all.Where(p => !ReferenceEquals(p, trainedOn[p.FullSequence])))
            {
                Assert.That(psm.PeptideFdrInfo.QValue, Is.EqualTo(2), psm.FullSequence);
                Assert.That(psm.PeptideFdrInfo.PEP_QValue, Is.EqualTo(2), psm.FullSequence);
            }
        }

        [Test]
        public static void PepTrainingCountsAcceptedPeptidesOnTheMatchesItTrainedOn()
        {
            // Ranked by PEP on each peptide's best-scoring match, X (PEP 0.3) sits below two decoys and has a q-value of
            // 3/11; the ten other targets are at 0.1. On X's lowest-PEP match (0.001) it would be first, and all eleven
            // targets would be at 1/11.
            var (all, _, _) = PeptideXFixture();
            Assert.That(PepAnalysisEngine.CountAcceptedPeptides(all, 0.2), Is.EqualTo(10));
        }
    }
}
