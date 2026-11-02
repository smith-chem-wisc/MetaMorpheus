using Chemistry;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test
{
    /// <summary>
    /// <see cref="FdrInfoSnapshot"/> must hand back every match's whole-search FDR after a per-file pass, including
    /// the per-hypothesis notch values of a notch-ambiguous match, whose hypotheses carry different q-value notches.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class FdrInfoSnapshotTests
    {
        private const int NumNotches = 2;

        [Test]
        public static void Clone_CopiesEveryValueIntoANewObject()
        {
            var original = new FdrInfo();
            var properties = typeof(FdrInfo).GetProperties();
            for (int i = 0; i < properties.Length; i++)
                properties[i].SetValue(original, 0.125 * (i + 1));

            FdrInfo clone = original.Clone();

            Assert.That(clone, Is.Not.SameAs(original));
            foreach (var property in properties)
                Assert.That(property.GetValue(clone), Is.EqualTo(property.GetValue(original)), property.Name);

            double qValue = original.QValue;
            clone.QValue = qValue + 1;
            Assert.That(original.QValue, Is.EqualTo(qValue), "writing to the clone leaves the original unchanged");
        }

        /// <summary>
        /// Two files. File A holds a notch-ambiguous match X (one hypothesis at notch 0, one at notch 1). Over both
        /// files X's two hypotheses get different q-value notches, because file B's decoy at notch 1 counts only for
        /// notch 1. Run per file, as the individual-file tables do, X's notch-1 hypothesis gets file A's value instead.
        /// After the snapshot is disposed every value of every match and of every hypothesis, and every object, must
        /// be the whole-search one again.
        /// </summary>
        [Test]
        public static void PerFilePass_InsideTheSnapshot_LeavesEveryHypothesisItsOwnWholeSearchNotchValues()
        {
            var commonParameters = new CommonParameters();
            var b1 = Match("B.mzML", 1, 50, 1, Peptide("AAAAAAAK"));
            var b2 = Match("B.mzML", 2, 45, 1, Peptide("CCCCCCCK"));
            var b3 = Match("B.mzML", 3, 40, 1, Peptide("DDDDDDDK", isDecoy: true));
            var a1 = Match("A.mzML", 1, 35, 0, Peptide("EEEEEEEK"));
            var x = Match("A.mzML", 2, 30, 0, Peptide("FFFFFFFK"));
            x.AddOrReplace(Peptide("GGGGGGGK"), 30, 1, true, new List<MatchedFragmentIon>());
            var a3 = Match("A.mzML", 3, 20, 0, Peptide("HHHHHHHK", isDecoy: true));
            var all = new List<SpectralMatch> { b1, b2, b3, a1, x, a3 };
            all.ForEach(m => m.ResolveAllAmbiguities());

            // Whole search, as PostSearchAnalysisTask.Run does before anything is written.
            RunFdr(all, commonParameters, "PSM");

            var xNotch0 = x.BestMatchingBioPolymersWithSetMods.Single(h => h.Notch == 0);
            var xNotch1 = x.BestMatchingBioPolymersWithSetMods.Single(h => h.Notch == 1);
            Assert.That(x.Notch, Is.Null, "X is notch-ambiguous");
            Assert.That(xNotch0.QValueNotch, Is.Not.Null);
            Assert.That(xNotch1.QValueNotch, Is.Not.Null);
            Assert.That(xNotch1.QValueNotch, Is.Not.EqualTo(xNotch0.QValueNotch),
                "the whole-search pass gives X's two hypotheses different q-value notches");
            Assert.That(xNotch1.PeptideQValueNotch, Is.Not.EqualTo(xNotch0.PeptideQValueNotch));

            var before = all.Select(State).ToList();
            var xNotch1Before = HypothesisValues(xNotch1);
            var originalFdrInfos = all.Select(m => (m.PsmFdrInfo, m.PeptideFdrInfo)).ToList();
            var originalFdrValues = originalFdrInfos.Select(f => (Values(f.PsmFdrInfo), Values(f.PeptideFdrInfo))).ToList();

            using (FdrInfoSnapshot.Take(all))
            {
                // The individual-file tables: PSM FDR per file, then peptide FDR on each file's best match per sequence.
                foreach (var file in all.GroupBy(m => m.FullFilePath))
                {
                    RunFdr(file.ToList(), commonParameters, "PSM");
                    RunFdr(file.OrderByDescending(m => m).GroupBy(m => m.FullSequence).Select(g => g.First()).ToList(), commonParameters, "peptide");
                }

                Assert.That(HypothesisValues(xNotch1), Is.Not.EqualTo(xNotch1Before),
                    "the per-file pass must change X's notch-1 hypothesis, or this test proves nothing");
                for (int i = 0; i < all.Count; i++)
                {
                    Assert.That((Values(originalFdrInfos[i].PsmFdrInfo), Values(originalFdrInfos[i].PeptideFdrInfo)), Is.EqualTo(originalFdrValues[i]),
                        "the per-file pass writes into copies, never into the whole-search FdrInfo objects");
                }
            }

            for (int i = 0; i < all.Count; i++)
            {
                Assert.That(all[i].PsmFdrInfo, Is.SameAs(originalFdrInfos[i].PsmFdrInfo));
                Assert.That(all[i].PeptideFdrInfo, Is.SameAs(originalFdrInfos[i].PeptideFdrInfo));
                Assert.That(State(all[i]), Is.EqualTo(before[i]), $"match {i} must have its whole-search FDR back, hypothesis by hypothesis");
            }
            Assert.That(x.BestMatchingBioPolymersWithSetMods, Is.EqualTo(new[] { xNotch0, xNotch1 }).Or.EqualTo(new[] { xNotch1, xNotch0 }),
                "X keeps both hypothesis objects");
            Assert.That(xNotch1.QValueNotch, Is.Not.EqualTo(xNotch0.QValueNotch), "and they still carry different q-value notches");
        }

        private static void RunFdr(List<SpectralMatch> matches, CommonParameters commonParameters, string analysisType) =>
            new FdrAnalysisEngine(matches, NumNotches, commonParameters, new List<(string, CommonParameters)> { ("A.mzML", commonParameters), ("B.mzML", commonParameters) },
                new List<string>(), analysisType, doPEP: false).Run();

        private static PeptideWithSetModifications Peptide(string sequence, bool isDecoy = false) =>
            new(sequence, [], p: new Protein(sequence, (isDecoy ? "DECOY_" : "") + sequence, isDecoy: isDecoy));

        private static SpectralMatch Match(string file, int scanNumber, double score, int notch, PeptideWithSetModifications peptide)
        {
            var scan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), scanNumber, 2, true,
                Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=" + scanNumber,
                double.NaN, null, null, double.NaN, null, DissociationType.HCD, 1, null);
            var ms2 = new Ms2ScanWithSpecificMass(scan, peptide.MonoisotopicMass.ToMz(1), 1, file, new CommonParameters());
            return new PeptideSpectralMatch(peptide, notch, score, scanNumber, ms2, new CommonParameters(), new List<MatchedFragmentIon>());
        }

        private static (double, double, double, double, double, double, double, double) Values(FdrInfo f) =>
            (f.CumulativeTarget, f.CumulativeDecoy, f.CumulativeTargetNotch, f.CumulativeDecoyNotch, f.QValue, f.QValueNotch, f.PEP, f.PEP_QValue);

        private record HypothesisState(int Notch, double? QValueNotch, double? CumulativeTargetNotch, double? CumulativeDecoyNotch,
            double? PeptideQValueNotch, double? PeptideCumulativeTargetNotch, double? PeptideCumulativeDecoyNotch, double? Pep);

        private static HypothesisState HypothesisValues(SpectralMatchHypothesis h) => new(h.Notch, h.QValueNotch, h.CumulativeTargetNotch,
            h.CumulativeDecoyNotch, h.PeptideQValueNotch, h.PeptideCumulativeTargetNotch, h.PeptideCumulativeDecoyNotch, h.PEP);

        /// <summary>Every FDR value a match and each of its hypotheses carry, as text, so a failure shows what changed.</summary>
        private static string State(SpectralMatch m) => $"PSM {Values(m.PsmFdrInfo)}; peptide {Values(m.PeptideFdrInfo)}; PsmCount {m.PsmCount}; "
            + string.Join("; ", m.BestMatchingBioPolymersWithSetMods.Select(HypothesisValues));
    }
}
