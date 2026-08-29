using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Chemistry;
using EngineLayer;
using EngineLayer.NonSpecificEnzymeSearch;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    /// <summary>
    /// Covers <see cref="NonSpecificEnzymeSearchEngine.ResolveFdrCategorySpecificPsms"/>, which decides
    /// which spectral match survives when a non-specific search has produced one per FDR category for
    /// the same spectrum.
    ///
    /// It had no tests. Mutation testing (#2777) found the consequence: the comparisons that pick the
    /// winner can be inverted, and the score cutoff can be shifted, without the suite noticing --
    /// 36 surviving mutants in this file, of which the cluster in this method is the largest. A method
    /// that chooses between candidate identifications and cannot tell the better from the worse is the
    /// kind of defect that produces a plausible answer rather than an error.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class NonSpecificEnzymeFdrCategoryTests
    {
        private static readonly CommonParameters CommonParams = new();

        /// <summary>
        /// One PSM per category for a single spectrum, with the scores given. A null score leaves that
        /// category empty for the spectrum, which is the shape the engine sees when a category matched
        /// nothing.
        /// </summary>
        private static List<SpectralMatch>[] OneSpectrumAcrossCategories(params double?[] scoresByCategory)
        {
            var protein = new Protein("MNKNNKNNNKNNNNKPEPTIDEKPEPTIDER", "target");
            var peptides = protein
                .Digest(CommonParams.DigestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>().ToList();
            Assert.That(peptides.Count, Is.GreaterThanOrEqualTo(scoresByCategory.Length),
                "premise: the fixture protein yields a distinct peptide per category");

            var dataFile = new TestDataFile(peptides.Cast<IBioPolymerWithSetMods>().ToList());

            var allPsms = new List<SpectralMatch>[scoresByCategory.Length];
            for (int category = 0; category < scoresByCategory.Length; category++)
            {
                if (scoresByCategory[category] == null)
                {
                    allPsms[category] = new List<SpectralMatch> { null };
                    continue;
                }

                var peptide = peptides[category];
                var scan = new Ms2ScanWithSpecificMass(
                    dataFile.GetOneBasedScan(2), peptide.MonoisotopicMass.ToMz(1), 1, null, CommonParams);

                // Scan number 1 for every category: these are competing identifications of ONE spectrum,
                // which is the situation the method exists to resolve.
                var psm = new PeptideSpectralMatch(peptide, 0, scoresByCategory[category]!.Value, 1, scan,
                    CommonParams, new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                allPsms[category] = new List<SpectralMatch> { psm };
            }
            return allPsms;
        }

        private static List<SpectralMatch> Resolve(List<SpectralMatch>[] allPsms) =>
            NonSpecificEnzymeSearchEngine.ResolveFdrCategorySpecificPsms(
                allPsms, 1, "TestTask", CommonParams,
                new List<(string, CommonParameters)>());

        /// <summary>
        /// With the q-values tied, the highest-scoring identification wins.
        ///
        /// The scores are deliberately 5, 9, 3 rather than ascending or descending, because the two ways
        /// this can be got wrong both look right against a sorted fixture. Taking the first candidate
        /// unconditionally yields 5; taking the last yields 3; only actually comparing yields 9.
        ///
        /// That shape is what kills the mutants #2777 lists at line 553:
        ///   `currentQValue == lowestQ && currentPsm.Score > bestPsm.Score`
        ///     -> `||`  makes every later candidate win, so the answer becomes 3
        ///     -> `!=`  makes none of them win, so the answer stays 5
        ///     -> `<`   inverts the comparison, so the answer becomes 3
        /// </summary>
        [Test]
        public static void HighestScoringMatchWinsWhenQValuesAreTied()
        {
            var allPsms = OneSpectrumAcrossCategories(5, 9, 3);

            var best = Resolve(allPsms);

            Assert.That(best, Has.Count.EqualTo(1), "one spectrum in, one identification out");
            Assert.That(best[0].Score, Is.EqualTo(9).Within(1e-9),
                "the best-scoring candidate must win a q-value tie");
        }

        /// <summary>
        /// The losers are not merely unreported -- they are removed from the categories, because the FDR
        /// recalculation at the end of the method runs over what is left. A candidate that stays behind
        /// contributes to an FDR it was not chosen for.
        /// </summary>
        [Test]
        public static void LosingCandidatesAreClearedFromTheirCategories()
        {
            var allPsms = OneSpectrumAcrossCategories(5, 9, 3);

            var best = Resolve(allPsms);

            var survivors = allPsms.SelectMany(category => category).Where(psm => psm != null).ToList();
            Assert.That(survivors, Has.Count.EqualTo(1), "only the winner should remain in the categories");
            Assert.That(survivors[0].Score, Is.EqualTo(best[0].Score).Within(1e-9));
        }

        /// <summary>
        /// A category that matched nothing for this spectrum is skipped rather than treated as a
        /// candidate with no score.
        /// </summary>
        [Test]
        public static void EmptyCategoriesAreSkipped()
        {
            var allPsms = OneSpectrumAcrossCategories(null, 7, null);

            var best = Resolve(allPsms);

            Assert.That(best, Has.Count.EqualTo(1));
            Assert.That(best[0].Score, Is.EqualTo(7).Within(1e-9));
        }

        /// <summary>
        /// The returned list is ordered by q-value ascending, then by score descending -- the contract
        /// the caller relies on to take the best identifications first. Asserted over several spectra,
        /// since a single-element list is ordered by accident.
        /// </summary>
        [Test]
        public static void ResultIsOrderedByQValueThenDescendingScore()
        {
            var protein = new Protein("MNKNNKNNNKNNNNKPEPTIDEKPEPTIDER", "target");
            var peptides = protein
                .Digest(CommonParams.DigestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>().ToList();
            var dataFile = new TestDataFile(peptides.Cast<IBioPolymerWithSetMods>().ToList());

            // One category, several spectra, scores out of order so the ordering has work to do.
            var psms = new List<SpectralMatch>();
            double[] scores = { 4, 11, 7, 2, 9 };
            for (int i = 0; i < scores.Length && i < peptides.Count; i++)
            {
                var peptide = peptides[i];
                var scan = new Ms2ScanWithSpecificMass(
                    dataFile.GetOneBasedScan(2), peptide.MonoisotopicMass.ToMz(1), 1, null, CommonParams);
                var psm = new PeptideSpectralMatch(peptide, 0, scores[i], i + 1, scan, CommonParams,
                    new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                psms.Add(psm);
            }

            var best = Resolve(new[] { psms });

            Assert.That(best, Is.Not.Empty);
            for (int i = 1; i < best.Count; i++)
            {
                double previousQ = best[i - 1].PsmFdrInfo.QValue;
                double currentQ = best[i].PsmFdrInfo.QValue;
                Assert.That(previousQ, Is.LessThanOrEqualTo(currentQ), "q-values must ascend");
                if (previousQ == currentQ)
                {
                    Assert.That(best[i - 1].Score, Is.GreaterThanOrEqualTo(best[i].Score),
                        "within one q-value, scores must descend");
                }
            }
        }
    }
}
