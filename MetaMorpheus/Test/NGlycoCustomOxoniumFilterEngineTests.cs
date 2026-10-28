using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The custom-oxonium gate inside GlycoSearchEngine.FindNGlycan, driven through the engine rather than through
    /// GlycoPeptides.CustomOxoniumFilter alone, so the tests fail if the gate is removed from the N-search path.
    /// </summary>
    /// <remarks>
    /// Uses the N-glycopeptide the TestNGlyco fragment test already matches: EEQYNSTYR carrying HexNAc(4)Hex(3)Fuc(1) in
    /// Glyco_3383.mgf (HCD). FindNGlycan is private, so it is called by reflection, as GlycoSearchEngineTest does with
    /// CreateGsm; the oxonium intensities are handed in, so each test sets exactly the evidence the gate reads.
    /// NonParallelizable because registering a custom monosaccharide mutates process-wide Glycan statics.
    /// </remarks>
    [TestFixture]
    [NonParallelizable]
    public class NGlycoCustomOxoniumFilterEngineTests
    {
        private const int Ion512 = 51219700; // 512.197 m/z * 1e5, not a built-in oxonium ion
        private const int OxoniumIndex_R138 = 4;
        private const int OxoniumIndex_HexNAc204 = 9;
        private const string NGlycanDatabase = "NGlycan.gdb";

        // Reset to startup state, as CustomOxoniumFilterTests does: the shipped MonosaccharidesCustom.tsv may register rows.
        private static void RestoreStartupMonosaccharides()
        {
            Glycan.ResetCustomMonosaccharides();
            string shipped = GlobalVariables.CustomMonosaccharidePath;
            if (File.Exists(shipped))
            {
                GlycanDatabase.LoadCustomMonosaccharides(shipped);
            }
        }

        /// <summary>
        /// Runs FindNGlycan for EEQYNSTYR against the Glyco_3383 scan and returns the matches it adds.
        /// </summary>
        /// <param name="oxoniumIonFilter"> The task's OxoniumIonFilt setting. </param>
        /// <param name="customIonObserved"> Whether the oxonium intensities report the registered custom ion (if any). </param>
        /// <param name="customIonIntensity"> The custom ion's intensity when reported, against a 138.055 reference of 100. </param>
        /// <param name="customOxoniumThreshold"> The engine's CustomOxoniumRelativeIntensityThreshold. </param>
        private static List<GlycoSpectralMatch> FindNGlycanForKnownGlycopeptide(bool oxoniumIonFilter, bool customIonObserved,
            double customIonIntensity = 500, double customOxoniumThreshold = GlycoPeptides.DefaultCustomOxoniumRelativeIntensityThreshold)
        {
            var commonParameters = new CommonParameters(dissociationType: DissociationType.HCD, trimMsMsPeaks: false,
                precursorMassTolerance: new PpmTolerance(10));

            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/Glyco_3383.mgf");
            var msDataFile = new MyFileManager(true).LoadFile(filePath, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();
            var scan = scans[0];

            var peptide = new Protein("TKPREEQYNSTYR", "accession")
                .Digest(new DigestionParams(minPeptideLength: 7), new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>()
                .Last(); // EEQYNSTYR

            // Built after any custom monosaccharide is registered, so the loaded N-glycans' Kind[] has its slot.
            var engine = new GlycoSearchEngine(new List<GlycoSpectralMatch>[scans.Length], scans, new List<PeptideWithSetModifications> { peptide },
                null, null, 0, commonParameters, null, "OGlycan.gdb", NGlycanDatabase, GlycoSearchType.NGlycanSearch,
                glycoSearchTopNum: 30, maxOGlycanNum: 3, oxoniumIonFilter: oxoniumIonFilter, nestedIds: null)
            {
                CustomOxoniumRelativeIntensityThreshold = customOxoniumThreshold
            };

            double[] oxoniumIonIntensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
            oxoniumIonIntensities[OxoniumIndex_R138] = 100;      // the reference the custom-ion test divides by
            oxoniumIonIntensities[OxoniumIndex_HexNAc204] = 1000;
            if (customIonObserved && Glycan.HasCustomOxoniumIons)
            {
                oxoniumIonIntensities[Glycan.AllOxoniumIons.Length] = customIonIntensity; // first custom slot, ratio customIonIntensity / 100
            }

            double possibleGlycanMassLow = commonParameters.PrecursorMassTolerance.GetMinimumValue(scan.PrecursorMass) - peptide.MonoisotopicMass;
            var possibleMatches = new List<GlycoSpectralMatch>();
            object[] args = { scan, 0, 3, peptide, 0, possibleGlycanMassLow, oxoniumIonIntensities, possibleMatches };

            MethodInfo findNGlycan = typeof(GlycoSearchEngine).GetMethod("FindNGlycan", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.That(findNGlycan, Is.Not.Null, "FindNGlycan was renamed or its signature changed.");
            findNGlycan.Invoke(engine, args);

            return (List<GlycoSpectralMatch>)args[7];
        }

        // HexNAc(4)Hex(3)Fuc(1), the glycan the fragment test matches in this scan.
        private static bool IsKnownGlycan(GlycoSpectralMatch match)
        {
            byte[] expected = GlycanDatabase.ParseComposition("HexNAc(4)Hex(3)Fuc(1)");
            return match.NGlycan != null && match.NGlycan.Single().Kind.SequenceEqual(expected);
        }

        private static void RegisterCustomMonosaccharideWithIon()
        {
            Glycan.ResetCustomMonosaccharides();
            Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
        }

        [Test]
        public void FindNGlycan_NoCustomIons_FindsTheKnownGlycopeptide()
        {
            // Baseline the other tests depend on: without custom monosaccharides and custom ions, the filter-on N-search identifies the glycopeptide.
            try
            {
                Glycan.ResetCustomMonosaccharides();

                var matches = FindNGlycanForKnownGlycopeptide(oxoniumIonFilter: true, customIonObserved: false); // N-Search should behave the same with or without the custom ion, since it is not registered.

                Assert.That(matches, Is.Not.Empty);
                Assert.That(matches.Any(IsKnownGlycan), Is.True);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void FindNGlycan_CustomIonObservedButGlycanLacksIt_FilterOn_RejectsEveryGlycan()
        {
            // The gate itself: the spectrum reports the custom ion, no glycan in NGlycan.gdb carries the custom
            // monosaccharide, so under the strict rule every candidate is rejected. Fails if the gate is removed.
            try
            {
                RegisterCustomMonosaccharideWithIon();

                var matches = FindNGlycanForKnownGlycopeptide(oxoniumIonFilter: true, customIonObserved: true);

                Assert.That(matches, Is.Empty);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void FindNGlycan_CustomIonAbsentAndGlycanLacksIt_FilterOn_KeepsTheGlycopeptide()
        {
            // A registered custom ion that the spectrum does not show must not reject glycans without its monosaccharide.
            try
            {
                RegisterCustomMonosaccharideWithIon();

                var matches = FindNGlycanForKnownGlycopeptide(oxoniumIonFilter: true, customIonObserved: false);

                Assert.That(matches.Any(IsKnownGlycan), Is.True);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void FindNGlycan_CustomIonObservedButGlycanLacksIt_FilterOff_KeepsTheGlycopeptide()
        {
            // With OxoniumIonFilt unchecked the custom ion is only scored, never a gate.
            try
            {
                RegisterCustomMonosaccharideWithIon();

                var matches = FindNGlycanForKnownGlycopeptide(oxoniumIonFilter: false, customIonObserved: true);

                Assert.That(matches.Any(IsKnownGlycan), Is.True);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [TestCase(0.02, false)] // ratio 0.05 is above 2%: the ion is observed, and every glycan lacking it is cut
        [TestCase(0.1, true)]   // ratio 0.05 is below 10%: the ion counts as absent, so the glycopeptide is kept
        public void FindNGlycan_CustomIonAtFivePercent_EngineThresholdDecidesTheGate(double threshold, bool expectKnownGlycan)
        {
            // The engine's CustomOxoniumRelativeIntensityThreshold must reach the N-search gate: the same spectrum gives
            // opposite verdicts at two thresholds.
            try
            {
                RegisterCustomMonosaccharideWithIon();

                var matches = FindNGlycanForKnownGlycopeptide(oxoniumIonFilter: true, customIonObserved: true,
                    customIonIntensity: 5, customOxoniumThreshold: threshold);

                if (expectKnownGlycan)
                {
                    Assert.That(matches.Any(IsKnownGlycan), Is.True);
                }
                else
                {
                    Assert.That(matches, Is.Empty);
                }
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }
    }
}
