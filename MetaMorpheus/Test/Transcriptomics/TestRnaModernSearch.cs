using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Readers;
using TaskLayer;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    /// <summary>
    /// Modern (index-based) search over a nucleic acid database. The classic-search equivalents live in
    /// TestRnaSearchEngine and SearchTaskWithRna; these assert the indexed path reaches the same answer.
    /// </summary>
    [TestFixture]
    public class TestRnaModernSearch
    {
        private static readonly string SixmerFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "GUACUG_NegativeMode_Sliced.mzML");

        private static CommonParameters RnaCommonParameters => new CommonParameters(
            dissociationType: DissociationType.CID,
            deconvolutionMaxAssumedChargeState: -20,
            deconvolutionIntensityRatio: 3,
            deconvolutionMassTolerance: new PpmTolerance(20),
            precursorMassTolerance: new PpmTolerance(10),
            productMassTolerance: new PpmTolerance(20),
            scoreCutoff: 5,
            totalPartitions: 1,
            maxThreadsToUsePerFile: 1,
            doPrecursorDeconvolution: true,
            useProvidedPrecursorInfo: false,
            digestionParams: new RnaDigestionParams());

        /// <summary>
        /// The index used to refuse RNA outright: IndexingEngine threw "Not yet implemented for Rna
        /// Digestion" and SearchTask cast the biopolymer list to Protein. It now digests and fragments
        /// oligos like any other biopolymer.
        /// </summary>
        [Test]
        public static void IndexingEngine_IndexesOligos()
        {
            var commonParameters = RnaCommonParameters;
            List<IBioPolymer> targets = [new RNA("GUACUG")];

            var indexEngine = new IndexingEngine(targets, [], [], null, null, null, 0, DecoyType.None,
                commonParameters, null, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)indexEngine.Run();

            Assert.That(results.PeptideIndex, Is.Not.Empty, "no oligos were digested into the index");
            Assert.That(results.PeptideIndex.First(), Is.InstanceOf<OligoWithSetMods>());
            Assert.That(results.FragmentIndex.Count(bin => bin != null), Is.GreaterThan(0), "no fragments were binned");
        }

        /// <summary>
        /// Peptide ids must ascend within a fragment bin: the search binary-searches a bin by mass, and
        /// that only holds because the index is sorted by mass. True for oligos as much as for peptides.
        /// </summary>
        [Test]
        public static void OligoIndex_BinsAreSortedByOligoId()
        {
            List<IBioPolymer> targets = [new RNA("GUACUGGUACUGAAUUCC"), new RNA("CAUGCAUGCAUGAAUU")];

            var indexEngine = new IndexingEngine(targets, [], [], null, null, null, 0, DecoyType.Reverse,
                RnaCommonParameters, null, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)indexEngine.Run();

            int multiEntryBins = 0;
            foreach (var bin in results.FragmentIndex.Where(b => b != null))
            {
                for (int i = 1; i < bin.Count; i++)
                {
                    Assert.That(bin[i], Is.GreaterThan(bin[i - 1]), "a fragment bin is not ascending by oligo id");
                }
                if (bin.Count > 1)
                {
                    multiEntryBins++;
                }
            }

            Assert.That(multiEntryBins, Is.GreaterThan(0), "no multi-entry bins; the test proves nothing");
            Assert.That(results.PeptideIndex.Select(p => p.MonoisotopicMass), Is.Ordered, "index is not sorted by mass");
        }

        /// <summary>
        /// The payoff: modern search finds the same sixmer classic search finds, and reports it as an
        /// OligoSpectralMatch rather than mislabelling it a PSM.
        /// </summary>
        [Test]
        public static void ModernSearch_FindsSimpleSixmer()
        {
            var commonParameters = RnaCommonParameters;
            var searchParameters = new RnaSearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Custom,
                CustomMdac = "Custom interval [-5,5]",
            };

            var dataFile = MsDataFileReader.GetDataFile(SixmerFilePath);
            var ms2Scans = MetaMorpheusTask.GetMs2Scans(dataFile, SixmerFilePath, commonParameters)
                .OrderBy(b => b.PrecursorMass)
                .ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(commonParameters.PrecursorMassTolerance,
                searchParameters.MassDiffAcceptorType, searchParameters.CustomMdac);

            List<IBioPolymer> targets = [new RNA("GUACUG")];

            var indexEngine = new IndexingEngine(targets, [], [], null, null, null, 0, DecoyType.None,
                commonParameters, null, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            var osms = new SpectralMatch[ms2Scans.Length];
            new ModernSearchEngine(osms, ms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0,
                commonParameters, null, massDiffAcceptor, 0, ["modern rna search"]).Run();

            var matches = osms.Where(p => p != null).OrderByDescending(p => p.Score).ToList();

            Assert.That(matches, Is.Not.Empty, "modern search found no matches for the GUACUG sixmer");

            var match = matches.First();
            Assert.That(match, Is.TypeOf<OligoSpectralMatch>(), "match is not an OligoSpectralMatch");
            Assert.That(match.BaseSequence, Is.EqualTo("GUACUG"));
            Assert.That(match.Score, Is.GreaterThan(22), "modern score is below what classic search reaches");
        }
    }
}
