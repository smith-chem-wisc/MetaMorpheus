using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Omics;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    public class TestRnaSearchEngine
    {
        public static RnaSearchParameters SearchParameters;
        public static CommonParameters CommonParameters;
        public static string SixmerFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "GUACUG_NegativeMode_Sliced.mzML");

        [OneTimeSetUp]
        public static void Setup()
        {
            SearchParameters = new RnaSearchParameters
            {
                DecoyType = DecoyType.Reverse,
                MassDiffAcceptorType = MassDiffAcceptorType.Custom,
                CustomMdac = "Custom interval [-5,5]",
                DisposeOfFileWhenDone = true
            };
            CommonParameters = new CommonParameters
            (
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
                digestionParams: new RnaDigestionParams()
            );
        }

        [Test]
        public static void FindsSimpleSixmer()
        {
            List<Modification> fixedMods = new();
            List<Modification> variableMods = new();
            var dataFile = MsDataFileReader.GetDataFile(SixmerFilePath);
            var ms2Scans = MetaMorpheusTask.GetMs2Scans(dataFile, SixmerFilePath, CommonParameters)
                .OrderBy(b => b.PrecursorMass)
                .ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance,
                SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            var osms = new SpectralMatch[ms2Scans.Length];

            List<IBioPolymer> targets = new() { new RNA("GUACUG"), };
            var test = new ClassicSearchEngine(osms, ms2Scans, variableMods, fixedMods, null, null, null, targets, massDiffAcceptor, CommonParameters, [], null, ["search"], false);
            test.Run();

            var oligoSpectralMatches = osms.Where(p => p != null)
                .OrderByDescending(p => p.Score).ToList();
            Assert.That(oligoSpectralMatches.Count, Is.GreaterThan(0), "No matches found for GUACUG sixmer.");
            var match = oligoSpectralMatches.First();

            Assert.That(match, Is.TypeOf<OligoSpectralMatch>(), "Match is not of type OligoSpectralMatch.");
            Assert.That(match.Score, Is.GreaterThan(22), "Score for GUACUG sixmer match is not greater than 22.");
            Assert.That(match.BaseSequence, Is.EqualTo("GUACUG"), "Base sequence does not match GUACUG.");
        }
    }
}
