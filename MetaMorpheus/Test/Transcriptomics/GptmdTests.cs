using EngineLayer.Gptmd;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using Transcriptomics.Digestion;
using Transcriptomics;
using System.IO;
using System.Xml.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics;

public class GptmdTests
{
    [Test]
    [TestCase("GUACUG", 28.030300, 6)]
    [TestCase("GUACUG", 14.015650, 6)]
    public static void TestRunningGptmd_AllResiduesWithCombos(string rnaSequence, double modMassToAdd, int expectedGptmdMods)
    {
        var allResultingIdentifications = new List<SpectralMatch>();
        IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>()
                    { new Tuple<double, double>(14.015650, 14.015650) };
        Tolerance precursorMassTolerance = new PpmTolerance(10);
        var gptmdMods = GlobalVariables.AllRnaModsKnown.Where(p => p.IdWithMotif.Contains("Methyl"))
            .ToList();

        List<Modification> variableModifications = new List<Modification>();
        var commonParams = new CommonParameters(digestionParams: new RnaDigestionParams("top-down"),
            listOfModsVariable: new List<(string, string)>(),
            listOfModsFixed: new List<(string, string)>(),
            deconvolutionMaxAssumedChargeState: -12);
        var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
        fsp.Add(("filepath", commonParams));

        // run the engine where nothing will be found
        var engine = new GptmdEngine(allResultingIdentifications, gptmdMods, combos,
            new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } },
            new CommonParameters(), fsp, new List<string>(), null);
        var gptmdResults = (GptmdResults)engine.Run();
        Assert.That(gptmdResults.Mods.Count, Is.EqualTo(0));

        // set up an oligo 
        var rna = new RNA(rnaSequence, "accession");
        var digestedOligo = rna.Digest(commonParams.DigestionParams, new(), variableModifications).First();
        Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 },
            new double[] { 1 }, false), 0, 1, true, Polarity.Positive,
            double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN,
            null, null, "scan=1", double.NaN, null,
            null, double.NaN, null, DissociationType.AnyActivationType,
            0),
            (digestedOligo.MonoisotopicMass + modMassToAdd).ToMz(1),
        1, "filepath", new CommonParameters());

        SpectralMatch newPsm = new OligoSpectralMatch(digestedOligo, 0, 0, 0, scan, commonParams, new List<MatchedFragmentIon>());
        newPsm.SetMs2Scan(scan.TheScan);


        newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
        allResultingIdentifications.Add(newPsm);

        engine = new GptmdEngine(allResultingIdentifications, gptmdMods, combos,
            new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(),
            null, new List<string>(), null);
        gptmdResults = (GptmdResults)engine.Run();
        Assert.That(gptmdResults.Mods.Count, Is.EqualTo(1));

        if (expectedGptmdMods == 0)
            Assert.That(gptmdResults.Mods.Count, Is.EqualTo(0));
        else
            Assert.That(gptmdResults.Mods["accession"].Count, Is.EqualTo(expectedGptmdMods));
    }

    [Test]
    [TestCase("GUACUG", 28.030300, 2)]
    [TestCase("GUACUG", 14.015650, 2)]
    public static void TestRunningGptmd_FewResiduesWithCombos(string rnaSequence, double modMassToAdd, int expectedGptmdMods)
    {
        var allResultingIdentifications = new List<SpectralMatch>();
        IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>()
                    { new Tuple<double, double>(14.015650, 14.015650) };
        Tolerance precursorMassTolerance = new PpmTolerance(10);
        var gptmdMods = GlobalVariables.AllRnaModsKnown.Where(p => p.IdWithMotif.Contains("Methylation on G"))
            .ToList();

        List<Modification> variableModifications = new List<Modification>();
        var commonParams = new CommonParameters(digestionParams: new RnaDigestionParams("top-down"),
            listOfModsVariable: new List<(string, string)>(),
            listOfModsFixed: new List<(string, string)>(),
            deconvolutionMaxAssumedChargeState: -12);
        var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
        fsp.Add(("filepath", commonParams));

        // run the engine where nothing will be found
        var engine = new GptmdEngine(allResultingIdentifications, gptmdMods, combos,
            new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } },
            new CommonParameters(), fsp, new List<string>(), null);
        var gptmdResults = (GptmdResults)engine.Run();
        Assert.That(gptmdResults.Mods.Count, Is.EqualTo(0));

        // set up an oligo 
        var rna = new RNA(rnaSequence, "accession");
        var digestedOligo = rna.Digest(commonParams.DigestionParams, new(), variableModifications).First();
        Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 },
            new double[] { 1 }, false), 0, 1, true, Polarity.Positive,
            double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN,
            null, null, "scan=1", double.NaN, null,
            null, double.NaN, null, DissociationType.AnyActivationType,
            0, null),
            (digestedOligo.MonoisotopicMass + modMassToAdd).ToMz(1),
        1, "filepath", new CommonParameters());

        SpectralMatch newPsm = new OligoSpectralMatch(digestedOligo, 0, 0, 0, scan, commonParams, new List<MatchedFragmentIon>());
        newPsm.SetMs2Scan(scan.TheScan);


        newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
        allResultingIdentifications.Add(newPsm);

        engine = new GptmdEngine(allResultingIdentifications, gptmdMods, combos,
            new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(),
            null, new List<string>(), null);
        gptmdResults = (GptmdResults)engine.Run();
        Assert.That(gptmdResults.Mods.Count, Is.EqualTo(1));

        if (expectedGptmdMods == 0)
            Assert.That(gptmdResults.Mods.Count, Is.EqualTo(0));
        else
            Assert.That(gptmdResults.Mods["accession"].Count, Is.EqualTo(expectedGptmdMods));
    }

    [Test]
    [TestCase("GUACUGAUGAUAUAYAU", 14.015650, 1)]
    [TestCase("GUACUGAUGACUAUUAAAUA", 14.015650, 2)]
    public static void TestRunningGptmd_FewResidues_NoCombos(string rnaSequence, double modMassToAdd, int expectedGptmdMods)
    {
        var allResultingIdentifications = new List<SpectralMatch>();
        IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
        Tolerance precursorMassTolerance = new PpmTolerance(10);
        var gptmdMods = GlobalVariables.AllRnaModsKnown.Where(p => p.IdWithMotif.Contains("Methylation on C"))
            .ToList();

        List<Modification> variableModifications = new List<Modification>();
        var commonParams = new CommonParameters(digestionParams: new RnaDigestionParams("top-down"),
            listOfModsVariable: new List<(string, string)>(),
            listOfModsFixed: new List<(string, string)>(),
            deconvolutionMaxAssumedChargeState: -12);
        var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
        fsp.Add(("filepath", commonParams));

        // run the engine where nothing will be found
        var engine = new GptmdEngine(allResultingIdentifications, gptmdMods, combos,
            new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } },
            new CommonParameters(), fsp, new List<string>(), null);
        var gptmdResults = (GptmdResults)engine.Run();
        Assert.That(gptmdResults.Mods.Count, Is.EqualTo(0));

        // set up an oligo 
        var rna = new RNA(rnaSequence, "accession");
        var digestedOligo = rna.Digest(commonParams.DigestionParams, new(), variableModifications).First();
        Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 },
            new double[] { 1 }, false), 0, 1, true, Polarity.Positive,
            double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN,
            null, null, "scan=1", double.NaN, null,
            null, double.NaN, null, DissociationType.AnyActivationType,
            0, null),
            (digestedOligo.MonoisotopicMass + modMassToAdd).ToMz(1),
        1, "filepath", new CommonParameters());

        SpectralMatch newPsm = new OligoSpectralMatch(digestedOligo, 0, 0, 0, scan, commonParams, new List<MatchedFragmentIon>());
        newPsm.SetMs2Scan(scan.TheScan);


        newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
        allResultingIdentifications.Add(newPsm);

        engine = new GptmdEngine(allResultingIdentifications, gptmdMods, combos,
            new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(),
            null, new List<string>(), null);
        gptmdResults = (GptmdResults)engine.Run();
        Assert.That(gptmdResults.Mods.Count, Is.EqualTo(1));

        if (expectedGptmdMods == 0)
            Assert.That(gptmdResults.Mods.Count, Is.EqualTo(0));
        else
            Assert.That(gptmdResults.Mods["accession"].Count, Is.EqualTo(expectedGptmdMods));
    }

    [Test]
    public void GptmdTask_Rna_AddsExpectedMethylationMods()
    {
        // Arrange
        string testDir = TestContext.CurrentContext.TestDirectory;
        string fastaPath = Path.Combine(testDir, "Transcriptomics", "TestData", "22mer-Am.fasta");
        string spectraPath = Path.Combine(testDir, "Transcriptomics", "TestData", "22merMethyl_Snip.mzML");
        string outputDir = Path.Combine(testDir, "GptmdTask_Rna_AddsExpectedMethylationMods");

        if (Directory.Exists(outputDir))
            Directory.Delete(outputDir, true); // Clean up previous test run
        Directory.CreateDirectory(outputDir);

        // Set up GptmdTask for RNA
        var gptmdTask = new GptmdTask
        {
            CommonParameters = new CommonParameters(
                digestionParams: new RnaDigestionParams("top-down"),
                dissociationType: DissociationType.CID,
                deconvolutionMaxAssumedChargeState: -20,
                deconvolutionIntensityRatio: 3,
                deconvolutionMassTolerance: new PpmTolerance(5),
                precursorMassTolerance: new PpmTolerance(10),
                productMassTolerance: new PpmTolerance(20),
                scoreCutoff: 5,
                totalPartitions: 1,
                maxThreadsToUsePerFile: 1,
                doPrecursorDeconvolution: true,
                useProvidedPrecursorInfo: false
            ),
            GptmdParameters = new GptmdParameters
            {
                ListOfModsGptmd = GlobalVariables.AllRnaModsKnown
                    .Where(m => m.IdWithMotif.Contains("Methylation"))
                    .Select(m => (m.ModificationType, m.IdWithMotif)).ToList()
            }
        };

        // Prepare DbForTask and run
        var dbs = new List<DbForTask> { new DbForTask(fastaPath, false) };
        var rawFiles = new List<string> { spectraPath };
        var results = gptmdTask.RunTask(outputDir, dbs, rawFiles, "test");

        // Find the output GPTMD XML file
        string gptmdXml = Directory.GetFiles(outputDir, "*GPTMD.xml").FirstOrDefault();
        Assert.That(gptmdXml, Is.Not.Null.And.Not.Empty, "GPTMD output XML should exist.");

        // load in the GPTMD xml file
        var rna = RnaDbLoader.LoadRnaXML(gptmdXml, true, DecoyType.None, false, GlobalVariables.AllRnaModsKnown, [], out _);

        // check that the expected modifications are present
        int expectedModIndex = 9;
        string expectedMod = "Methylation on A";

        var modifications = rna.First().OneBasedPossibleLocalizedModifications;
        Assert.That(modifications, Is.Not.Null.And.Not.Empty, "Modifications should not be null or empty.");
        Assert.That(modifications, Does.ContainKey(expectedModIndex), $"Modifications should contain index {expectedModIndex}.");
        var modList = modifications[expectedModIndex];
        Assert.That(modList, Is.Not.Null.And.Not.Empty, "Modifications list should not be null or empty.");
        Assert.That(modList.Select(p => p.IdWithMotif), Does.Contain(expectedMod), $"Modifications should contain expected modification: {expectedMod}.");

        if (Directory.Exists(outputDir))
            Directory.Delete(outputDir, true);
    }
}
