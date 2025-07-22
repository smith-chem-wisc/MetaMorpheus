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
        // Set up file paths
        string testDir = TestContext.CurrentContext.TestDirectory;
        string fastaPath = Path.Combine(testDir, "Transcriptomics", "TestData", "16mer2.fasta");
        string xmlReferencePath = Path.Combine(testDir, "Transcriptomics", "TestData", "16mer2.xml");
        string spectraPath = Path.Combine(testDir, "Transcriptomics", "TestData", "RnaStandard_Subset.mzML");
        string outputDir = Path.Combine(testDir, "GptmdTask_Rna_AddsExpectedMethylationMods");
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

        // Parse the output XML and the reference XML
        var outputDoc = XDocument.Load(gptmdXml);
        var refDoc = XDocument.Load(xmlReferencePath);

        // Get the entry node
        var outputEntry = outputDoc.Descendants("entry").FirstOrDefault();
        var refEntry = refDoc.Descendants("entry").FirstOrDefault();
        Assert.That(outputEntry, Is.Not.Null, "Output XML should have an entry.");
        Assert.That(refEntry, Is.Not.Null, "Reference XML should have an entry.");

        // Check accession and sequence
        string outAcc = outputEntry.Element("accession")?.Value;
        string refAcc = refEntry.Element("accession")?.Value;
        Assert.That(outAcc, Is.EqualTo(refAcc));
        string outSeq = outputEntry.Element("sequence")?.Value;
        string refSeq = refEntry.Element("sequence")?.Value;
        Assert.That(outSeq, Is.EqualTo(refSeq));

        // Check mod positions and types
        var outMods = outputEntry.Elements("feature")
            .Where(f => (string)f.Attribute("type") == "modified residue")
            .Select(f => new
            {
                Description = (string)f.Attribute("description"),
                Position = int.Parse(f.Element("location").Element("position").Attribute("position").Value)
            }).OrderBy(f => f.Position).ToList();

        var refMods = refEntry.Elements("feature")
            .Where(f => (string)f.Attribute("type") == "modified residue")
            .Select(f => new
            {
                Description = (string)f.Attribute("description"),
                Position = int.Parse(f.Element("location").Element("position").Attribute("position").Value)
            }).OrderBy(f => f.Position).ToList();

        Assert.That(outMods.Count, Is.EqualTo(refMods.Count), "Number of mods should match reference.");
        for (int i = 0; i < refMods.Count; i++)
        {
            Assert.That(outMods[i].Position, Is.EqualTo(refMods[i].Position), $"Mod position {i} should match.");
            Assert.That(outMods[i].Description, Is.EqualTo(refMods[i].Description), $"Mod description {i} should match.");
        }
    }
}
