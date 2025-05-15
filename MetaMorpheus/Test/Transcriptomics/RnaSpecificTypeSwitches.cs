using System;
using System.Collections.Generic;
using System.IO;
using EngineLayer;
using EngineLayer.Indexing;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics;

[TestFixture]
public class RnaSpecificTypeSwitches
{
    private static IDigestionParams BottomUpParams => new DigestionParams();
    private static IDigestionParams TopDownParams => new DigestionParams("top-down");
    private static IDigestionParams OligoParams => new RnaDigestionParams();

    [Test]
    [TestCase(AnalyteType.Peptide, nameof(BottomUpParams))]
    [TestCase(AnalyteType.Proteoform, nameof(TopDownParams))]
    [TestCase(AnalyteType.Oligo, nameof(OligoParams))]
    public static void MetaMorpheusTask_DetectAnalyteType(AnalyteType expected, string digestionParamsName)
    {
        // Get the right params. 
        IDigestionParams digestionParams = digestionParamsName switch
        {
            nameof(BottomUpParams) => BottomUpParams,
            nameof(TopDownParams) => TopDownParams,
            nameof(OligoParams) => OligoParams,
            _ => throw new ArgumentException("Invalid digestionParamsName", nameof(digestionParamsName))
        };

        if (GlobalVariables.AnalyteType == expected)
            GlobalVariables.AnalyteType = AnalyteType.Oligo;

        var commonParams = new CommonParameters(digestionParams: digestionParams);
        MetaMorpheusTask.DetermineAnalyteType(commonParams);
        Assert.That(GlobalVariables.AnalyteType, Is.EqualTo(expected), $"Expected {expected} but got {GlobalVariables.AnalyteType}");
    }

    [Test]
    public static void IndexingEngine_ThrowsOnRNA()
    {
        var indexEngine = new IndexingEngine([], [], [], null, null, null,
            1, DecoyType.Reverse, new CommonParameters(digestionParams: new RnaDigestionParams()), [], 14, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

        Assert.Throws<MetaMorpheusException>(() => indexEngine.Run());
    }

    [Test]
    public static void IndexingEngine_ToString_HasInitiatorMethionineWhenProteinMode()
    {
        var protEngine = new IndexingEngine([], [], [], null, null, null,
            1, DecoyType.Reverse, new CommonParameters(), [], 14, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
        var protString = protEngine.ToString();

        var oligoEngine = new IndexingEngine([], [], [], null, null, null,
            1, DecoyType.Reverse, new CommonParameters(digestionParams: new RnaDigestionParams()), [], 14, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
        var oligoString = oligoEngine.ToString();

        Assert.That(protString, Does.Contain("initiatorMethionineBehavior"));
        Assert.That(oligoString, Does.Not.Contain("initiatorMethionineBehavior"));

        Assert.That(protString, Does.Contain("specificProtease"));
        Assert.That(oligoString, Does.Not.Contain("specificProtease"));
    }
}
