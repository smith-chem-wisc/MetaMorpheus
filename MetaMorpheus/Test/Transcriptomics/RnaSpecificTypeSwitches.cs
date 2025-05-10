using System;
using EngineLayer;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using Transcriptomics.Digestion;

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
}
