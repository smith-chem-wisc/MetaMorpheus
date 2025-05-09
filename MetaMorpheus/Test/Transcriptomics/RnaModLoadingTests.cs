using NUnit.Framework;
using EngineLayer;

namespace Test.Transcriptomics;

[TestFixture]
public class RnaModLoadingTests
{
    [Test]
    public static void ModsLoadWithoutCrash()
    {
        // This test is to ensure that the mods load without crashing.
        // It does not check for correctness of the mods.

        Assert.That(GlobalVariables.AllRnaModsKnown, Is.Not.Null);
        Assert.That(GlobalVariables.AllRnaModsKnown, Is.Not.Empty);

        Assert.That(GlobalVariables.AllRnaModTypesKnown, Is.Not.Null);
        Assert.That(GlobalVariables.AllRnaModTypesKnown, Is.Not.Empty);

        Assert.That(GlobalVariables.AllRnaModsKnownDictionary, Is.Not.Null);
        Assert.That(GlobalVariables.AllRnaModsKnownDictionary, Is.Not.Empty);
    }

    [Test]
    public static void RnaDbLoading_ModsLoadFromXml()
    {

    }
}