using NUnit.Framework;
using EngineLayer;
using Omics.Modifications;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Transcriptomics;
using UsefulProteomicsDatabases.Transcriptomics;
using UsefulProteomicsDatabases;
using Omics;
using TaskLayer;

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
}