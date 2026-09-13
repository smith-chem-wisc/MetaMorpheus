using System;
using EngineLayer;
using NUnit.Framework;
using UsefulProteomicsDatabases;

// Namespace Test, not Test.GuiFunctionsCore: a [SetUpFixture] applies to its own namespace and those
// beneath it, and the linked tests live in Test.GuiTests and Test.MetaDraw. Getting this wrong is
// silent - 175 of 187 tests still passed, and only the twelve that actually read global state failed.
namespace Test;

/// <summary>
/// The minimum of Test/SetUpTests.cs that applies here: current directory, elements, and the global
/// variables every view model reads. Deliberately not a copy of the whole thing - the results-writing
/// and MessageBoxHelper parts of that fixture belong to the Windows-only suite.
/// </summary>
[SetUpFixture]
public class SetUpCoreTests
{
    [OneTimeSetUp]
    public static void GlobalSetup()
    {
        Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        Loaders.LoadElements();
        GlobalVariables.SetUpGlobalVariables();
    }
}
