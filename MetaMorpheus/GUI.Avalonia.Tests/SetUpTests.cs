using EngineLayer;
using NUnit.Framework;

namespace Test.Avalonia;

/// <summary>
/// Loads the modification, protease, crosslinker and glycan databases once for the whole assembly,
/// mirroring Test/SetUpTests.cs. Four fixtures used to do this per test, re-parsing ~4,300
/// modifications each time.
/// </summary>
[SetUpFixture]
public class SetUpTests
{
    [OneTimeSetUp]
    public void SetUp() => GlobalVariables.SetUpGlobalVariables();
}
