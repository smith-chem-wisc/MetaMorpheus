using NUnit.Framework;
using UsefulProteomicsDatabases;

[SetUpFixture]
public class HarnessSetUp
{
    [OneTimeSetUp]
    public void Setup()
    {
        Loaders.LoadElements();
        EngineLayer.GlobalVariables.SetUpGlobalVariables();
    }
}
