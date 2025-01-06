using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Uses test cases found in EverythingRunnerEngineTestCase.cs
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ReproducibilityTests
    {
        public static Array GetTestCases() => Enum.GetValues(typeof(EverythingRunnerEngineTestCases));



    }
}
