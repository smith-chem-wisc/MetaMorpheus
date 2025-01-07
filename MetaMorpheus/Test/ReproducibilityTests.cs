using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
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

        //public static async Task DrainResources(CancellationToken cancellationToken)
        //{
            //await Task.Run(() =>
            //{
            //    while(!cancellationToken.IsCancellationRequested)
            //    {
            //        // Tie up the CPU
            //        int[] unsortedArray = new int[1000000];
            //        Random random = new Random();
            //        for (int i = 0; i < unsortedArray.Length; i++)
            //        {
            //            unsortedArray[i] = random.Next(0, 1000000);
            //        }
            //        Array.Sort(unsortedArray);
            //    }
            //});
        //}

        //[Test]
        //public static void ReproducibilityTest()
        //{
        //    EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpGPTMD, out var testCase);
        //    string outputFolder = testCase.OutputDirectory;
        //    var allResultsFile = Path.Combine(outputFolder, "allResults.txt");
        //}

    }
}
