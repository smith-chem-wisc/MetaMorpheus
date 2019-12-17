using NUnit.Framework;
using System;
using System.Diagnostics;
using System.IO;

namespace Test
{
    [TestFixture]
    public static class TestCmd
    {
        [Test]
        //[Ignore("Ignored on AppVeyor")]
        public static void TestCommandLineMicroVignette()
        {
            Stopwatch s = new Stopwatch();
            s.Start();

            // run the micro vignette via command-line
            MetaMorpheusCommandLine.Program.Main(new string[] {
                "-v",
                "-o" + Path.Combine(TestContext.CurrentContext.TestDirectory, @"CommandLineMicroVignette") } );

            s.Stop();
            Console.WriteLine("Command-line microvignette took: " + s.ToString());
        }
    }
}
