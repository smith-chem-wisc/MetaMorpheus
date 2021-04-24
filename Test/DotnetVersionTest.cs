using EngineLayer;
using NUnit.Framework;
using System;
using System.IO;

namespace Test
{
    [TestFixture]
    internal class DotnetVersionTest
    {
        [Test]
        public void FrameworkMatchesReadme()
        {
            string readmeDirectory = Directory.GetParent(Environment.CurrentDirectory).Parent.Parent.Parent.FullName;
            string dotnetVersionInReadme = DotnetVersion.GetDotnetVersionFromReadme(File.ReadAllText(Path.Combine(readmeDirectory, "README.md")));
            Assert.IsTrue(DotnetVersion.IsSameAsVersion(dotnetVersionInReadme));
        }
    }
}