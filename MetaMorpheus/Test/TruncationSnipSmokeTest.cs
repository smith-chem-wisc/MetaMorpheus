using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test
{
    /// <summary>
    /// Phase 1.4 gate: exercise the mzLib snip utility once on Windows to confirm it produces a
    /// readable mzML. Snippets cut with this API become the real-data unit-test fixtures in later
    /// phases (01_Architecture.md #21, 02_Implementation.md Phase 4).
    /// </summary>
    [TestFixture]
    public class TruncationSnipSmokeTest
    {
        [Test]
        public void ExportSnipAsMzML_ProducesReadableSnippet()
        {
            string source = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "slicedTDYeast.mzML");
            MsDataFile file = MsDataFileReader.GetDataFile(source).LoadAllStaticData();
            int totalScans = file.GetAllScansList().Count;
            Assert.That(totalScans, Is.GreaterThan(0));

            string snipPath = file.ExportSnipAsMzML(1, totalScans);
            try
            {
                Assert.That(File.Exists(snipPath), Is.True);

                var snipped = MsDataFileReader.GetDataFile(snipPath).LoadAllStaticData().GetAllScansList();
                Assert.That(snipped.Count, Is.GreaterThan(0));
                Assert.That(snipped.First().OneBasedScanNumber, Is.EqualTo(1)); // snip renumbers scans 1-based
            }
            finally
            {
                if (File.Exists(snipPath))
                {
                    File.Delete(snipPath);
                }
            }
        }
    }
}
