using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using TaskLayer;

namespace Test
{
    [TestFixture]
    class XLSearchOutputTest
    {
        [Test]
        public void WriteTsvTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLTestData\OutputTest1");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLTestData\BSA_DSS_23747.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLTestData\BSA.fasta");

            XLSearchTask xLSearch = new XLSearchTask();
            xLSearch.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");

            var resultsPath = File.ReadAllLines(Path.Combine(outputFolder, @"allPsms.tsv"));
            string expected = "C:\\MetaMorpheus\\Test\\bin\\Debug\\netcoreapp2.0\\XLTestData\\BSA_DSS_23747.mzML\t9\t808.380676269531\t4\t3229.49359921061\tCross\t\tDECOY_3336842(0)\tTVEVFEAKPFK(8)\tTVEVFEAKPFK\t1293.696888652\t4.01633312976919\t4.01633312976919\t30\t\t3336842(0)\tPDPNTLCDEFKADEK(11)\tPDPNTLC[Common Fixed:Carbamidomethyl on C]DEFKADEK\t1777.782879933\t4.01654288568969\t4.01654288568969\t42\t\t8.03287601545888\t158.013830625609\t-\t.";
            Assert.That(resultsPath[1].Equals(expected));

        }
    }
}
