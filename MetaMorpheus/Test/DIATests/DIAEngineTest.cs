using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;
using TaskLayer;
using System.IO;

namespace Test.DIATests
{
    public class DIAEngineTest
    {
        [Test]
        public static void DIAScanWindowMapTest()
        {
            string dataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData\\DIA\\18300_REP2_500ng_HumanLysate_SWATH_1_RT25.63-25.81.mzML");
            var myFileManager = new MyFileManager(true);
            var dataFile = myFileManager.LoadFile(dataPath, new CommonParameters());
            var allMs2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var diaScanWindowMap = DIAEngine.ConstructMs2Groups(allMs2Scans);
            Assert.That(diaScanWindowMap.Count, Is.EqualTo(34));
            foreach(var group in diaScanWindowMap)
            {
                Assert.That(group.Value.Count, Is.EqualTo(3));
            }
        }

        [Test]
        public static void PseudoMs2ScansTest()
        {
            string dataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData\\DIA\\18300_REP2_500ng_HumanLysate_SWATH_1_RT25.63-25.81.mzML");
            var myFileManager = new MyFileManager(true);
            var dataFile = myFileManager.LoadFile(dataPath, new CommonParameters());
            var diaParams = new DIAparameters();
            var pseudoMs2Scans = DIAEngine.GetPseudoMs2Scans(dataFile, new CommonParameters(), diaParams);
            Assert.That(pseudoMs2Scans.Count, Is.GreaterThan(0));
        }
    }
}