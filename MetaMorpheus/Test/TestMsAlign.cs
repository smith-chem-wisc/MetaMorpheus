using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using EngineLayer;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Readers;
using TaskLayer;

namespace Test
{
    internal class TestMsAlign
    {

        [Test]
        public static void TestMsAlign_GetMs2Scans()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "JurkatTopDownRep2Fract1_ms2.msalign");
            var fileManager = new MyFileManager(true);
            var dataFile = fileManager.LoadFile(path, new CommonParameters());

            dataFile.LoadAllStaticData();
            
            foreach (var scanWithMass in MetaMorpheusTask.GetMs2Scans(dataFile, path, new CommonParameters()))
            {
                var scan = scanWithMass.TheScan;
                int expectedLength = scanWithMass.TheScan.MassSpectrum.XArray.Length;
                Assert.That(scanWithMass.PrecursorMonoisotopicPeakMz, Is.GreaterThan(0));
                Assert.That(scanWithMass.PrecursorCharge, Is.GreaterThan(0));
                Assert.That(scanWithMass.PrecursorIntensity, Is.GreaterThan(0));
                Assert.That(scanWithMass.PrecursorEnvelopePeakCount, Is.GreaterThan(0));
                Assert.That(scanWithMass.PrecursorFractionalIntensity, Is.GreaterThanOrEqualTo(-1));
                Assert.That(scanWithMass.FullFilePath, Is.EqualTo(path));

                Assert.That(scanWithMass.ExperimentalFragments.Length, Is.EqualTo(expectedLength));
                for (int i = 0; i < scanWithMass.ExperimentalFragments.Length; i++)
                {
                    Assert.That(scanWithMass.ExperimentalFragments[i], Is.Not.Null);
                    Assert.That(scanWithMass.ExperimentalFragments[i].Peaks, Is.Not.Null);
                    Assert.That(scanWithMass.ExperimentalFragments[i].Peaks.Count, Is.EqualTo(1));
                    Assert.That(scanWithMass.ExperimentalFragments[i].Peaks[0].mz, Is.EqualTo(scan.MassSpectrum.XArray[i].ToMz(scanWithMass.ExperimentalFragments[i].Charge)));
                    Assert.That(scanWithMass.ExperimentalFragments[i].Peaks[0].intensity, Is.EqualTo(scan.MassSpectrum.YArray[i]));
                }
            }
        }
    }
}
