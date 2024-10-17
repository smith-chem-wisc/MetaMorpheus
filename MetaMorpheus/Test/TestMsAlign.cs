using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using TaskLayer;

namespace Test
{
    internal class TestMsAlign
    {

        [Test]
        public static void TESTNAME()
        {
            var path =
                @"B:\Users\Nic\Chimeras\TopDown_Analysis\Jurkat\DeconResults\TopFD\02-17-20_jurkat_td_rep2_fract4_ms2.msalign";
            var fileManager = new MyFileManager(true);
            var dataFile = fileManager.LoadFile(path, new CommonParameters());

            dataFile.LoadAllStaticData();
            
            foreach (var scanWithMass in MetaMorpheusTask.GetMs2Scans(dataFile, path, new CommonParameters()))
            {
                int expectedLength = scanWithMass.TheScan.MassSpectrum.XArray.Length;

                Assert.That(scanWithMass.ExperimentalFragments.Length, Is.EqualTo(expectedLength));
            }
        }
    }
}
