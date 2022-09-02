using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using GuiFunctions;
using IO.MzML;
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot;
using OxyPlot.Wpf;

namespace Test
{
    public class AAAChimeraSpectraGraphGenerator
    {
        [Test]
        public static void GenerateChimericSpectraPlot()
        {
            string outpath = @"C:\Users\Nic\Downloads\chimeraImage.png";
            string proteoformPath = @"D:\Projects\Top Down MetaMorpheus\For paper KHB\Cali_MOxAndBioMetArtModsGPTMD_Search\Task2-SearchTask\AllProteoforms.psmtsv";
            string spectraPath = @"D:\Projects\Top Down MetaMorpheus\For paper KHB\Calibrated Spectra\FXN5_tr1_032017-calib.mzML";
            List<PsmFromTsv> psms = PsmTsvReader.ReadTsv(proteoformPath, out List<string> warnings);
            List<MsDataScan> scans = Mzml.LoadAllStaticData(spectraPath).GetAllScansList();

            List<PsmFromTsv> chimeras = psms
                .Where(p => p.Ms2ScanNumber == 1314 && p.QValue <= 0.01).ToList();
            MsDataScan chimericScan = scans.First(p => p.OneBasedScanNumber == 1314);

          
        }

    }
}
