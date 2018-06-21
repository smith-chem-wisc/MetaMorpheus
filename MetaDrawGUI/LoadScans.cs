using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;
using System.IO;
using MzLibUtil;
using Chemistry;
using TaskLayer;

namespace MetaDrawGUI
{
    class LoadScans
    {
        private string outputFolder;
        private string currentRawDataFilename;

        public LoadScans(string startingRawFilename, string outputFolder)
        {
            this.outputFolder = outputFolder;
            currentRawDataFilename = startingRawFilename;
        }

        public Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass { get; set; }        

        public MsDataFile Run()
        {
            SearchParameters searchParameters = new SearchParameters();

            CommonParameters commonParameters = new CommonParameters();

            MyFileManager myFileManager = new MyFileManager(searchParameters.DisposeOfFileWhenDone);

            var msDataScans = myFileManager.LoadFile(currentRawDataFilename, commonParameters.TopNpeaks, commonParameters.MinRatio, commonParameters.TrimMs1Peaks, commonParameters.TrimMsMsPeaks);

            return msDataScans;
        }


       
    }
}
