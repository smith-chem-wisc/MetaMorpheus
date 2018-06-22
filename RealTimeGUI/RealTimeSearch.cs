using System;
using System.Threading;
using System.Collections.Generic;

using Thermo.Interfaces.ExactiveAccess_V1;
using Thermo.Interfaces.InstrumentAccess_V1.MsScanContainer;
using IMsScan = Thermo.Interfaces.InstrumentAccess_V2.MsScanContainer.IMsScan;

using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using TaskLayer;
using EngineLayer;

namespace RealTimeGUI
{
    public class RealTimeSearch
    {

        private string outputFolder;
        private List<DbForTask> currentXmlDbFilenameList;


        public RealTimeSearch(List<DbForTask> startingXmlDbFilenameList, string outputFolder)
        {
            this.outputFolder = outputFolder;

            currentXmlDbFilenameList = startingXmlDbFilenameList;
        }

        public static event EventHandler<StringEventArgs> WarnHandler;

        public void Run()
        {
            var stopWatch = new Stopwatch();
            stopWatch.Start();

            StringBuilder allResultsText = new StringBuilder();

            {
                if (!currentXmlDbFilenameList.Any())
                {
                    Warn("Cannot proceed. No protein database files selected.");
                    return;
                }

                // Actual task running code


            }
            stopWatch.Stop();
            var resultsFileName = Path.Combine(outputFolder, "allResults.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine("MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion);
                file.WriteLine("Total time: " + stopWatch.Elapsed);
                file.Write(allResultsText.ToString());
            }
        }

        #region Private Methods

        private void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
        }

        #endregion Private Methods

    }
}
