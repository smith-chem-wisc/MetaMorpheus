//using MassSpectrometry;
//using Proteomics.Fragmentation;
//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using System.Windows.Controls;

//namespace EngineLayer
//{
//    public class ParentChildSpectrumMatchPlot : PeptideSpectrumMatchPlot
//    {
//        public ParentChildScanPlotsView ParentChildScanPlotsView { get; private set; }

//        public ParentChildSpectrumMatchPlot(OxyPlot.Wpf.PlotView parentPlotView, Canvas sequenceDrawingCanvas, PsmFromTsv psm, MsDataScan scan,
//            DynamicDataConnection dataFile, ParentChildScanPlotsView parentChildScanPlotsView) : base(parentPlotView, sequenceDrawingCanvas, psm, scan)
//        {
//            ParentChildScanPlotsView = parentChildScanPlotsView;
//            PlotParentAndChildScans(dataFile);
//        }

//        protected void PlotParentAndChildScans(DynamicDataConnection dataFile)
//        {
//            // clear old parent/child scans
//            ParentChildScanPlotsView.Plots.Clear();

//            // draw parent scan
//            MsDataScan parentScan = dataFile.GetOneBasedScanFromDynamicConnection(SpectrumMatch.Ms2ScanNumber);

//            string parentAnnotation = "Scan: " + parentScan.OneBasedScanNumber
//                    + " Dissociation Type: " + parentScan.DissociationType
//                    + " MsOrder: " + parentScan.MsnOrder
//                    + " Selected Mz: " + parentScan.SelectedIonMZ.Value.ToString("0.##")
//                    + " Retention Time: " + parentScan.RetentionTime.ToString("0.##");

//            ParentChildScanPlotsView.AddNewRow(this, parentAnnotation, null);

//            // draw child scans
//            HashSet<int> scansDrawn = new HashSet<int>();
//            var allChildScanMatchedIons = SpectrumMatch.ChildScanMatchedIons;

//            if (SpectrumMatch.BetaPeptideChildScanMatchedIons != null)
//            {
//                allChildScanMatchedIons = allChildScanMatchedIons.Concat(SpectrumMatch.BetaPeptideChildScanMatchedIons).ToDictionary(p => p.Key, q => q.Value);
//            }

//            foreach (var childScanMatchedIons in allChildScanMatchedIons)
//            {
//                int childScanNumber = childScanMatchedIons.Key;

//                if (scansDrawn.Contains(childScanNumber))
//                {
//                    continue;
//                }

//                scansDrawn.Add(childScanNumber);

//                List<MatchedFragmentIon> matchedIons = childScanMatchedIons.Value;

//                MsDataScan childScan = dataFile.GetOneBasedScanFromDynamicConnection(childScanNumber);

//                string childAnnotation = "Scan: " + childScan.OneBasedScanNumber
//                    + " Dissociation Type: " + childScan.DissociationType
//                    + " MsOrder: " + childScan.MsnOrder
//                    + " Selected Mz: " + childScan.SelectedIonMZ.Value.ToString("0.##")
//                    + " RetentionTime: " + childScan.RetentionTime.ToString("0.##");

//                Canvas canvas = new Canvas();
//                OxyPlot.Wpf.PlotView plotView = new OxyPlot.Wpf.PlotView();
//                var psmPlot = new PeptideSpectrumMatchPlot(plotView, canvas, SpectrumMatch, childScan, matchedIons);
//                psmPlot.Model.Title = null;
//                psmPlot.Model.Subtitle = null;

//                ParentChildScanPlotsView.AddNewRow(psmPlot, childAnnotation, canvas);
//            }
//        }
//    }
//}
