using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using Proteomics.ProteolyticDigestion;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics.RootFinding;
using MzLibUtil.NoiseEstimation;
using static System.Net.WebRequestMethods;
using System.Runtime.InteropServices;
using System.Xml.Linq;
using ThermoFisher.CommonCore.Data.Business;
using EngineLayer.DIA;

//namespace EngineLayer
//{
//    public  class DIA_Function
//    {
//        //public (double,double) XYData { get; set; }
//        public  static void FindAllMzTracePeakCurves(List<MsDataScan> scanCollection,double minMZ,CommonParameters commonParameters)
//        {
//            double preRT = 0;

//        //Loop for each scan in the ScanCollection
//         //ArrayList<ForkJoinTask<ArrayList<PeakCurve>>> ftemp = new ArrayList<>();
//         //ForkJoinPool fjp=new ForkJoinPool(NoCPUs);
//         //int idx_end = scanCollection.GetScanNoArray(MSlevel).size();
//         //var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
//         var ms2Scans = scanCollection.Where(p => p.MsnOrder == 2).ToArray();
//         var ms3Scans = scanCollection.Where(p => p.MsnOrder == 3).ToArray();
//         int idx_end = ms2Scans.Length;

//         int[] ia = new int[idx_end + 1];
//         ia[0]=0;
//         for (int idx = 0; idx<idx_end; idx++) {
//                int scanNO = scanCollection[idx].OneBasedScanNumber;
//                MsDataScan sd =scanCollection[idx];
//                ia[idx + 1]=sd.MassSpectrum.Size+ia[idx];
//            }

//        Boolean[] included = new Boolean[ia[ia.Length - 1]];
//        //if (step == -1)
//        //    step = fjp.getParallelism()* 32;
//        long peakCurvesCount = 0;
//        for (int idx = 0; idx<idx_end; idx++) {
//            int scanNO = scanCollection.[idx].OneBasedScanNumber;
//            MsDataScan scanData = scanCollection[idx];

//                //If we are doing targeted peak detection and the RT of current scan is not in the range of targeted list, jump to the next scan 
//                //if (TargetedOnly && !FoundInInclusionRTList(scanData.RetentionTime)){
//                //continue;
//                 //}
//        if (idx == 0)
//        {
//            preRT = scanData.RetentionTime - 0.01;
//        }
//        for (int i = 0; i < scanData.MassSpectrum.Size; i++)
//        {
//            (double,double) peak = (scanData.MassSpectrum.XArray[i], scanData.MassSpectrum.YArray[i]);
//        ////If we are doing targeted peak detection and the RT and m/z of current peak is not in the range of targeted list, jump to the next peak 
//        //if (TargetedOnly && !FoundInInclusionMZList(scanData.RetentionTime, peak.getX()))
//        //{
//        //    continue;
//        //}

//        if (peak.Item1 < minMZ)
//        {
//            continue;
//        }

//    //Check if the current peak has been included in previously developed peak curves
//    //                if (!IncludedHashMap.contains(scanNO + "_" + peak.getX())) {//The peak hasn't been included
//     //int id_scanNO_peak = int_id(ia, idx, i);
//   // if (!included[id_scanNO_peak])
//    //{//The peak hasn't been included
//     //The current peak will be the starting peak of a new peak curve
//     //Add it to the hash table

//        //                    IncludedHashMap.add(scanNO + "_" + peak.getX());
//        //included[id_scanNO_peak] = true;

//        double startmz = peak.Item1;
//        double startint = peak.Item2;

//        //Find the maximum peak within PPM window as the starting peak
//        for (int j = i + 1; j < scanData.MassSpectrum.Size; j++)
//        {
//            (double,double) currentpeak = (scanData.MassSpectrum.XArray[j], scanData.MassSpectrum.YArray[j]); ;
//            //int id_scanNO_currentpeak = int_id(ia, idx, j);
//            //if (!included[id_scanNO_currentpeak])
//            //{
                
//                if (InstrumentParameter.CalcPPM(currentpeak.getX(), startmz) <= PPM)
//                {
//                    //included[id_scanNO_currentpeak] = true;
                  
//                    if (currentpeak.Item2 >= startint)
//                    {
//                        startmz = currentpeak.Item1;
//                        startint = currentpeak.Item2;
//                    }
//                }
//                else
//                {
//                    break;
//                }
//            //}
//       // }

//        //Initialize a new peak curve
//        PeakCurve Peakcurve = new PeakCurve();
//        //Add a background peak
//        Peakcurve.AddPeak(preRT, startmz, commonParameters.MinimumAllowedIntensityRatioToBasePeak.Value);
//        //Add the starting peak
//        Peakcurve.AddPeak(scanData.RetentionTime, startmz, startint);
//        Peakcurve.StartScan = scanNO;

//        int missedScan = 0;
//        float endrt = scanData.RetentionTime;
//        int endScan = scanData.ScanNum;
//        float bk = 0f;

//        //Starting from the next scan, find the following peaks given the starting peak
//        for (int idx2 = idx + 1; idx2 < scanCollection.GetScanNoArray(MSlevel).size() && (missedScan < parameter.NoMissedScan /*|| (TargetedOnly && Peakcurve.RTWidth()<parameter.MaxCurveRTRange)*/); idx2++)
//        {
//            int scanNO2 = scanCollection.GetScanNoArray(MSlevel).get(idx2);
//            ScanData scanData2 = scanCollection.GetScan(scanNO2);

//            endrt = scanData2.RetentionTime;
//            endScan = scanData2.ScanNum;
//            bk = scanData2.background;
//            float currentmz = 0f;
//            float currentint = 0f;

//            //If the scan is empty
//            if (scanData2.PointCount() == 0)
//            {
//                if (parameter.FillGapByBK)
//                {
//                    Peakcurve.AddPeak(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background);
//                }
//                missedScan++;
//                continue;
//            }

//            //Find the m/z index 
//            int mzidx = scanData2.GetLowerIndexOfX(Peakcurve.TargetMz);
//            for (int pkidx = mzidx; pkidx < scanData2.Data.size(); pkidx++)
//            {
//                XYData currentpeak = scanData2.Data.get(pkidx);
//                if (currentpeak.getX() < parameter.MinMZ)
//                {
//                    continue;
//                }
//                //Check if the peak has been included or not
//                final int int_id_scanNO2_currentpeak = int_id(ia, idx2, pkidx);
//                //                            if (!included.get(int_id_scanNO2_currentpeak)) {
//                if (!included[int_id_scanNO2_currentpeak])
//                {
//                    if (InstrumentParameter.CalcPPM(currentpeak.getX(), Peakcurve.TargetMz) > PPM)
//                    {
//                        if (currentpeak.getX() > Peakcurve.TargetMz)
//                        {
//                            break;
//                        }
//                    }
//                    else
//                    {
//                        //////////The peak is in the ppm window, select the highest peak
//                        included[int_id_scanNO2_currentpeak] = true;
//                        //                                    IncludedHashMap.add(scanNO2 + "_" + currentpeak.getX());
//                        if (currentint < currentpeak.getY())
//                        {
//                            currentmz = currentpeak.getX();
//                            currentint = currentpeak.getY();
//                        }
//                    }
//                }
//            }

//            //No peak in the PPM window has been found
//            if (currentmz == 0f)
//            {
//                if (parameter.FillGapByBK)
//                {
//                    Peakcurve.AddPeak(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background);
//                }
//                missedScan++;
//            }
//            else
//            {
//                missedScan = 0;
//                Peakcurve.AddPeak(scanData2.RetentionTime, currentmz, currentint);
//            }
//        }
//        Peakcurve.AddPeak(endrt, Peakcurve.TargetMz, bk);
//        Peakcurve.EndScan = endScan;

//        //First check if the peak curve is in targeted list
//        if (FoundInInclusionList(Peakcurve.TargetMz, Peakcurve.StartRT(), Peakcurve.EndRT()))
//        {
//            //                        LCMSPeakBase.UnSortedPeakCurves.add(Peakcurve);
//            ++peakCurvesCount;
//            ftemp.add(fjp.submit(new PeakCurveSmoothingUnit(Peakcurve, parameter)));
//            //Then check if the peak curve passes the criteria
//        }
//        else if (Peakcurve.GetRawSNR() > LCMSPeakBase.SNR && Peakcurve.GetPeakList().size() >= parameter.MinPeakPerPeakCurve + 2)
//        {
//            //                        LCMSPeakBase.UnSortedPeakCurves.add(Peakcurve);
//            ++peakCurvesCount;
//            ftemp.add(fjp.submit(new PeakCurveSmoothingUnit(Peakcurve, parameter)));
//        }
//        else
//        {
//            Peakcurve = null;
//        }
//    }
//}
//preRT = scanData.RetentionTime;
///** the if statement below does PeakCurveSmoothing() and ClearRawPeaks()
// */
//final boolean last_iter = idx + 1 == idx_end;
//if (ftemp.size() == step || last_iter)
//{
//    final List<ForkJoinTask< ArrayList < PeakCurve >>> ftemp_sublist_view =
//            last_iter ?
//            ftemp
//            : ftemp.subList(0, step / 2);
//    for (final Future< ArrayList < PeakCurve >> f : ftemp_sublist_view)
//    {
//        try { LCMSPeakBase.UnSortedPeakCurves.addAll(f.get()); }
//        catch (InterruptedException| ExecutionException e){ throw new RuntimeException(e); }
//}
//ftemp_sublist_view.clear();
//if (!last_iter && fjp.getActiveThreadCount() < fjp.getParallelism())
//{
//    //                    System.out.println("PeakCurveSmoothingUnit: fjp.getActiveThreadCount()\t"+fjp.getActiveThreadCount()+"\t"+step);
//    step *= 2;
//}
//            }
//        }
//        assert ftemp.isEmpty();
////System.out.print("PSM removed (PeakCurve generation):" + PSMRemoved );
//fjp.shutdown();

//try
//{
//    fjp.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//}
//catch (InterruptedException e)
//{
//    LogManager.getRootLogger().info("interrupted..");
//}
//int i = 1;
////Assign peak curve index
//for (PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves)
//{
//    peakCurve.Index = i++;
//}

//System.gc();
////        LogManager.getRootLogger().info(LCMSPeakBase.UnSortedPeakCurves.size() + " Peak curves found (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
//LogManager.getRootLogger().info(peakCurvesCount + " Peak curves found (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
//    }
        
//    private boolean FoundInInclusionRTList(float rt)
//{
//    return Math.abs(InclusionRT.Data.get(InclusionRT.GetClosetIndexOfX(rt)).getX() - rt) < parameter.MaxCurveRTRange;
//}
        
//    }
//}
