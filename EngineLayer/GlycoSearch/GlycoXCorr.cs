using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;

namespace EngineLayer.GlycoSearch
{
    public static class GlycoXCorr
    {

        private static void BinIons(Ms2ScanWithSpecificMass scan, double binSize, double invBinSize, double binOffset,
            ref int xCorrArraySize, ref int kojakBins, ref double[] pdCorrelationData_XArray, ref float[] pdCorrelationData_YArray,
            ref int iHighestIon, ref double dHighestIntensity)
        {


            //iHighestIon = 0;
            //dHighestIntensity = 0;


            //xCorrArraySize = (int)((scan.TheScan.MassSpectrum.XArray.Last() + 100.0) / binSize);
            //kojakBins = (int)(scan.TheScan.MassSpectrum.XArray.Last() + 100.0);

            //pdCorrelationData_XArray = new double[xCorrArraySize];
            //pdCorrelationData_YArray = new double[xCorrArraySize];

            //double invBinSize = 1.0 / binSize;
            double dPrecursor = scan.PrecursorMass;
            int i = 0;
            while (i < scan.TheScan.MassSpectrum.XArray.Length)
            {
                double dIon = scan.TheScan.MassSpectrum.XArray[i];
                double dIntensity = scan.TheScan.MassSpectrum.YArray[i];
                i++;

                if (dIntensity > 0.0)
                {
                    if (dIon < (dPrecursor + 50.0))
                    {
                        int iBinIon = (int)(dIon * invBinSize + binOffset);
                        dIntensity = Math.Sqrt(dIntensity);
                        if (iBinIon > iHighestIon)
                        {
                            iHighestIon = iBinIon;
                        }
                        if (iBinIon < xCorrArraySize && dIntensity > pdCorrelationData_YArray[iBinIon])
                        {
                            pdCorrelationData_YArray[iBinIon] = (float)dIntensity;
                            pdCorrelationData_XArray[iBinIon] = dIon;
                            if (pdCorrelationData_YArray[iBinIon] > dHighestIntensity)
                            {
                                dHighestIntensity = pdCorrelationData_YArray[iBinIon];
                            }
                        }

                    }
                }
            }
        }

        private static void MakeCorrData(double scale, ref int iHighestIon, ref double dHighestIntensity, ref float[] pdCorrelationData_YArray, ref int xCorrArraySize, ref double[] pdTempRawData)
        {
            int iNumWindows = 10;
            int iWindowSize = (int)Math.Ceiling((double)iHighestIon / iNumWindows);

            double[] dMaxWindowInten = new double[10];

            double dMaxOverallInten = 0.0;
            double dTmp1 = 1.0;
            double dTmp2;

            if (dHighestIntensity > 0.000001)
            {
                dTmp1 = 100.0 / dHighestIntensity;
            }

            int x = 0;
            int c = 0;
            for (int i = 0; i < pdCorrelationData_YArray.Length; i++)
            {
                dTmp2 = pdCorrelationData_YArray[i] * dTmp1;
                pdTempRawData[i] = dTmp2;
                pdCorrelationData_YArray[i] = 0;

                if (x < iNumWindows)
                {
                    if (dMaxWindowInten[x] < dTmp2)
                    {
                        dMaxWindowInten[x] = dTmp2;
                    }
                    c++;
                    if (c == iWindowSize)
                    {
                        c = 0;
                        x++;
                    }
                }
            }

            dMaxOverallInten = 100;
            dTmp2 = 0.05 * dMaxOverallInten;
            int iBin;
            for (int i = 0; i < iNumWindows; i++)
            {
                if (dMaxWindowInten[i] > 0.0)
                {
                    dTmp1 = scale / dMaxWindowInten[i];
                    for (int ii = 0; ii < iWindowSize; ii++)
                    {
                        iBin = i * iWindowSize + ii;
                        if (iBin < xCorrArraySize)
                        {
                            if (pdTempRawData[iBin] > dTmp2)
                            {
                                pdCorrelationData_YArray[iBin] = (float)(pdTempRawData[iBin] * dTmp1);
                            }
                        }
                    }

                }
            }

            var _test_pdCorr = new List<(int, double)>();
            for (int i = 0; i < xCorrArraySize; i++)
            {
                if (pdCorrelationData_YArray[i] > 0.01)
                {
                    _test_pdCorr.Add((i, pdCorrelationData_YArray[i]));
                }
            }
        }

        /// <summary>
        /// The function is for XCorr high res data, based on 'Kojak Crosslink KAnalysis.cpp'. 
        /// https://github.com/mhoopmann/kojak
        /// The binSize and binOffset is always setted. For high res data, binSize = 0.01, binOffset = 1.
        /// </summary>
        /// <param name="scan"></param>
        /// <param name="binSize"></param>
        /// <param name="binOffset"></param>
        /// <returns></returns>
        public static Dictionary<int, sbyte[]> kojakXCorr(Ms2ScanWithSpecificMass scan, double binSize, double binOffset)
        {
            int iHighestIon = 0;
            double dHighestIntensity = 0;

            int xCorrArraySize = (int)((scan.TheScan.MassSpectrum.XArray.Last() + 100.0) / binSize);
            int kojakBins = (int)(scan.TheScan.MassSpectrum.XArray.Last() + 100.0);

            double[] pdCorrelationData_XArray = new double[xCorrArraySize];
            float[] pdCorrelationData_YArray = new float[xCorrArraySize];

            double invBinSize = 1.0 / binSize;
            BinIons(scan, binSize, invBinSize, binOffset, ref xCorrArraySize, ref kojakBins, ref pdCorrelationData_XArray, ref pdCorrelationData_YArray, ref iHighestIon, ref dHighestIntensity);

            double[] pdTempRawData = new double[xCorrArraySize];
            double[] pdTmpFastXcorrData = new double[xCorrArraySize];
            float[] pfFastXcorrData = new float[xCorrArraySize];

            //char[] kojakSparseArray = new char*[kojakBins];
            //var kojakSparseArray = new List<sbyte[]>();
            var kojakSparseArray = new Dictionary<int, sbyte[]>();
            //for (int i = 0; i < kojakBins; i++)
            //{
            //    kojakSparseArray.Add(new List<sbyte[]>());
            //}

            MakeCorrData(50.0, ref iHighestIon, ref dHighestIntensity, ref pdCorrelationData_YArray, ref xCorrArraySize, ref pdTempRawData);

            double dSum = 0.0;
            for (int i = 0; i < 75; i++)
            {
                dSum += pdCorrelationData_YArray[i];
            }
            for (int i = 75; i < xCorrArraySize + 75; i++)
            {
                if (i < xCorrArraySize && pdCorrelationData_YArray[i] > 0)
                {
                    dSum += pdCorrelationData_YArray[i];
                }
                if (i >= 151 && pdCorrelationData_YArray[i - 151] > 0)
                {
                    dSum -= pdCorrelationData_YArray[i - 151];
                }
                pdTmpFastXcorrData[i - 75] = (dSum - pdCorrelationData_YArray[i - 75]) * 0.0066666667;
            }

            int xCorrSparseArraySize = 1;

            double dTmp0 = pdCorrelationData_YArray[0] - pdTmpFastXcorrData[0];
            double dTmp1 = pdCorrelationData_YArray[1] - pdTmpFastXcorrData[1];
            double dTmp2 = pdCorrelationData_YArray[2] - pdTmpFastXcorrData[2];

            pfFastXcorrData[0] = (float)((dTmp0 + dTmp2) * 0.5);
            pfFastXcorrData[1] = (float)(dTmp1 + (dTmp0 + dTmp2) * 0.5);
            for (int i = 2; i < xCorrArraySize - 1; i++)
            {
                dTmp0 = dTmp1;
                dTmp1 = dTmp2;
                dTmp2 = pdCorrelationData_YArray[i + 1] - pdTmpFastXcorrData[i + 1];
                pfFastXcorrData[i] = (float)(dTmp1 + (dTmp0 + dTmp2) * 0.5);
            }
            pfFastXcorrData[xCorrArraySize - 1] = (float)(dTmp2 + dTmp1 * 0.5);

            var _test_pfFast = new List<(int, double)>();
            for (int i = 0; i < xCorrArraySize; i++)
            {
                if (pfFastXcorrData[i] > 0.01)
                {
                    _test_pfFast.Add((i, pfFastXcorrData[i]));
                }
            }

            double dTmp;
            int iTmp;
            int jx;
            for (int i = 0; i < xCorrArraySize; i++)
            {
                if (pfFastXcorrData[i] > 0.5 || pfFastXcorrData[i] < -0.5)
                {
                    dTmp = binSize * i;
                    iTmp = (int)dTmp;
                    if (!kojakSparseArray.Keys.Contains(iTmp))
                    {
                        kojakSparseArray[iTmp] = new sbyte[(int)invBinSize + 1];
                        for (int j = 0; j < (int)invBinSize + 1; j++)
                        {
                            kojakSparseArray[iTmp][j] = 0;
                        }
                    }
                    jx = (int)((dTmp - iTmp) * invBinSize);
                    if (pfFastXcorrData[i] > 127)
                    {
                        kojakSparseArray[iTmp][jx] = 127;
                    }
                    else if (pfFastXcorrData[i] < -128)
                    {
                        kojakSparseArray[iTmp][jx] = -128;
                    }
                    else if (pfFastXcorrData[i] > 0)
                    {
                        kojakSparseArray[iTmp][jx] = (sbyte)(pfFastXcorrData[i] + 0.5);
                    }
                    else
                    {
                        kojakSparseArray[iTmp][jx] = (sbyte)(pfFastXcorrData[i] - 0.5);
                    }
                }
            }

            var _test_KojakSparse = new List<(int, int, sbyte)>();
            for (int i = 0; i < kojakBins; i++)
            {
                if (kojakSparseArray.Keys.Contains(i))
                {
                    for (int j = 0; j < (int)invBinSize + 1; j++)
                    {
                        if (kojakSparseArray[i][j] != 0)
                        {
                            _test_KojakSparse.Add((i, j, kojakSparseArray[i][j]));
                        }
                    }
                }
            }

            return kojakSparseArray;
        }


        public static double KojakScoring(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, double binSize, double binOffset,
            Dictionary<int, sbyte[]> kojakSparseArray)
        {
            int kojakBins = (int)(scan.TheScan.MassSpectrum.XArray.Last() + 100.0);

            double invBinSize = 1.0 / binSize;
            double mz;

            List<(int, int, double)> keypos = new List<(int, int, double)>();
            for (int i = 0; i < theoreticalProducts.Count; i++)
            {
                for (int j = 1; j < scan.PrecursorCharge; j++)
                {
                    var tp = theoreticalProducts[i];
                    mz = binSize * (int)(tp.NeutralMass.ToMz(j) * invBinSize + binOffset);
                    var key = (int)mz;
                    var pos = (int)((mz - key) * invBinSize);
                    keypos.Add((key, pos, mz));
                }
            }

            int con;
            int match = 0;
            int conFrag = 0;
            double diff;

            double dXcorr = 0.0;
            var _test_keypos = new List<(int, int, double)>();
            for (int i = 0; i < keypos.Count; i++)
            {
                diff = 0;
                con = 0;

                var key = keypos[i].Item1;
                var pos = keypos[i].Item2;
                if (keypos[i].Item3 < 0)
                {
                    mz = binSize * (int)((diff - keypos[i].Item3) * invBinSize + binOffset);
                    key = (int)mz;
                    if (key > kojakBins)
                    {
                        if (con > conFrag)
                        {
                            conFrag = con;
                        }
                        con = 0;
                        break;
                    }

                    if (!kojakSparseArray.Keys.Contains(key))
                    {
                        if (con > conFrag)
                        {
                            conFrag = con;
                        }
                        con = 0;
                        continue;
                    }

                    pos = (int)((mz - key) * invBinSize);
                    dXcorr += kojakSparseArray[key][pos];
                    _test_keypos.Add((key, pos, kojakSparseArray[key][pos]));
                    if (kojakSparseArray[key][pos] > 5)
                    {
                        match++;
                        con++;
                    }
                    else
                    {
                        if (con > conFrag)
                        {
                            conFrag = con;
                        }
                        con = 0;
                    }
                }
                else
                {
                    if (key > kojakBins)
                    {
                        if (con > conFrag)
                        {
                            conFrag = con;
                        }
                        con = 0;
                        break;
                    }

                    if (!kojakSparseArray.Keys.Contains(key))
                    {
                        if (con > conFrag)
                        {
                            conFrag = con;
                        }
                        con = 0;
                        continue;
                    }


                    dXcorr += kojakSparseArray[key][pos];
                    _test_keypos.Add((key, pos, kojakSparseArray[key][pos]));
                    if (kojakSparseArray[key][pos] > 5)
                    {
                        match++;
                        con++;
                    }
                    else
                    {
                        if (con > conFrag)
                        {
                            conFrag = con;
                        }
                        con = 0;
                    }
                }
            }

            if (dXcorr <= 0.0)
            {
                dXcorr = 0.0;
            }
            else
            {
                dXcorr *= 0.005;
            }

            return dXcorr;
        }

        /// <summary>
        /// Similar to the 'KojakScoring' above, the difference is ......I forget now......
        /// Please keep it here untill I figure out which one is more robust.
        /// </summary>
        /// <param name="scan"></param>
        /// <param name="theoreticalProducts"></param>
        /// <param name="binSize"></param>
        /// <param name="binOffset"></param>
        /// <param name="kojakSparseArray"></param>
        /// <returns>Kojak Score</returns>
        public static double KojakScoring2(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, double binSize, double binOffset,
    Dictionary<int, sbyte[]> kojakSparseArray)
        {
            int kojakBins = (int)(scan.TheScan.MassSpectrum.XArray.Last() + 100.0);

            double invBinSize = 1.0 / binSize;
            double mz;

            List<List<List<(int, int, double)>>> xxx = new List<List<List<(int, int, double)>>>();
            for (int charge = 1; charge < scan.PrecursorCharge; charge++)
            {
                List<List<(int, int, double)>> xx = new List<List<(int, int, double)>>();
                foreach (var items in theoreticalProducts.GroupBy(p => p.ProductType))
                {
                    List<(int, int, double)> x = new List<(int, int, double)>();
                    foreach (var item in items)
                    {
                        mz = binSize * (int)(item.NeutralMass.ToMz(charge) * invBinSize + binOffset);
                        var key = (int)mz;
                        var pos = (int)((mz - key) * invBinSize);
                        x.Add((key, pos, mz));
                    }
                    xx.Add(x);
                }
                xxx.Add(xx);
            }


            int con;
            int match = 0;
            int conFrag = 0;
            double diff;

            double dXcorr = 0.0;
            var _test_keypos = new List<(int, int, double)>();

            foreach (var xx in xxx)
            {
                diff = 0;
                foreach (var x in xx)
                {
                    con = 0;
                    foreach (var _x in x)
                    {
                        var key = _x.Item1;
                        var pos = _x.Item2;
                        if (_x.Item3 < 0)
                        {
                            mz = binSize * (int)((diff - _x.Item3) * invBinSize + binOffset);
                            key = (int)mz;
                            if (key > kojakBins)
                            {
                                if (con > conFrag)
                                {
                                    conFrag = con;
                                }
                                con = 0;
                                break;
                            }

                            if (!kojakSparseArray.Keys.Contains(key))
                            {
                                if (con > conFrag)
                                {
                                    conFrag = con;
                                }
                                con = 0;
                                continue;
                            }

                            pos = (int)((mz - key) * invBinSize);
                            dXcorr += kojakSparseArray[key][pos];
                            _test_keypos.Add((key, pos, kojakSparseArray[key][pos]));
                            if (kojakSparseArray[key][pos] > 5)
                            {
                                match++;
                                con++;
                            }
                            else
                            {
                                if (con > conFrag)
                                {
                                    conFrag = con;
                                }
                                con = 0;
                            }
                        }
                        else
                        {
                            if (key > kojakBins)
                            {
                                if (con > conFrag)
                                {
                                    conFrag = con;
                                }
                                con = 0;
                                break;
                            }

                            if (!kojakSparseArray.Keys.Contains(key))
                            {
                                if (con > conFrag)
                                {
                                    conFrag = con;
                                }
                                con = 0;
                                continue;
                            }


                            dXcorr += kojakSparseArray[key][pos];
                            _test_keypos.Add((key, pos, kojakSparseArray[key][pos]));
                            if (kojakSparseArray[key][pos] > 5)
                            {
                                match++;
                                con++;
                            }
                            else
                            {
                                if (con > conFrag)
                                {
                                    conFrag = con;
                                }
                                con = 0;
                            }
                        }
                    }
                }
            }

            if (dXcorr <= 0.0)
            {
                dXcorr = 0.0;
            }
            else
            {
                dXcorr *= 0.005;
            }

            return dXcorr;
        }


        /// <summary>
        /// Xcorr based on 'HiXCorr Bioinfomatics 2015'.
        /// The code is not tested, so it is unknow if it works or not.
        /// The code is left here for possible future usage. 
        /// There is not unit test for it, but it is fine now.
        /// </summary>
        public static double CalcXcorr(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, double fragment_bin_tol)
        {
            int peak_e = scan.NumPeaks;
            int peak_t = theoreticalProducts.Count;

            var experimental_mz = scan.TheScan.MassSpectrum.XArray;
            var experimental_intensity = scan.TheScan.MassSpectrum.YArray;

            var theoretical_mz = theoreticalProducts.Select(p => p.NeutralMass).ToArray();
            // var theoretical_intensity = new List<double>();

            double intenSum = 0;
            double XCorrSum = 0;
            int i = 0;
            while (experimental_mz[i] < 75 * fragment_bin_tol)
            {
                intenSum += experimental_intensity[i];
                i++;
            }

            int left = 0;
            int right = i;
            for (int j = 0; j < peak_t; j++)
            {
                while (right < experimental_mz.Length && experimental_mz[right] < theoretical_mz[j] + (75 * fragment_bin_tol))
                {
                    intenSum += experimental_intensity[right];
                    right++;
                }

                while (left < experimental_mz.Length && experimental_mz[left] < theoretical_mz[j] - (75 * fragment_bin_tol))
                {
                    intenSum -= experimental_intensity[left];
                    left++;
                }
                //XCorrSum += theoretical_intensity[j] * intenSum;
                XCorrSum += intenSum;
            }

            double XCorrAvg = XCorrSum / (2 * 75);
            double XCorr = 0;

            int ix = 0;
            int jx = 0;
            while (ix < peak_e && jx < peak_t)
            {
                if (Math.Abs(experimental_mz[ix] - theoretical_mz[jx]) < fragment_bin_tol)
                {
                    //XCorr += experimental_intensity[ix] * theoretical_intensity[jx];
                    XCorr += experimental_intensity[ix];
                    ix++;
                    jx++;
                }
                else if (experimental_mz[ix] > theoretical_mz[jx])
                {
                    jx++;
                }
                else
                {
                    ix++;
                }
            }
            XCorr -= XCorrAvg;

            return XCorr;
        }

    }
}
