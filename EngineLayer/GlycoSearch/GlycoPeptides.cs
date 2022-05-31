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
    public static class GlycoPeptides
    {
        //Motif should be writen with required form
        public static List<int> GetPossibleModSites(PeptideWithSetModifications peptide, string[] motifs)
        {
            List<int> possibleModSites = new List<int>();

            List<Modification> modifications = new List<Modification>();

            foreach (var mtf in motifs)
            {
                if (ModificationMotif.TryGetMotif(mtf, out ModificationMotif aMotif))
                {
                    Modification modWithMotif = new Modification(_target: aMotif, _locationRestriction: "Anywhere.");
                    modifications.Add(modWithMotif);
                }
            }

            foreach (var modWithMotif in modifications)
            {
                for (int r = 0; r < peptide.Length; r++)
                {
                    if (peptide.AllModsOneIsNterminus.Keys.Contains(r + 2))
                    {
                        continue;
                    }

                    //FullSequence is used here to avoid duplicated modification on same sites?
                    if (ModificationLocalization.ModFits(modWithMotif, peptide.BaseSequence, r + 1, peptide.Length, r + 1))
                    {
                        possibleModSites.Add(r + 2);
                    }
                }
            }

            return possibleModSites;
        }

        public static bool MotifExist(string baseSeq, string[] motifs)
        {
            List<Modification> modifications = new List<Modification>();

            foreach (var mtf in motifs)
            {
                if (ModificationMotif.TryGetMotif(mtf, out ModificationMotif aMotif))
                {
                    Modification modWithMotif = new Modification(_target: aMotif, _locationRestriction: "Anywhere.");
                    modifications.Add(modWithMotif);
                }
            }

            foreach (var modWithMotif in modifications)
            {
                for (int r = 0; r < baseSeq.Length; r++)
                {
                    //Modification is not considered.                  
                    if (ModificationLocalization.ModFits(modWithMotif, baseSeq, r + 1, baseSeq.Length, r + 1))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        public static double[] ScanOxoniumIonFilter(Ms2ScanWithSpecificMass theScan, MassDiffAcceptor massDiffAcceptor, DissociationType dissociationType)
        {
            double[] oxoniumIonsintensities = new double[Glycan.AllOxoniumIons.Length];

            if (dissociationType != DissociationType.HCD && dissociationType != DissociationType.CID && dissociationType != DissociationType.EThcD)
            {
                return oxoniumIonsintensities;
            }

            for (int i = 0; i < Glycan.AllOxoniumIons.Length; i++)
            {
                var oxoMass = ((double)Glycan.AllOxoniumIons[i] / 1E5).ToMass(1);
                var envelope = theScan.GetClosestExperimentalIsotopicEnvelope(oxoMass);
                if (massDiffAcceptor.Accepts(envelope.MonoisotopicMass, oxoMass) >= 0)
                {
                    oxoniumIonsintensities[i] = envelope.TotalIntensity;
                }

            }

            //Normalize by 204. What will happen if 204 is 0.
            if (oxoniumIonsintensities[9] != 0)
            {
                var x204 = oxoniumIonsintensities[9];
                for (int i = 0; i < Glycan.AllOxoniumIons.Length; i++)
                {

                    oxoniumIonsintensities[i] = oxoniumIonsintensities[i] / x204;
                }
            }

            return oxoniumIonsintensities;
        }

        //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons. 
        //Rules are coded in the function.    
        public static bool OxoniumIonsAnalysis(double[] oxoniumIonsintensities, GlycanBox glycanBox)
        {
            //If a glycopeptide spectrum does not have 292.1027 or 274.0921, then remove all glycans that have NeuAc from the search.
            if (oxoniumIonsintensities[10] <= 0 && oxoniumIonsintensities[12] <= 0)
            {
                if (glycanBox.Kind[2] != 0 )
                {
                    return false;
                }
            }

            //If a glycopeptide spectrum does not have 290.0870 or 308.0976, then remove all glycans that have NeuGc from the search.
            if (oxoniumIonsintensities[11] <= 0 && oxoniumIonsintensities[13] <= 0)
            {
                if (glycanBox.Kind[3] != 0)
                {
                    return false;
                }
            }

            //pGlyco3 paper also mentioned the ion 243.0264304315 which is important for phospho-Hex (phoH) exist.

            //If a spectrum has 366.1395, remove glycans that do not have HexNAc(1)Hex(1) or more. Here use the total glycan of glycanBox to calculate. 
            if (oxoniumIonsintensities[14] > 0)
            {
                if (glycanBox.Kind[0] < 1 && glycanBox.Kind[1] < 1)
                {
                    return false;
                }
            }

            //Other rules:
            //A spectrum needs to have 204.0867 to be considered as a glycopeptide.              
            //Ratio of 138.055 to 144.0655 can seperate O/N glycan.

            return true;
        }

        //The function may be dreprecated.
        //TO THINK: it now seems like unnecessary to filter the YCore before match the glycopeptide.
        public static int YCoreIonsFilter(Ms2ScanWithSpecificMass theScan, PeptideWithSetModifications theScanBestPeptide, GlycoType glycoType, MassDiffAcceptor massDiffAcceptor)
        {
            HashSet<double> ycores = new HashSet<double>();

            if(glycoType != GlycoType.OGlycoPep)
            {
                foreach (var gmass in Glycan.NYCoreIons)
                {
                    double mz = theScanBestPeptide.MonoisotopicMass + gmass;
                    ycores.Add(mz);
                }           
            }
            else
            {
                foreach (var gmass in Glycan.OYCoreIons)
                {
                    double mz = theScanBestPeptide.MonoisotopicMass + gmass;
                    ycores.Add(mz);
                }
            }

            int count = 0;
            foreach (var c in ycores)
            {         
                var envelope = theScan.GetClosestExperimentalIsotopicEnvelope(c);
                if (massDiffAcceptor.Accepts(envelope.MonoisotopicMass, c) >= 0)
                {
                    count++;
                }
            }

            return count;
        }

        //Get specific glycan Y ions such 'NN' or 'NH'.
        public static Product GetIndicatorYIon(double peptideMonomassWithNoGlycan, string glycanString)
        {
            Product product = new Product(ProductType.M, FragmentationTerminus.Both, peptideMonomassWithNoGlycan + (double)Glycan.GetMass(glycanString) / 1E5, 0, 0, 0);
            return product;
        }

        public static bool MatchIndicatorYIon(Ms2ScanWithSpecificMass scan, Product theoreticalProduct, CommonParameters commonParameters)
        {
            List<Product> products = new List<Product>();
            products.Add(theoreticalProduct);
            var x = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters).Count();

            foreach (var childScan in scan.ChildScans)
            {
                x += MetaMorpheusEngine.MatchFragmentIons(childScan, products, commonParameters).Count();
            }

            return x > 0;
        }

        //Find Glycan index or glycanBox index.
        public static int BinarySearchGetIndex(double[] massArray, double targetMass)
        {
            var iD = Array.BinarySearch(massArray, targetMass);
            if (iD < 0) { iD = ~iD; }
            else
            {
                while (iD - 1 >= 0 && massArray[iD - 1] >= targetMass - 0.00000001)
                {
                    iD--;
                }
            }
            return iD;
        }

        public static bool DissociationTypeContainETD(DissociationType dissociationType)
        {
            if (dissociationType == DissociationType.ETD || dissociationType == DissociationType.EThcD)
            {
                return true;
            }

            return false;
        }

        public static bool DissociationTypeContainHCD(DissociationType dissociationType)
        {
            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID || dissociationType == DissociationType.EThcD)
            {
                return true;
            }

            return false;
        }

        public static double CalculateSinglePeptideScore(List<MatchedFragmentIon> matchedFragmentIons, List<Product> theoreticalProducts, CommonParameters commonParameters, double g = 0.70, double w = 0.35)
        {
            //The function is using pGlyco score equation.
            double score = 0;
            double score_p = 0;

            double match_peps_count = 0;
            var theo_peps_count = theoreticalProducts.Count(p => p.ProductType != ProductType.M && p.ProductType != ProductType.D);

            foreach (var m in matchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M && p.NeutralTheoreticalProduct.ProductType != ProductType.D))
            {
                score_p += Math.Log10(m.Intensity) * (1 - Math.Pow(m.MassErrorPpm / commonParameters.ProductMassTolerance.Value, 4));
                match_peps_count++;
            }

            score = (1 - w)*score_p * Math.Pow(match_peps_count / theo_peps_count, g);

            return score;
        }


        /// <summary>
        /// pGlyco2 score function.
        /// Ref 'pGlyco 2.0 enables precision N-glycoproteomics with comprehensive quality control and one-step mass spectrometry for intact glycopeptide identification'
        /// </summary>
        public static double CalculatePeptideScore(List<MatchedFragmentIon> matchedFragmentIons, List<Product> theoreticalProducts, CommonParameters commonParameters, double a = 0.56, double b = 1.42, double g = 0.70, double w = 0.35)
        {
            //The function is using pGlyco score equation.
            double score = 0;
            double score_g = 0;
            double score_p = 0;

            double match_y_count = 0;
            var theo_y_count = theoreticalProducts.Count(p => p.ProductType == ProductType.Y);

            double match_ycores_count = 0;
            var theo_ycores_count = theoreticalProducts.Count(p => p.ProductType == ProductType.Ycore );

            double match_peps_count = 0;
            var theo_peps_count = theoreticalProducts.Count(p => p.ProductType != ProductType.Ycore && p.ProductType != ProductType.Y && p.ProductType != ProductType.D);

            foreach (var m in matchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.Y))
            {
                score_g += Math.Log(m.Intensity) * (1 - Math.Pow(m.MassErrorPpm / commonParameters.ProductMassTolerance.Value, 4));
                match_y_count++;
            }

            foreach (var m in matchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.Ycore))
            {
                score_g += Math.Log(m.Intensity) * (1 - Math.Pow(m.MassErrorPpm / commonParameters.ProductMassTolerance.Value, 4));
                match_ycores_count++;
            }

            foreach(var m in matchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType != ProductType.Ycore
                && p.NeutralTheoreticalProduct.ProductType != ProductType.Y
                && p.NeutralTheoreticalProduct.ProductType != ProductType.D))
            {
                score_p += Math.Log(m.Intensity) * (1 - Math.Pow(m.MassErrorPpm / commonParameters.ProductMassTolerance.Value, 4));
                match_peps_count++;
            }


            if (theo_y_count == 0)
            {
                score = w * Math.Pow(match_ycores_count / theo_ycores_count, b)
                    + (1 - w) * score_p * Math.Pow(match_peps_count / theo_peps_count, g);
            }
            else
            {
                score = w * score_g * Math.Pow(match_y_count / theo_y_count, a)
                    * Math.Pow(match_ycores_count / theo_ycores_count, b)
                    + (1 - w) * score_p * Math.Pow(match_peps_count / theo_peps_count, g);
            }
            return score;
        }

        /// <summary>
        /// The function is for XCorr high res data, based on 'Kojak Crosslink'.
        /// </summary>
        public static void BinIons(Ms2ScanWithSpecificMass scan, double binSize, double invBinSize, double binOffset,
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

        public static void MakeCorrData(double scale, ref int iHighestIon, ref double dHighestIntensity, ref float[] pdCorrelationData_YArray, ref int xCorrArraySize, ref double[] pdTempRawData)
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

                if (x<iNumWindows)
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
            for (int i = 75; i < xCorrArraySize+75; i++)
            {
                if (i<xCorrArraySize && pdCorrelationData_YArray[i]>0)
                {
                    dSum += pdCorrelationData_YArray[i];
                }
                if (i>= 151 && pdCorrelationData_YArray[i-151] >0)
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
            for (int i = 2; i < xCorrArraySize-1; i++)
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
                        for (int j = 0; j < (int)invBinSize+1; j++)
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
                        if (kojakSparseArray[i][j]!=0)
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
                if (keypos[i].Item3 <0)
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
            while (experimental_mz[i] < 75* fragment_bin_tol)
            {
                intenSum += experimental_intensity[i];
                i++;
            }

            int left = 0;
            int right = i;
            for (int j = 0; j < peak_t; j++)
            {
                while ( right < experimental_mz.Length && experimental_mz[right] < theoretical_mz[j] + (75*fragment_bin_tol))
                {
                    intenSum += experimental_intensity[right];
                    right++;
                }

                while (left < experimental_mz.Length && experimental_mz[left] < theoretical_mz[j] - (75*fragment_bin_tol))
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
                if ( Math.Abs(experimental_mz[ix] - theoretical_mz[jx])< fragment_bin_tol)
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

        #region Glycopeptide and fragmentation functions

        //TO THINK: if one peptide modified with >1 nglycans, how Yions will be generated. Here we just simplified to treat it as two seperate generation.
        public static List<Product> GetGlycanYIons(PeptideWithSetModifications peptide, Glycan[] glycans)
        {
            double possiblePeptideMass = peptide.MonoisotopicMass;
            List<Product> YIons = new List<Product>();
            HashSet<int> masses = new HashSet<int>();

            foreach (var glycan in glycans)
            {
                foreach (var ion in glycan.Ions)
                {
                    if (masses.Contains(ion.IonMass))
                    {
                        continue;
                    }
                    masses.Add(ion.IonMass);
                    if ( Glycan.NYCoreIntIons.Contains(ion.IonMass))
                    {
                        Product product = new Product(ProductType.Ycore, FragmentationTerminus.Both, possiblePeptideMass + (double)ion.IonMass / 1E5, 0, 0, (double)ion.IonMass / 1E5);
                        YIons.Add(product);                    
                    }
                    else
                    {
                        Product product = new Product(ProductType.Y, FragmentationTerminus.Both, possiblePeptideMass + (double)ion.IonMass / 1E5, 0, 0, (double)ion.IonMass / 1E5);
                        YIons.Add(product);
                    }

                }
            }
            return YIons;
        }

        public static List<Product> GetOxoniumIons(Glycan[] glycans)
        {
            List<Product> OxoIons = new List<Product>();

            HashSet<int> masses = new HashSet<int>();

            foreach (var g in glycans)
            {
                foreach (var d in g.DiagnosticIons)
                {
                    if (masses.Contains(d))
                    {
                        continue;
                    }
                    masses.Add(d);

                    var dmass = (double)d / 1E5;
                    Product D = new Product(ProductType.D, FragmentationTerminus.None, dmass, (int)Math.Round(dmass.ToMz(1), 0), 0, 0);
                    OxoIons.Add(D);
                }
            }
            return OxoIons;
        }

        public static PeptideWithSetModifications GlyGetTheoreticalPeptide(int[] modPos, PeptideWithSetModifications peptide, Modification[] modifications)
        {
            if (modPos.Length == 0)
            {
                return null;
            }

            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            for (int i = 0; i < modPos.Length; i++)
            {
                testMods.Add(modPos[i], modifications[i]);
            }

            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

        public static PeptideWithSetModifications GlyGetTheoreticalPeptide(int[] modPos, PeptideWithSetModifications peptide, GlycanBox glycanBox, Modification[] globalMods)
        {
            Modification[] modifications = new Modification[glycanBox.ModCount];
            for (int i = 0; i < glycanBox.ModCount; i++)
            {
                modifications[i] = globalMods[glycanBox.ModIds.ElementAt(i)];
            }

            return GlyGetTheoreticalPeptide(modPos, peptide, modifications);
        }

        public static PeptideWithSetModifications GlyGetTheoreticalPeptide(Route theModPositions, PeptideWithSetModifications peptide, Modification[] globalMods)
        {
            Modification[] modifications = new Modification[theModPositions.Mods.Count];
            for (int i = 0; i < theModPositions.Mods.Count; i++)
            {
                modifications[i] = globalMods[theModPositions.Mods[i].GlycanID];
            }

            return GlyGetTheoreticalPeptide(theModPositions.Mods.Select(p=>p.ModSite).ToArray(), peptide, modifications);
        }

        public static List<Product> GlyGetTheoreticalFragments(GlycoType glycanType, DissociationType dissociationType, PeptideWithSetModifications peptide, PeptideWithSetModifications modPeptide, List<int> NPos, Glycan[] glycans)
        {
            List<Product> theoreticalProducts = new List<Product>();
            HashSet<double> masses = new HashSet<double>();
            List<Product> products = new List<Product>();

            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID)
            {
                peptide.Fragment(dissociationType, FragmentationTerminus.Both, products);

                if (glycanType == GlycoType.NGlycoPep || glycanType == GlycoType.MixedGlycoPep)
                {
                    List<Product> shiftProducts = new List<Product>();

                    foreach (var pd in products)
                    {
                        if (pd.ProductType == ProductType.b)
                        {
                            int count = NPos.Count(p => pd.AminoAcidPosition >= p - 1 );
                            if (count > 0)
                            {
                                Product b = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(b);
                            }
                        }
                        else if (pd.ProductType == ProductType.y)
                        {
                            int count = NPos.Count(p => pd.AminoAcidPosition < p );
                            if (count > 0)
                            {
                                Product y = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(y);
                            }
                        }
                    }

                    products.AddRange(shiftProducts);
                }

                products.AddRange(GetOxoniumIons(glycans));

                products.AddRange(GetGlycanYIons(peptide, glycans));

            }
            else if (dissociationType == DissociationType.ETD)
            {
                modPeptide.Fragment(dissociationType, FragmentationTerminus.Both, products);
            }
            else if (dissociationType == DissociationType.EThcD)
            {
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);

                if (glycanType == GlycoType.NGlycoPep || glycanType == GlycoType.MixedGlycoPep)
                {
                    List<Product> shiftProducts = new List<Product>();

                    foreach (var pd in products)
                    {
                        if (pd.ProductType == ProductType.b)
                        {
                            int count = NPos.Count(p => p > pd.AminoAcidPosition);
                            if (count > 0)
                            {
                                Product b = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(b);
                            }
                        }
                        else if (pd.ProductType == ProductType.y)
                        {
                            int count = NPos.Count(p => p < pd.AminoAcidPosition);
                            if (count > 0)
                            {
                                Product y = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(y);
                            }
                        }
                    }

                    products.AddRange(shiftProducts);
                }

                List<Product> etdProducts = new List<Product>();

                modPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, etdProducts);

                products.AddRange(etdProducts.Where(p => p.ProductType != ProductType.y));

                products.AddRange(GetOxoniumIons(glycans));

                products.AddRange(GetGlycanYIons(peptide, glycans));

            }          

            foreach (var fragment in products)
            {
                if (!masses.Contains(fragment.NeutralMass))
                {
                    masses.Add(fragment.NeutralMass);
                    theoreticalProducts.Add(fragment);
                }
            }

            return theoreticalProducts;
        }


        public static List<Product> GlyGetPepHCDFragments(GlycoType glycanType, DissociationType dissociationType, PeptideWithSetModifications peptide, List<int> NPos, Glycan[] glycans)
        {
            List<Product> products = new List<Product>();

            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID || dissociationType == DissociationType.EThcD)
            {
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);

                if (glycanType == GlycoType.NGlycoPep || glycanType == GlycoType.MixedGlycoPep)
                {
                    List<Product> shiftProducts = new List<Product>();

                    foreach (var pd in products)
                    {
                        if (pd.ProductType == ProductType.b)
                        {
                            int count = NPos.Count(p => pd.AminoAcidPosition >= p - 1);
                            if (count > 0)
                            {
                                Product b = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, Glycan.HexNAcMass);
                                shiftProducts.Add(b);
                            }
                        }
                        else if (pd.ProductType == ProductType.y)
                        {
                            int count = NPos.Count(p => pd.AminoAcidPosition < p);
                            if (count > 0)
                            {
                                Product y = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass +  Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, Glycan.HexNAcMass);
                                shiftProducts.Add(y);
                            }
                        }
                    }

                    products.AddRange(shiftProducts);
                }

                products.AddRange(GetGlycanYIons(peptide, glycans));

            }

            //TO THINK: How necessary to remove the products with shared mass?
            //List<Product> theoreticalProducts = new List<Product>();
            //HashSet<double> masses = new HashSet<double>();
            //foreach (var fragment in products)
            //{
            //    if (!masses.Contains(fragment.NeutralMass))
            //    {
            //        masses.Add(fragment.NeutralMass);
            //        theoreticalProducts.Add(fragment);
            //    }
            //}

            return products;
        }
        //Find FragmentMass for the fragments that contain localization Information.
        public static List<double> GetLocalFragmentGlycan(List<Product> products, int modInd, int childBoxInd, LocalizationGraph localizationGraph)
        {
            List<double> newFragments = new List<double>();

            var local_c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] - 1 && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1] - 1).ToList();
            foreach (var c in local_c_fragments)
            {
                var newMass = c.NeutralMass + localizationGraph.ChildModBoxes[childBoxInd].Mass;
                newFragments.Add(newMass);
            }

            var local_z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1]).ToList();

            foreach (var z in local_z_fragments)
            {
                var newMass = z.NeutralMass + (localizationGraph.ModBox.Mass - localizationGraph.ChildModBoxes[childBoxInd].Mass);
                newFragments.Add(newMass);
            }

            var local_b_fragments = products.Where(p => p.ProductType == ProductType.b && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] - 1 && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1] - 1).ToList();
            foreach (var b in local_b_fragments)
            {
                newFragments.Add(b.NeutralMass);
                var newMass = b.NeutralMass + ((GlycanBox)localizationGraph.ChildModBoxes[childBoxInd]).NGlycanCount * Glycan.HexNAcMass;
                newFragments.Add(newMass);
            }

            var local_y_fragments = products.Where(p => p.ProductType == ProductType.y && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1]).ToList();
            foreach (var y in local_y_fragments)
            {
                newFragments.Add(y.NeutralMass);
                var newMass = y.NeutralMass + (((GlycanBox)localizationGraph.ModBox).NGlycanCount - ((GlycanBox)localizationGraph.ChildModBoxes[childBoxInd]).NGlycanCount) * Glycan.HexNAcMass;
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        //Find FragmentMass for the fragments that doesn't contain localization Information. For example, "A|TAABBS|B", c1 and c7, z1 and z7, z8 ion don't contain localization information.
        public static List<double> GetUnlocalFragmentGlycan(List<Product> products, int[] modPoses, ModBox modBox)
        {
            var mass = modBox.Mass;

            List<double> newFragments = new List<double>();

            var c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition < modPoses.First() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(c_fragments);

            var c_fragments_shift = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= modPoses.Last() - 1).Select(p => p.NeutralMass);

            foreach (var c in c_fragments_shift)
            {
                var newMass = c + mass;
                newFragments.Add(newMass);
            }

            var z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition > modPoses.Last() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(z_fragments);

            var z_fragments_shift = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition <= modPoses.First() - 1).Select(p => p.NeutralMass);

            foreach (var z in z_fragments_shift)
            {
                var newMass = z + mass;
                newFragments.Add(newMass);
            }

            var b_fragments = products.Where(p => p.ProductType == ProductType.b && p.AminoAcidPosition < modPoses.First() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(b_fragments);

            var b_fragments_shift = products.Where(p => p.ProductType == ProductType.b && p.AminoAcidPosition >= modPoses.Last() - 1).Select(p => p.NeutralMass);
            double nlocalMass = ((GlycanBox)modBox).NGlycanCount * Glycan.HexNAcMass;
            foreach (var b in b_fragments_shift)
            {
                newFragments.Add(b);
                var newMass = b + nlocalMass;
                newFragments.Add(newMass);
            }

            var y_fragments = products.Where(p => p.ProductType == ProductType.y && p.AminoAcidPosition > modPoses.Last() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(y_fragments);

            var y_fragments_shift = products.Where(p => p.ProductType == ProductType.y && p.AminoAcidPosition <= modPoses.First() - 1).Select(p => p.NeutralMass);

            foreach (var y in y_fragments_shift)
            {
                newFragments.Add(y);
                var newMass = y + nlocalMass;
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        //According to pGlyco3, under HCD N-glycopeptide contain >=2 core ions. O-Glycopeptide contains >=1 core ions. 
        //Such a filter can remove a lot of unnecessary searches.
        public static bool ScanTrimannosylCoreFilter(List<MatchedFragmentIon> matchedFragmentIons, GlycoType glycoType)
        {
            if (glycoType == GlycoType.OGlycoPep && matchedFragmentIons.Count(p=>p.Annotation.Contains("Ycore")) >= 1)
            {
                return true;
            }
            else if ((glycoType == GlycoType.NGlycoPep || glycoType == GlycoType.MixedGlycoPep) 
                && matchedFragmentIons.Count(p => p.Annotation.Contains("Ycore")) >= 2)
            {
                return true;
            }
            return false;
        }

        #endregion


        #region O-Glyco related functions, plan to be deprecated.

        //TO THINK: filter reasonable fragments here. The final solution is to change mzLib.Proteomics.PeptideWithSetModifications.Fragment
        public static List<Product> OGlyGetTheoreticalFragments(DissociationType dissociationType, PeptideWithSetModifications peptide, PeptideWithSetModifications modPeptide)
        {
            List<Product> theoreticalProducts = new List<Product>();        
            HashSet<double> masses = new HashSet<double>();

            List<Product> products = new List<Product>();
            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID)
            {
                List<Product> diag = new List<Product>();
                modPeptide.Fragment(dissociationType, FragmentationTerminus.Both, diag);
                peptide.Fragment(dissociationType, FragmentationTerminus.Both, products);
                products = products.Concat(diag.Where(p => p.ProductType != ProductType.b && p.ProductType != ProductType.y)).ToList();
            }
            else if(dissociationType == DissociationType.ETD)
            {
                modPeptide.Fragment(dissociationType, FragmentationTerminus.Both, products);
            }
            else if(dissociationType == DissociationType.EThcD)
            {
                List<Product> diag = new List<Product>();
                modPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, diag);
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
                products = products.Concat(diag.Where(p => p.ProductType != ProductType.b && p.ProductType != ProductType.y)).ToList();

                List<Product> etdProducts = new List<Product>();
                modPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, etdProducts);
                products = products.Concat(etdProducts.Where(p => p.ProductType != ProductType.y)).ToList();
            }

            foreach (var fragment in products)
            {
                if (!masses.Contains(fragment.NeutralMass))
                {
                    masses.Add(fragment.NeutralMass);
                    theoreticalProducts.Add(fragment);
                }           
            }

            return theoreticalProducts;
        }

        //Find FragmentMass for the fragments that contain localization Information.
        public static List<double> GetLocalFragment(List<Product> products, int modInd, int childBoxInd, LocalizationGraph localizationGraph)
        {
            List<double> newFragments = new List<double>();
            var local_c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] - 1 && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1] - 1).ToList();

            foreach (var c in local_c_fragments)
            {
                var newMass = c.NeutralMass + localizationGraph.ChildModBoxes[childBoxInd].Mass;
                newFragments.Add(newMass);
            }

            var local_z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1]).ToList();

            foreach (var z in local_z_fragments)
            {
                var newMass = z.NeutralMass + (localizationGraph.ModBox.Mass - localizationGraph.ChildModBoxes[childBoxInd].Mass);
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        //Find FragmentMass for the fragments that doesn't contain localization Information. For example, "A|TAABBS|B", c1 and c7, z1 and z7, z8 ion don't contain localization information.
        public static List<double> GetUnlocalFragment(List<Product> products, int[] modPoses, ModBox OGlycanBox)
        {
            var mass = OGlycanBox.Mass;

            List<double> newFragments = new List<double>();
            var c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition < modPoses.First() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(c_fragments);

            var c_fragments_shift = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= modPoses.Last() - 1).Select(p => p.NeutralMass);

            foreach (var c in c_fragments_shift)
            {
                var newMass = c + mass;
                newFragments.Add(newMass);
            }

            var z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition > modPoses.Last() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(z_fragments);

            var z_fragments_shift = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition <= modPoses.First() - 1).Select(p => p.NeutralMass);

            foreach (var z in z_fragments_shift)
            {
                var newMass = z + mass;
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        #endregion


        #region Functions are not used now, could be useful in the future. 

        //NGlycopeptide usually contain Y ions with different charge states, especially in sceHCD data. 
        //The purpose of this function is to try match all Y ion with different charges. The usage of this function requires further investigation. 
        //Not sure about OGlycopeptide. 
        public static List<MatchedFragmentIon> GlyMatchOriginFragmentIons(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();

            // if the spectrum has no peaks
            if (scan.ExperimentalFragments != null && !scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum

            for (int id = 0; id < theoreticalProducts.Count; id++)
            {
                var product = theoreticalProducts[id];
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                if (product.ProductType == ProductType.Y || product.ProductType == ProductType.Ycore)
                {
                    for (int i = 1; i <= scan.PrecursorCharge; i++)
                    {

                        var closestExperimentalMz = scan.GetClosestExperimentalFragmentMz(product.NeutralMass.ToMz(i), out double? intensity);

                        if (closestExperimentalMz.HasValue && commonParameters.ProductMassTolerance.Within(closestExperimentalMz.Value, product.NeutralMass.ToMz(i)))
                        {
                            matchedFragmentIons.Add(new MatchedFragmentIon(ref product, closestExperimentalMz.Value, intensity.Value, i));
                        }
                    }
                }

                else
                {
                    // get the closest peak in the spectrum to the theoretical peak
                    var closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(product.NeutralMass);

                    // is the mass error acceptable?
                    if (closestExperimentalMass != null && commonParameters.ProductMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, product.NeutralMass) && closestExperimentalMass.Charge <= scan.PrecursorCharge)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(ref product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                            closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge));
                    }
                }
            }

            return matchedFragmentIons;
        }


        //The function here is to calculate permutation localization which could be used to compare with Graph-Localization.
        public static List<int[]> GetPermutations(List<int> allModPos, int[] glycanBoxId)
        {
            var length = glycanBoxId.Length;
            var indexes = Enumerable.Range(0, length).ToArray();
            int[] orderGlycan = new int[length];

            List<int[]> permutateModPositions = new List<int[]>();

            var combinations = Glycan.GetKCombs(allModPos, length);

            foreach (var com in combinations)
            {
                var permutation = Glycan.GetPermutations(com, length);

                HashSet<string> keys = new HashSet<string>();

                foreach (var per in permutation)
                {
                    Array.Sort(indexes);

                    var orderedPer = per.ToArray();
                    Array.Sort(orderedPer, indexes);

                    for (int i = 0; i < length; i++)
                    {
                        orderGlycan[i] = glycanBoxId[indexes[i]];
                    }
                    var key = string.Join(",", orderGlycan.Select(p => p.ToString()));
                    if (!keys.Contains(key))
                    {
                        keys.Add(key);
                        permutateModPositions.Add(per.ToArray());
                    }
                }
            }

            return permutateModPositions;
        }

        #endregion

    }
}
