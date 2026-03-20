using Chemistry;
using MassSpectrometry;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using Transcriptomics.Digestion;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {
        protected static readonly Dictionary<DissociationType, List<double>> complementaryIonConversionDictionary = new Dictionary<DissociationType, List<double>>
        {
            { DissociationType.LowCID, new List<double>(){ Constants.ProtonMass } },
            { DissociationType.HCD, new List<double>(){ Constants.ProtonMass } },
            { DissociationType.ETD,new List<double>() {2 * Constants.ProtonMass } }, //presence of zplusone (zdot) makes this two instead of one
            { DissociationType.CID,new List<double>() {Constants.ProtonMass } },
            { DissociationType.EThcD,new List<double>() {Constants.ProtonMass, 2 * Constants.ProtonMass } },
            //TODO: refactor such that complementary ions are generated specifically for their complementary pair.
            //TODO: create a method to auto-determine the conversion
        };

        public readonly CommonParameters CommonParameters;
        protected readonly List<(string FileName, CommonParameters Parameters)> FileSpecificParameters;
        protected readonly List<string> NestedIds;

        protected MetaMorpheusEngine(CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds)
        {
            CommonParameters = commonParameters;
            FileSpecificParameters = fileSpecificParameters;
            NestedIds = nestedIds;
        }

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        public static double CalculatePeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons, bool fragmentsCanHaveDifferentCharges = false)
        {
            if(fragmentsCanHaveDifferentCharges)
            {
                return CalculateAllChargesPeptideScore(thisScan, matchedFragmentIons);
            }

            double score = 0;

            if (thisScan.MassSpectrum.XcorrProcessed)
            {
                // XCorr
                foreach (var fragment in matchedFragmentIons)
                {
                    switch (fragment.NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.aDegree:
                        case ProductType.aStar:
                        case ProductType.bWaterLoss:
                        case ProductType.bAmmoniaLoss:
                        case ProductType.yWaterLoss:
                        case ProductType.yAmmoniaLoss:
                            score += 0.01 * fragment.Intensity;
                            break;
                        case ProductType.D: //count nothing for diagnostic ions.
                            break;
                        default:
                            score += 1 * fragment.Intensity;
                            break;
                    }
                }
            }
            else
            {
                // Morpheus score
                for (int i = 0; i < matchedFragmentIons.Count; i++)
                {
                    switch (matchedFragmentIons[i].NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.D:
                            break;
                        default:
                            score += 1 + matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent;
                            break;
                    }
                    
                }
            }

            return score;
        }

        //Used only when user wants to generate spectral library.
        //Normal search only looks for one match ion for one fragment, and if it accepts it then it doesn't try to look for different charge states of that same fragment. 
        //The score will be the number of matched ions and plus some fraction calculated by intensity(matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent).
        //Like b1, b2, b3 will have score 3.xxx;But when generating library, we need look for match ions with all charges.So we will have b1,b2,b3, b1^2, b2^3. If using 
        //the normal scoring function, the score will be 5.xxxx, which is not proper. The score for b1 and b1^2 should also be 1 plus some some fraction calculated by intensity, 
        //because they are matching the same fragment ion just with different charges. So b1, b2, b3, b1^2, b2^3 should be also 3.xxx(but a little higher than b1, b2, b3 as 
        //the fraction part) rather than 5.xxx.
        private static double CalculateAllChargesPeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons)
        {
            double score = 0;

            // Morpheus score
            List<String> ions = new List<String>();
            for (int i = 0; i < matchedFragmentIons.Count; i++)
            {
                String ion = $"{ matchedFragmentIons[i].NeutralTheoreticalProduct.ProductType.ToString()}{  matchedFragmentIons[i].NeutralTheoreticalProduct.FragmentNumber}";
                if (ions.Contains(ion))
                {
                    score += matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent;
                }
                else
                {
                    score += 1 + matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent;
                    ions.Add(ion);
                }
            }

            return score;

        }

        public static List<MatchedFragmentIon> MatchFragmentIons(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters, bool matchAllCharges = false)
        {
            if (matchAllCharges)
            {
                return MatchFragmentIonsOfAllCharges(scan, theoreticalProducts, commonParameters);
            }

            var matchedFragmentIons = new List<MatchedFragmentIon>();

            if (scan.TheScan.MassSpectrum.XcorrProcessed && scan.TheScan.MassSpectrum.XArray.Length != 0)
            {

                for (int i = 0; i < theoreticalProducts.Count; i++)
                {
                    var product = theoreticalProducts[i];
                    // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                    if (double.IsNaN(product.NeutralMass))
                    {
                        continue;
                    }

                    // Magic number represents mzbinning space. 
                    double theoreticalFragmentMz = Math.Round(product.NeutralMass.ToMz(1) / 1.0005079, 0) * 1.0005079;
                    var closestMzIndex = scan.TheScan.MassSpectrum.GetClosestPeakIndex(theoreticalFragmentMz);

                    if (commonParameters.ProductMassTolerance.Within(scan.TheScan.MassSpectrum.XArray[closestMzIndex], theoreticalFragmentMz))
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, theoreticalFragmentMz, scan.TheScan.MassSpectrum.YArray[closestMzIndex], 1));
                    }
                }

                return matchedFragmentIons;
            }

            // if the spectrum has no peaks
            if (scan.ExperimentalFragments != null && !scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum
            for (int i = 0; i < theoreticalProducts.Count; i++)
            {
                var product = theoreticalProducts[i];
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // get the closest peak in the spectrum to the theoretical peak
                var closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(product.NeutralMass);

                // is the mass error acceptable?
                if (closestExperimentalMass != null
                    && commonParameters.ProductMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, product.NeutralMass)
                    && Math.Abs(closestExperimentalMass.Charge) <= Math.Abs(scan.PrecursorCharge))//TODO apply this filter before picking the envelope
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                        closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge));
                }
            }
            if (commonParameters.AddCompIons)
            {
                foreach (double massShift in complementaryIonConversionDictionary[commonParameters.DissociationType])
                {
                    double protonMassShift = massShift.ToMass(1);

                    for (int i = 0; i < theoreticalProducts.Count; i++)
                    {
                        var product = theoreticalProducts[i];
                        // unknown fragment mass or diagnostic ion or precursor; skip those
                        if (double.IsNaN(product.NeutralMass) || product.ProductType == ProductType.D || product.ProductType == ProductType.M)
                        {
                            continue;
                        }

                        double compIonMass = scan.PrecursorMass + protonMassShift - product.NeutralMass;

                        // get the closest peak in the spectrum to the theoretical peak
                        IsotopicEnvelope closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(compIonMass);

                        // is the mass error acceptable?
                        if (commonParameters.ProductMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, compIonMass) && closestExperimentalMass.Charge <= scan.PrecursorCharge)
                        {
                            //found the peak, but we don't want to save that m/z because it's the complementary of the observed ion that we "added". Need to create a fake ion instead.
                            double mz = (scan.PrecursorMass + protonMassShift - closestExperimentalMass.MonoisotopicMass).ToMz(closestExperimentalMass.Charge);

                            matchedFragmentIons.Add(new MatchedFragmentIon(product, mz, closestExperimentalMass.TotalIntensity, closestExperimentalMass.Charge));
                        }
                    }
                }
            }

            return matchedFragmentIons;
        }
        
        //Used only when user wants to generate spectral library.
        //Normal search only looks for one match ion for one fragment, and if it accepts it then it doesn't try to look for different charge states of that same fragment. 
        //But for library generation, we need find all the matched peaks with all the different charges.
        private static List<MatchedFragmentIon> MatchFragmentIonsOfAllCharges(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();
            var ions = new List<string>();

            // if the spectrum has no peaks
            if (scan.ExperimentalFragments != null && !scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum
            foreach (Product product in theoreticalProducts)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                //get the range we can accept 
                var minMass = commonParameters.ProductMassTolerance.GetMinimumValue(product.NeutralMass);
                var maxMass = commonParameters.ProductMassTolerance.GetMaximumValue(product.NeutralMass);
                var closestExperimentalMassList = scan.GetClosestExperimentalIsotopicEnvelopeList(minMass, maxMass);
                if (closestExperimentalMassList != null)
                {
                    foreach (var x in closestExperimentalMassList)
                    {
                        String ion = $"{product.ProductType.ToString()}{ product.FragmentNumber}^{x.Charge}-{product.NeutralLoss}";
                        if (x != null 
                            && !ions.Contains(ion) 
                            && commonParameters.ProductMassTolerance.Within(x.MonoisotopicMass, product.NeutralMass) 
                            && Math.Abs(x.Charge) <= Math.Abs(scan.PrecursorCharge))//TODO apply this filter before picking the envelope
                        {
                            Product temProduct = product;
                            matchedFragmentIons.Add(new MatchedFragmentIon(temProduct, x.MonoisotopicMass.ToMz(x.Charge),
                                x.Peaks.First().intensity, x.Charge));

                            ions.Add(ion);
                        }
                    }
                }
            }

            return matchedFragmentIons;
        }
        protected abstract MetaMorpheusEngineResults RunSpecific();

        public MetaMorpheusEngineResults Run()
        {
            DetermineAnalyteType(CommonParameters);
            StartingSingleEngine();
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            this.CommonParameters.SetCustomProductTypes();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            FinishedSingleEngine(myResults);
            return myResults;
        }

        public Task<MetaMorpheusEngineResults> RunAsync() => Task.Run(Run);

        /// <summary>
        /// Determines and sets the analyte type based on CommonParameters digestion settings.
        /// This method is called automatically by MetaMorpheusEngine.Run() to handle:
        /// - RNA mode (RnaDigestionParams → Oligo)
        /// - Top-down mode (protease == "top-down" → Proteoform)  
        /// - Bottom-up/default mode (→ Peptide)
        /// 
        /// IMPORTANT: This recalculates the analyte type at runtime and may differ from the GUI mode
        /// set by GuiGlobalParamsViewModel.IsRnaMode. This is intentional to support:
        /// - File-specific parameters with different modes
        /// - Mixed mode workflows
        /// 
        /// For GUI initialization, rely on GuiGlobalParamsViewModel.IsRnaMode which sets 
        /// GlobalVariables.AnalyteType. Do NOT call this method during GUI task window initialization.
        /// </summary>
        /// <param name="commonParameters"></param>
        public static void DetermineAnalyteType(CommonParameters commonParameters)
        {
            // Comment made while DetermineAnalyteType happened at the task layer
            // TODO: note that this will not function well if the user is using file-specific settings, but it's assumed
            // that bottom-up and top-down data is not being searched in the same task. 

            // Update: Now that it is in the engine layer, analyte type specific operations will be okay at the engine layer, meaning searching top-down and bottom-up with file specific params will execute the proper control flow. However, a problem still exists in PostSearchAnalysis where that analyte type will be set to whatever the main parameters are. 

            if (commonParameters == null || commonParameters.DigestionParams == null)
                return;

            GlobalVariables.AnalyteType = commonParameters.DetermineAnalyteType();
        }

        #region Event Helpers

        public string GetId()
        {
            return string.Join(",", NestedIds);
        }

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, NestedIds));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, NestedIds));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        private void StartingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void FinishedSingleEngine(MetaMorpheusEngineResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }

        #endregion
    }
}
