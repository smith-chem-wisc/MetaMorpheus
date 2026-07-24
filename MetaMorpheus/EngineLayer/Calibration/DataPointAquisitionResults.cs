using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.Statistics;
using System.IO;

namespace EngineLayer.Calibration
{
    /// <summary>
    /// Returns PSMs that can be used for calibration based on tolerance limits passed from the calibration task
    /// </summary>
    public class DataPointAquisitionResults : MetaMorpheusEngineResults
    {
        public DataPointAquisitionResults(
            MetaMorpheusEngine dataPointAcquisitionEngine,
            List<SpectralMatch> psms,
            List<LabeledDataPoint> ms1List,
            List<LabeledDataPoint> ms2List,
            int numMs1MassChargeCombinationsConsidered,
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
            int numMs2MassChargeCombinationsConsidered,
            int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
            CommonParameters commonParameters = null)
            : base(dataPointAcquisitionEngine)
        {
            Psms = psms;

            Ms1List = ms1List;
            Ms2List = ms2List;

            var ms1Range = Ms1List.Select(b => b.ExperimentalMz - b.TheoreticalMz).ToArray();
            var ms2Range = Ms2List.Select(b => b.ExperimentalMz - b.TheoreticalMz).ToArray();

            var ms1PpmRange = Ms1List.Select(b => (b.ExperimentalMz - b.TheoreticalMz) / b.TheoreticalMz).ToArray();
            var ms2PpmRange = Ms2List.Select(b => (b.ExperimentalMz - b.TheoreticalMz) / b.TheoreticalMz).ToArray();

            // Calibration must see instrument m/z drift, not isotope-assignment error. A most-abundant
            // search admits PSMs whose deconvoluted monoisotopic peak is off by whole isotopologues (the
            // apex notch set), so the raw (ScanPrecursorMass - peptideMonoisotopic) difference carries
            // those offsets; left in, they inflate the IQR and make calibration write runaway precursor
            // tolerances (e.g. 1940 ppm). GetObservedMonoisotopicMass removes exactly the offset the search
            // allowed, using the isotope spacing from the deconvolution parameters. In the default
            // monoisotopic mode it is the identity, so baseline calibration behaviour is preserved exactly.
            var precursorErrors = psms.Select(p =>
            {
                double theoreticalMass = p.BioPolymerWithSetModsMonoisotopicMass.Value;
                double observedMass = p.GetObservedMonoisotopicMass(theoreticalMass, commonParameters);
                return (observedMass - theoreticalMass) / theoreticalMass * 1e6;
            }).ToList();
            PsmPrecursorIqrPpmError = precursorErrors.InterquartileRange();
            PsmPrecursorMedianPpmError = precursorErrors.Median();

            var productErrors = psms.Where(p => p.MatchedFragmentIons != null).SelectMany(p => p.MatchedFragmentIons).Select(p => (p.Mz.ToMass(p.Charge) - p.NeutralTheoreticalProduct.NeutralMass) / p.NeutralTheoreticalProduct.NeutralMass * 1e6).ToList();
            PsmProductIqrPpmError = productErrors.InterquartileRange();
            PsmProductMedianPpmError = productErrors.Median();
        }

        public List<LabeledDataPoint> Ms1List { get; }
        public List<LabeledDataPoint> Ms2List { get; }

        public readonly double PsmPrecursorMedianPpmError;
        public readonly double PsmProductMedianPpmError;
        public readonly double PsmPrecursorIqrPpmError;
        public readonly double PsmProductIqrPpmError;
        public readonly List<SpectralMatch> Psms;

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("MS1 calibration datapoint count: " + Ms1List.Count);
            sb.AppendLine("MS1 ppm error median: " + Math.Round(PsmPrecursorMedianPpmError, 3));
            sb.AppendLine("MS1 ppm error interquartile range: " + Math.Round(PsmPrecursorIqrPpmError, 3));

            sb.AppendLine("MS2 calibration datapoint count: " + Ms2List.Count);
            sb.AppendLine("MS2 ppm error median: " + Math.Round(PsmProductMedianPpmError, 3));
            sb.AppendLine("MS2 ppm error interquartile range: " + Math.Round(PsmProductIqrPpmError, 3));
            return sb.ToString();
        }
    }
}