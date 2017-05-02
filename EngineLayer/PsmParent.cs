using Proteomics;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public abstract class PsmParent
    {

        #region Public Fields

        public readonly string fileName;
        public readonly int scanNumber;
        public readonly int precursorScanNumber;
        public readonly double scanRetentionTime;
        public readonly int scanExperimentalPeaks;
        public readonly double totalIonCurrent;
        public readonly int scanPrecursorCharge;
        public readonly double scanPrecursorMZ;
        public readonly double scanPrecursorMass;
        public readonly double score;
        public double[] quantIntensity;
        public double apexMz;
        public double quantRT;
        public double mostAbundantMass;

        public Dictionary<ProductType, double[]> matchedIonsListPositiveIsMatch;
        public List<double> LocalizedScores;

        #endregion Public Fields

        #region Internal Fields

        internal readonly int notch;
        internal double? precursorScanBestMass;

        #endregion Internal Fields

        #region Protected Constructors

        protected PsmParent(Ms2ScanWithSpecificMass scan, double score, int notch)
        {
            this.fileName = scan.FileNameWithoutExtension;
            this.scanNumber = scan.TheScan.OneBasedScanNumber;
            this.precursorScanNumber = scan.TheScan.OneBasedPrecursorScanNumber;
            this.scanRetentionTime = scan.TheScan.RetentionTime;
            this.scanExperimentalPeaks = scan.TheScan.MassSpectrum.Size;
            this.totalIonCurrent = scan.TheScan.TotalIonCurrent;
            this.scanPrecursorCharge = scan.PrecursorCharge;
            this.scanPrecursorMZ = scan.PrecursorMz;
            this.scanPrecursorMass = scan.PrecursorMass;
            this.score = score;
            this.notch = notch;
            quantIntensity = new double[1];
        }

        #endregion Protected Constructors

        #region Public Methods

        public abstract CompactPeptide GetCompactPeptide(Dictionary<ModificationWithMass, ushort> modsDictionary);

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append(Path.GetFileNameWithoutExtension(fileName) + '\t');
            sb.Append(scanNumber.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanRetentionTime.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(totalIonCurrent.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(precursorScanNumber.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorMZ.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(score.ToString("F3", CultureInfo.InvariantCulture) + '\t');
            sb.Append(notch.ToString("F3", CultureInfo.InvariantCulture) + '\t');
            sb.Append(string.Join("|", quantIntensity) + '\t');
            sb.Append(quantRT.ToString("F5", CultureInfo.InvariantCulture) + '\t');

            sb.Append("[");
            foreach (var kvp in matchedIonsListPositiveIsMatch)
                sb.Append("[" + string.Join(",", kvp.Value.Where(b => b > 0).Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
            sb.Append("]" + '\t');

            sb.Append(string.Join(";", matchedIonsListPositiveIsMatch.Select(b => b.Value.Count(c => c > 0))) + '\t');

            sb.Append("[" + string.Join(",", LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]" + '\t');

            sb.Append((LocalizedScores.Max() - score).ToString("F3", CultureInfo.InvariantCulture) + '\t');

            if (LocalizedScores.IndexOf(LocalizedScores.Max()) == 0)
                sb.Append("N");
            else if (LocalizedScores.IndexOf(LocalizedScores.Max()) == LocalizedScores.Count - 1)
                sb.Append("C");
            else
                sb.Append("");

            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("fileName" + '\t');
            sb.Append("scanNumber" + '\t');
            sb.Append("scanRetentionTime" + '\t');
            sb.Append("scanExperimentalPeaks" + '\t');
            sb.Append("totalIonCurrent" + '\t');
            sb.Append("precursorScanNumber" + '\t');
            sb.Append("scanPrecursorCharge" + '\t');
            sb.Append("scanPrecursorMZ" + '\t');
            sb.Append("scanPrecursorMass" + '\t');
            sb.Append("score" + '\t');
            sb.Append("notch" + '\t');
            sb.Append("quantificationIntensity" + '\t');
            sb.Append("quantificationRT" + '\t');

            sb.Append("matched ions" + '\t');
            sb.Append("matched ion counts" + '\t');
            sb.Append("localized scores" + '\t');
            sb.Append("improvement" + '\t');
            sb.Append("terminal localization");
            return sb.ToString();
        }

        #endregion Internal Methods

    }
}