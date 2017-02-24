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
        public readonly double scanPrecursorIntensity;
        public readonly int scanPrecursorCharge;
        public readonly double scanPrecursorMZ;
        public readonly double scanPrecursorMass;
        public readonly double score;
        public double apexIntensity;

        public Dictionary<ProductType, double[]> matchedIonsList;
        public List<double> LocalizedScores;

        #endregion Public Fields

        #region Internal Fields

        internal readonly int notch;

        #endregion Internal Fields

        #region Protected Constructors

        protected PsmParent(string fileName, double scanRetentionTime, double scanPrecursorIntensity, double scanPrecursorMass, int scanNumber, int precursorScanNumber, int scanPrecursorCharge, int scanExperimentalPeaks, double totalIonCurrent, double scanPrecursorMZ, double score, int notch)
        {
            this.fileName = fileName;
            this.scanNumber = scanNumber;
            this.precursorScanNumber = precursorScanNumber;
            this.scanRetentionTime = scanRetentionTime;
            this.scanExperimentalPeaks = scanExperimentalPeaks;
            this.totalIonCurrent = totalIonCurrent;
            this.scanPrecursorIntensity = scanPrecursorIntensity;
            this.scanPrecursorCharge = scanPrecursorCharge;
            this.scanPrecursorMZ = scanPrecursorMZ;
            this.scanPrecursorMass = scanPrecursorMass;
            this.score = score;
            this.notch = notch;
        }

        #endregion Protected Constructors

        #region Public Methods

        public abstract CompactPeptide GetCompactPeptide(List<ModificationWithMass> variableModifications, List<ModificationWithMass> localizeableModifications, List<ModificationWithMass> fixedModifications);

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append(Path.GetFileNameWithoutExtension(fileName) + '\t');
            sb.Append(scanNumber.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanRetentionTime.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(totalIonCurrent.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorIntensity.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorMZ.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(score.ToString("F3", CultureInfo.InvariantCulture) + '\t');
            sb.Append(notch.ToString("F3", CultureInfo.InvariantCulture) + '\t');
            sb.Append(apexIntensity.ToString("F3", CultureInfo.InvariantCulture) + '\t');

            sb.Append("[");
            foreach (var kvp in matchedIonsList)
                sb.Append("[" + string.Join(",", kvp.Value.Where(b => b > 0).Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
            sb.Append("]" + '\t');

            sb.Append(string.Join(";", matchedIonsList.Select(b => b.Value.Count(c => c > 0))) + '\t');

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
            sb.Append("scanPrecursorIntensity" + '\t');
            sb.Append("scanPrecursorCharge" + '\t');
            sb.Append("scanPrecursorMZ" + '\t');
            sb.Append("scanPrecursorMass" + '\t');
            sb.Append("score" + '\t');
            sb.Append("notch" + '\t');
            sb.Append("apexIntensity" + '\t');

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