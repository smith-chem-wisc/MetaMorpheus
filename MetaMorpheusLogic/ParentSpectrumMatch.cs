using OldInternalLogic;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace MetaMorpheusLogic
{
    public abstract class ParentSpectrumMatch
    {
        internal Dictionary<ProductType, double[]> matchedIonsList;
        internal int scanNumber;
        internal int scanPrecursorCharge;
        internal double scanPrecursorMass;

        public List<double> LocalizedScores { get; internal set; }
        protected CompactPeptide compactPeptide;
        public double Score { get; protected set; }

        internal abstract CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications);

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            //sb.Append(spectraFileIndex.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanNumber.ToString(CultureInfo.InvariantCulture) + '\t');
            //sb.Append(scanRT.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            //sb.Append(scanPrecursorMZ.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            // sb.Append(scanPrecursorIntensity.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            //sb.Append(scanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            // sb.Append(TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            //sb.Append(ScoreFromSearch.ToString("F3", CultureInfo.InvariantCulture) + '\t');

            sb.Append("[");
            foreach (var kvp in matchedIonsList)
                sb.Append("[" + string.Join(",", kvp.Value.Where(b => b > 0).Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
            sb.Append("]" + '\t');

            sb.Append(string.Join(";", matchedIonsList.Select(b => b.Value.Where(c => c > 0).Count())) + '\t');

            sb.Append("[" + string.Join(",", LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]" + '\t');

            sb.Append((LocalizedScores.Max() - Score).ToString("F3", CultureInfo.InvariantCulture) + '\t');
            if (LocalizedScores.IndexOf(LocalizedScores.Max()) == 0)
                sb.Append("N");
            else if (LocalizedScores.IndexOf(LocalizedScores.Max()) == LocalizedScores.Count - 1)
                sb.Append("C");
            else
                sb.Append("");

            return sb.ToString();
        }
    }
}