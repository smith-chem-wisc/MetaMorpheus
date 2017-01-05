using Chemistry;
using MassSpectrometry;
using MetaMorpheus;
using Spectra;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace IndexSearchAndAnalyze
{
    public class NewPsm
    {
        public double scanPrecursorMass { get; private set; }
        public int scanNumber { get; private set; }
        public double ScoreFromSearch { get; private set; }
        public int spectraFileIndex { get; private set; }
        public CompactPeptide peptide { get; internal set; }
        public double scanRT { get; private set; }
        public double scanPrecursorMZ { get; private set; }
        public int scanPrecursorCharge { get; private set; }
        public double scanPrecursorIntensity { get; private set; }
        public int scanExperimentalPeaks { get; private set; }
        public double TotalIonCurrent { get; private set; }

        public List<double> LocalizedScores;

        public Dictionary<ProductType, double[]> matchedIonsList;

        public double ScoreFromMatch;

        public NewPsm(double scanPrecursorMZ, int scanNumber, double scanRT,int scanPrecursorCharge, int scanExperimentalPeaksCount, double totalIonCurrent, double precursorIntensity, int spectraFileIndex, CompactPeptide theBestPeptide, double score)
        {
            this.scanPrecursorMZ = scanPrecursorMZ;
            this.scanNumber = scanNumber;
            this.scanPrecursorCharge = scanPrecursorCharge;
            this.scanRT = scanRT;
            this.scanPrecursorMass = scanPrecursorMZ.ToMass(scanPrecursorCharge);
            this.scanPrecursorIntensity = precursorIntensity;
            this.scanExperimentalPeaks = scanExperimentalPeaksCount;
            this.TotalIonCurrent = totalIonCurrent;
            this.ScoreFromSearch = score;
            this.spectraFileIndex = spectraFileIndex;
            this.peptide = theBestPeptide;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(spectraFileIndex.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanNumber.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanRT.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorMZ.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorIntensity.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(scanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(ScoreFromSearch.ToString("F3", CultureInfo.InvariantCulture) + '\t');

            sb.Append("[");
            foreach (var kvp in matchedIonsList)
                sb.Append("[" + string.Join(",", kvp.Value.Where(b => b > 0).Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
            sb.Append("]" + '\t');

            sb.Append(string.Join(";", matchedIonsList.Select(b => b.Value.Where(c => c > 0).Count())) + '\t');

            sb.Append("[" + string.Join(",", LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]" + '\t');

            sb.Append((LocalizedScores.Max() - ScoreFromMatch).ToString("F3", CultureInfo.InvariantCulture) + '\t');
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