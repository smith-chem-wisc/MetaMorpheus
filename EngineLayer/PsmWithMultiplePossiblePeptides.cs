//using System.Collections.Generic;
//using System.Globalization;
//using System.IO;
//using System.Linq;
//using System.Text;

//namespace EngineLayer
//{
//    public class PsmWithMultiplePossiblePeptides
//    {

//        #region Public Fields

//        public readonly int scanNumber;

//        public readonly int precursorScanNumber;

//        public readonly double scanRetentionTime;

//        public readonly int scanExperimentalPeaks;

//        public readonly double totalIonCurrent;

//        public readonly int scanPrecursorCharge;

//        public readonly double scanPrecursorMZ;

//        public readonly double scanPrecursorMass;

//        public double[] quantIntensity;

//        public double apexMz;

//        public double quantRT;

//        public double mostAbundantMass;

//        #endregion Public Fields

//        #region Public Constructors

//        public PsmWithMultiplePossiblePeptides(Ms2ScanWithSpecificMass scan, Psm match)
//        {
//            this.Score = match.score;
//            this.FileName = scan.FileNameWithoutExtension;
//            this.scanNumber = scan.TheScan.OneBasedScanNumber;
//            this.precursorScanNumber = scan.TheScan.OneBasedPrecursorScanNumber;
//            this.scanRetentionTime = scan.TheScan.RetentionTime;
//            this.scanExperimentalPeaks = scan.TheScan.MassSpectrum.Size;
//            this.totalIonCurrent = scan.TheScan.TotalIonCurrent;
//            this.scanPrecursorCharge = scan.PrecursorCharge;
//            this.scanPrecursorMZ = scan.PrecursorMz;
//            this.scanPrecursorMass = scan.PrecursorMass;
//            ListOfindividualPsmsWithUniquePeptides = new List<Psm> { match };
//            quantIntensity = new double[0];
//        }

//        #endregion Public Constructors

//        #region Public Properties

//        public List<Psm> ListOfindividualPsmsWithUniquePeptides { get; }

//        public double Score { get; }
//        public string FileName { get; }
//        public bool IsDecoy { get { return ListOfindividualPsmsWithUniquePeptides.Any(b => b.IsDecoy); } }
//        public bool IsContaminant { get { return ListOfindividualPsmsWithUniquePeptides.Any(b => b.IsContaminant); } }

//        #endregion Public Properties

//        #region Public Methods

//        public override string ToString()
//        {
//            var sb = new StringBuilder();

//            sb.Append(Path.GetFileNameWithoutExtension(FileName) + '\t');
//            sb.Append(scanNumber.ToString(CultureInfo.InvariantCulture) + '\t');
//            sb.Append(scanRetentionTime.ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            sb.Append(scanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            sb.Append(totalIonCurrent.ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            sb.Append(precursorScanNumber.ToString(CultureInfo.InvariantCulture) + '\t');
//            sb.Append(scanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            sb.Append(scanPrecursorMZ.ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            sb.Append(scanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            sb.Append(Score.ToString("F3", CultureInfo.InvariantCulture) + '\t');
//            sb.Append(string.Join("|", quantIntensity) + '\t');
//            sb.Append(quantRT.ToString("F5", CultureInfo.InvariantCulture) + '\t');

//            StringBuilder thisSb = new StringBuilder();
//            foreach (var psm in ListOfindividualPsmsWithUniquePeptides)
//            {
//                thisSb.Append("[");
//                foreach (var kvp in psm.matchedIonsListPositiveIsMatch)
//                    thisSb.Append("[" + string.Join(",", kvp.Value.Where(b => b > 0).Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
//                thisSb.Append("] ");
//            }
//            if (thisSb.Length > 32000)
//                sb.Append("too many\t");
//            else
//                sb.Append(thisSb.ToString() + '\t');

//            thisSb = new StringBuilder();
//            foreach (var psm in ListOfindividualPsmsWithUniquePeptides)
//            {
//                thisSb.Append("[");
//                thisSb.Append(string.Join(";", psm.matchedIonsListPositiveIsMatch.Select(b => b.Value.Count(c => c > 0))));
//                thisSb.Append("] ");
//            }
//            if (thisSb.Length > 32000)
//                sb.Append("too many\t");
//            else
//                sb.Append(thisSb.ToString() + '\t');

//            //sb.Append(string.Join(";", matchedIonsListPositiveIsMatch.Select(b => b.Value.Count(c => c > 0))) + '\t');

//            //sb.Append("[" + string.Join(",", LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]" + '\t');

//            //sb.Append((LocalizedScores.Max() - score).ToString("F3", CultureInfo.InvariantCulture) + '\t');

//            //if (LocalizedScores.IndexOf(LocalizedScores.Max()) == 0)
//            //    sb.Append("N");
//            //else if (LocalizedScores.IndexOf(LocalizedScores.Max()) == LocalizedScores.Count - 1)
//            //    sb.Append("C");
//            //else
//            //    sb.Append("");

//            //var s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.Accession));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.FullName));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PeptideDescription));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //s = string.Join(" or ", PeptidesWithSetModifications.Select(b => "[" + b.OneBasedStartResidueInProtein + " to " + b.OneBasedEndResidueInProtein + "]"));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PreviousAminoAcid));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.NextAminoAcid));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.BaseSequence));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Sequence));
//            //if (s.Length > 32000)
//            //    s = "too many";
//            //sb.Append(s + "\t");

//            //sb.Append(NumVariableMods.ToString(CultureInfo.InvariantCulture) + '\t');
//            //sb.Append(MissedCleavages.ToString(CultureInfo.InvariantCulture) + '\t');
//            //sb.Append(PeptideMonoisotopicMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            //sb.Append((ScanPrecursorMass - PeptideMonoisotopicMass).ToString("F5", CultureInfo.InvariantCulture) + '\t');
//            //sb.Append(((ScanPrecursorMass - PeptideMonoisotopicMass) / PeptideMonoisotopicMass * 1e6).ToString("F5", CultureInfo.InvariantCulture) + '\t');

//            if (IsDecoy)
//                sb.Append("D");
//            else if (IsContaminant)
//                sb.Append("C");
//            else
//                sb.Append("T");
//            return sb.ToString();
//        }

//        #endregion Public Methods

//        #region Internal Methods

//        internal static string GetTabSeparatedHeader()
//        {
//            var sb = new StringBuilder();
//            sb.Append("fileName" + '\t');
//            sb.Append("scanNumber" + '\t');
//            sb.Append("scanRetentionTime" + '\t');
//            sb.Append("scanExperimentalPeaks" + '\t');
//            sb.Append("totalIonCurrent" + '\t');
//            sb.Append("precursorScanNumber" + '\t');
//            sb.Append("scanPrecursorCharge" + '\t');
//            sb.Append("scanPrecursorMZ" + '\t');
//            sb.Append("scanPrecursorMass" + '\t');
//            sb.Append("score" + '\t');
//            sb.Append("notch" + '\t');
//            sb.Append("quantificationIntensity" + '\t');
//            sb.Append("quantificationRT" + '\t');

//            sb.Append("matched ions" + '\t');
//            sb.Append("matched ion counts" + '\t');
//            sb.Append("localized scores" + '\t');
//            sb.Append("improvement" + '\t');
//            sb.Append("terminal localization");

//            sb.Append("Protein Accession" + '\t');
//            sb.Append("Protein FullName" + '\t');
//            sb.Append("Peptide Description" + '\t');
//            sb.Append("Start and End ResidueInProtein" + '\t');
//            sb.Append("PreviousAminoAcid" + '\t');
//            sb.Append("NextAminoAcid" + '\t');
//            sb.Append("BaseSequence" + '\t');
//            sb.Append("FullSequence" + '\t');
//            sb.Append("numVariableMods" + '\t');
//            sb.Append("MissedCleavages" + '\t');
//            sb.Append("PeptideMonoisotopicMass" + '\t');
//            sb.Append("MassDiff (Da)" + '\t');
//            sb.Append("MassDiff (ppm)" + '\t');
//            sb.Append("Decoy/Contaminant/Target");
//            return sb.ToString();
//        }

//        #endregion Internal Methods

//    }
//}