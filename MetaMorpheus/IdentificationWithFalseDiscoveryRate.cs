using System.Globalization;
using System.Text;

namespace OldInternalLogic
{
    public class IdentificationWithFalseDiscoveryRate<PP>
    {
        public PeptideSpectrumMatch Identification { get; private set; }
        public int CumulativeTarget { get; private set; }
        public int CumulativeDecoy { get; private set; }
        public double QValue { get; set; }

        public static readonly string Header = "\tDecoy?\tCumulative Target\tCumulative Decoy\tQ-Value (%)";

        public IdentificationWithFalseDiscoveryRate(PeptideSpectrumMatch identification, int cumulativeTarget, int cumulativeDecoy, double qValue)
        {
            Identification = identification;
            CumulativeTarget = cumulativeTarget;
            CumulativeDecoy = cumulativeDecoy;
            QValue = qValue;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder(Identification.ToString() + '\t');

            sb.Append(Identification.isDecoy.ToString() + '\t');
            sb.Append(CumulativeTarget.ToString() + '\t');
            sb.Append(CumulativeDecoy.ToString() + '\t');
            sb.Append((QValue * 100.0).ToString(CultureInfo.InvariantCulture));

            return sb.ToString();
        }
    }
}