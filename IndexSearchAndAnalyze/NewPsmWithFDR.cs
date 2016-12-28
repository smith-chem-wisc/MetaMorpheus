using MetaMorpheus;
using System.Globalization;
using System.Text;

namespace IndexSearchAndAnalyze
{
    public class NewPsmWithFDR
    {
        private int cumulative_target;
        private int cumulative_decoy;
        public double QValue { get; set; }
        public PSMwithPeptide thisPSM;
        public bool isDecoy
        {
            get
            {
                return thisPSM.isDecoy;
            }
        }

        public NewPsmWithFDR(PSMwithPeptide thisPSM, int cumulative_target, int cumulative_decoy, double temp_q_value)
        {
            this.thisPSM = thisPSM;
            this.cumulative_target = cumulative_target;
            this.cumulative_decoy = cumulative_decoy;
            this.QValue = temp_q_value;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(thisPSM.ToString() + '\t');
            sb.Append(isDecoy.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulative_target.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulative_decoy.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(QValue.ToString("F5", CultureInfo.InvariantCulture));

            return sb.ToString();
        }
    }
}