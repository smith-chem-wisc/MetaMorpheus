using System.Globalization;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class NewPsmWithFDR
    {
        #region Public Fields

        public PSMwithTargetDecoyKnown thisPSM;

        #endregion Public Fields

        #region Private Fields

        private int cumulative_target;
        private int cumulative_decoy;

        #endregion Private Fields

        #region Internal Constructors

        internal NewPsmWithFDR(PSMwithTargetDecoyKnown thisPSM, int cumulative_target, int cumulative_decoy, double temp_q_value)
        {
            this.thisPSM = thisPSM;
            this.cumulative_target = cumulative_target;
            this.cumulative_decoy = cumulative_decoy;
            this.QValue = temp_q_value;
        }

        #endregion Internal Constructors

        #region Public Properties

        public double QValue { get; set; }

        public bool isDecoy
        {
            get
            {
                return thisPSM.isDecoy;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public static string GetTabSeparatedHeader()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(PSMwithTargetDecoyKnown.GetTabSeparatedHeader() + '\t');
            sb.Append("Decoy" + '\t');
            sb.Append("cumulative_target" + '\t');
            sb.Append("cumulative_decoy" + '\t');
            sb.Append("QValue");
            return sb.ToString();
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

        #endregion Public Methods
    }
}