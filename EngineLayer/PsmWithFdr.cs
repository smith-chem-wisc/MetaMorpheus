using System.Globalization;
using System.Text;

namespace EngineLayer
{
    public class NewPsmWithFdr
    {

        #region Public Fields

        public PsmParent thisPSM;

        #endregion Public Fields

        #region Private Fields

        private int cumulativeTarget;
        private int cumulativeDecoy;
        private int cumulativeTargetNotch;
        private int cumulativeDecoyNotch;

        #endregion Private Fields

        #region Public Constructors

        public NewPsmWithFdr(PsmParent thisPsm)
        {
            this.thisPSM = thisPsm;
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PsmParent.GetTabSeparatedHeader() + '\t');
                sb.Append("cumulative_target" + '\t');
                sb.Append("cumulative_decoy" + '\t');
                sb.Append("QValue" + '\t');
                sb.Append("cumulative_target_notch" + '\t');
                sb.Append("cumulative_decoy_notch" + '\t');
                sb.Append("QValue_notch");
                return sb.ToString();
            }
        }

        public double QValue { get; set; }
        public double QValueNotch { get; set; }

        public bool IsDecoy
        {
            get
            {
                return thisPSM.Pli.IsDecoy;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(thisPSM.ToString() + '\t');
            sb.Append(cumulativeTarget.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulativeDecoy.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(QValue.ToString("F6", CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulativeTargetNotch.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulativeDecoyNotch.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(QValueNotch.ToString("F6", CultureInfo.InvariantCulture));

            return sb.ToString();
        }

        public void SetValues(int cumulativeTarget, int cumulativeDecoy, double tempQValue, int cumulativeTargetNotch, int cumulativeDecoyNotch, double tempQValueNotch)
        {
            this.cumulativeTarget = cumulativeTarget;
            this.cumulativeDecoy = cumulativeDecoy;
            QValue = tempQValue;

            this.cumulativeTargetNotch = cumulativeTargetNotch;
            this.cumulativeDecoyNotch = cumulativeDecoyNotch;
            QValueNotch = tempQValueNotch;
        }

        #endregion Public Methods

    }
}