using System.Globalization;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class NewPsmWithFdr
    {
        #region Public Fields

        public PSMwithTargetDecoyKnown thisPSM;

        #endregion Public Fields

        #region Private Fields

        private readonly int cumulativeTarget;
        private readonly int cumulativeDecoy;

        #endregion Private Fields

        #region Public Constructors

        public NewPsmWithFdr(PSMwithTargetDecoyKnown thisPsm, int cumulativeTarget, int cumulativeDecoy, double tempQValue)
        {
            this.thisPSM = thisPsm;
            this.cumulativeTarget = cumulativeTarget;
            this.cumulativeDecoy = cumulativeDecoy;
            qValue = tempQValue;
        }

        #endregion Public Constructors

        #region Public Properties

        public double qValue { get; set; }

        public bool IsDecoy
        {
            get
            {
                return thisPSM.IsDecoy;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PSMwithTargetDecoyKnown.TabSeparatedHeader + '\t');
                sb.Append("Decoy" + '\t');
                sb.Append("cumulative_target" + '\t');
                sb.Append("cumulative_decoy" + '\t');
                sb.Append("QValue");
                return sb.ToString();
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(thisPSM.ToString() + '\t');
            sb.Append(IsDecoy.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulativeTarget.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulativeDecoy.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(qValue.ToString("F5", CultureInfo.InvariantCulture));

            return sb.ToString();
        }

        #endregion Public Methods
    }
}