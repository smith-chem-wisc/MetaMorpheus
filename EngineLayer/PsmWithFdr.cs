using System.Globalization;
using System.Text;

namespace EngineLayer
{
    public class NewPsmWithFdr
    {

        #region Public Fields

        public PsmWithMultiplePossiblePeptides thisPSM;

        #endregion Public Fields

        #region Private Fields

        private readonly int cumulativeTarget;
        private readonly int cumulativeDecoy;

        #endregion Private Fields

        #region Public Constructors

        public NewPsmWithFdr(PsmWithMultiplePossiblePeptides thisPsm, int cumulativeTarget, int cumulativeDecoy, double tempQValue)
        {
            this.thisPSM = thisPsm;
            this.cumulativeTarget = cumulativeTarget;
            this.cumulativeDecoy = cumulativeDecoy;
            qValue = tempQValue;
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PsmWithMultiplePossiblePeptides.TabSeparatedHeader + '\t');
                sb.Append("cumulative_target" + '\t');
                sb.Append("cumulative_decoy" + '\t');
                sb.Append("QValue");
                return sb.ToString();
            }
        }

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

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(thisPSM.ToString() + '\t');
            sb.Append(cumulativeTarget.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(cumulativeDecoy.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(qValue.ToString("F5", CultureInfo.InvariantCulture));

            return sb.ToString();
        }

        #endregion Public Methods

    }
}