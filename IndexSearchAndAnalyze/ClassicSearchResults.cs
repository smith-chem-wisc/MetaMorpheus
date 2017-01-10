using System.Text;

namespace IndexSearchAndAnalyze
{
    public class ClassicSearchResults : MyResults
    {
        internal ClassicSpectrumMatch[][] outerPsms { get; private set; }

        internal ClassicSearchResults(ClassicSearchParams searchParams, ClassicSpectrumMatch[][] outerPsms) : base(searchParams)
        {
            this.outerPsms = outerPsms;
        }

        public override string ToString()
        {
            var sp = (ClassicSearchParams)s;
            StringBuilder sb = new StringBuilder();
            sb.Append("ClassicSearchResults: ");
            sb.Append(base.ToString());

            return sb.ToString();
        }
    }
}