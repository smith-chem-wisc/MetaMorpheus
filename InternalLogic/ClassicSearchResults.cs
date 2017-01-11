using System.Text;

namespace InternalLogic
{
    public class ClassicSearchResults : MyResults
    {
        public ClassicSpectrumMatch[][] outerPsms { get; private set; }

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