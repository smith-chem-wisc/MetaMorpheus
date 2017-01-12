using System.Text;

namespace InternalLogicEngineLayer
{
    public class ClassicSearchResults : MyResults
    {
        public ClassicSpectrumMatch[][] outerPsms { get; private set; }

        internal ClassicSearchResults(ClassicSearchEngine searchParams, ClassicSpectrumMatch[][] outerPsms) : base(searchParams)
        {
            this.outerPsms = outerPsms;
        }

        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("\t\tSome search result");
            return sb.ToString();
        }
    }
}