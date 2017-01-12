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

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("ClassicSearchResults:");
            sb.Append(base.ToString());
            return sb.ToString();
        }
    }
}