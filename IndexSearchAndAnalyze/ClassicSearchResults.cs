using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IndexSearchAndAnalyze
{
    public class ClassicSearchResults : MyResults
    {
        public ClassicSearchResults(MyParams s, List<NewPsm> newPsms) : base(s)
        {
            this.newPsms =  newPsms;
        }

        public List<NewPsm> newPsms { get; private set; }
        public override string ToString()
        {
            var sp = (ClassicSearchParams)s;
            StringBuilder sb = new StringBuilder();
            sb.Append("ClassicSearchResults: ");
            sb.Append(base.ToString());
            sb.AppendLine();
            sb.Append("Total psms: " + newPsms.Count);
            sb.AppendLine();

            return sb.ToString();
        }
    }
}