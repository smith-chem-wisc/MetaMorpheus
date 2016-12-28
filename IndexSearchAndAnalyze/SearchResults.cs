using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IndexSearchAndAnalyze
{
    public class SearchResults : MyResults
    {
        public SearchResults(NewPsm[][] newPsms)
        {
            this.newPsms = newPsms;
        }

        public NewPsm[][] newPsms { get; private set; }
    }
}
