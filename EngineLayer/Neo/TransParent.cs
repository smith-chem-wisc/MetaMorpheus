using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class TransParent : Parent
    {
        #region Public Constructors

        public TransParent(string id, string seq, List<int> start, int length, ParentInfo.terminal terminal) : base(id, seq)
        {
            this.start = start;
            this.peptideLength = length;
            this.terminal = terminal;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<int> start { get; set; }
        public int peptideLength { get; set; }
        public ParentInfo.terminal terminal { get; set; }

        #endregion Public Properties
    }
}