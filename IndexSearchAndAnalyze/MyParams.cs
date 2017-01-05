using System;

namespace IndexSearchAndAnalyze
{
    public abstract class MyParams
    {
        public Action<string> outputAction;
        public Action<int> progressAction;

        public MyParams(Action<string> a1, Action<int> a2)
        {
            this.outputAction = a1;
            this.progressAction = a2;
        }

        internal abstract void Validate();
    }
}