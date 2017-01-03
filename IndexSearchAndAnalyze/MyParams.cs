using System;

namespace IndexSearchAndAnalyze
{
    public abstract class MyParams
    {
        public event EventHandler<string> outputHandler;

        public event EventHandler<int> progressHandler;

        internal abstract void Validate();

        public void OnOutput(string e)
        {
            outputHandler?.Invoke(this, e);
        }

        public void OnProgress(int e)
        {
            progressHandler?.Invoke(this, e);
        }
    }
}