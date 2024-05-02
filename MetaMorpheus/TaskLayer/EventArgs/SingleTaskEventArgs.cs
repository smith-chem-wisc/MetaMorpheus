using System;

namespace TaskLayer
{
    [Serializable]
    public class SingleTaskEventArgs : EventArgs
    {
        public SingleTaskEventArgs(string displayName)
        {
            this.DisplayName = displayName;
        }

        public string DisplayName { get; private set; }
    }
}