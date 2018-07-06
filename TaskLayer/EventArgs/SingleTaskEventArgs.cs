using System;

namespace TaskLayer
{
    public class SingleTaskEventArgs : EventArgs
    {
        public SingleTaskEventArgs(string displayName)
        {
            DisplayName = displayName;
        }

        public string DisplayName { get; private set; }
    }
}