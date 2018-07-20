using TaskLayer;

namespace MetaMorpheusGUI
{
    public class InRunTask : ForTreeView
    {
        public readonly MetaMorpheusTask Task;

        public InRunTask(string displayName, MetaMorpheusTask task) : base(displayName, displayName)
        {
            Task = task;
        }
    }
}