using TaskLayer;

namespace MetaMorpheusGUI
{
    public class InRunTask : ForTreeView
    {
        public readonly MetaMorpheusTask task;

        public InRunTask(string displayName, MetaMorpheusTask task) : base(displayName, displayName)
        {
            this.task = task;
        }
    }
}