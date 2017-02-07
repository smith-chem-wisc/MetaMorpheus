using TaskLayer;

namespace MetaMorpheusGUI
{
    internal class MetaMorpheusTaskForDataGrid
    {
        public string Task { get { return metaMorpheusTask.GetType().Name; } }
        public bool IsMySelected { get; set; }
        public readonly MetaMorpheusTask metaMorpheusTask;

        public MetaMorpheusTaskForDataGrid(MetaMorpheusTask theTask)
        {
            this.metaMorpheusTask = theTask;
        }
    }
}