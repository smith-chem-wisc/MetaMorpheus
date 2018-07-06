using TaskLayer;

namespace MetaMorpheusGUI
{
    internal class PreRunTask
    {
        public readonly MetaMorpheusTask MetaMorpheusTask;

        public PreRunTask(MetaMorpheusTask theTask)
        {
            MetaMorpheusTask = theTask;
        }

        public string DisplayName { get; set; }

        public PreRunTask Clone()
        {
            return (PreRunTask)MemberwiseClone();
        }
    }
}