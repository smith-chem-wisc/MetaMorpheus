using TaskLayer;

namespace MetaMorpheusGUI
{
    internal class PreRunTask
    {
        #region Public Fields

        public readonly MetaMorpheusTask metaMorpheusTask;

        #endregion Public Fields

        #region Public Constructors

        public PreRunTask(MetaMorpheusTask theTask)
        {
            metaMorpheusTask = theTask;
        }

        #endregion Public Constructors

        #region Public Properties

        public string DisplayName { get; set; }

        public bool IsRnaTask => metaMorpheusTask.CommonParameters.DetermineAnalyteType() == EngineLayer.AnalyteType.Oligo;

        public bool IsModeAgnosticTask => metaMorpheusTask is SpectralAveragingTask;

        #endregion Public Properties

        #region Public Methods

        public PreRunTask Clone()
        {
            return (PreRunTask)this.MemberwiseClone();
        }

        #endregion Public Methods
    }
}