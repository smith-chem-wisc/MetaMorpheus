using EngineLayer;
using System.Collections.Generic;
using TaskLayer;

namespace MetaMorpheusGUI
{
    class ModListForCalibrationTask
    {

        #region Private Fields

        private readonly ModList uu;

        #endregion Private Fields

        #region Public Constructors

        public ModListForCalibrationTask(ModList uu)
        {
            this.uu = uu;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Fixed { get; set; }
        public bool Variable { get; set; }
        public bool Localize { get; set; }

        public string FileName
        {
            get
            {
                return uu.FileName;
            }
        }

        public string Description { get { return uu.Description; } }

        public IEnumerable<MetaMorpheusModification> Mods
        {
            get
            {
                return uu.Mods;
            }
        }

        public ModList ModList { get { return uu; } }

        #endregion Public Properties

    }
}