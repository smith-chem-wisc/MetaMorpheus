using Chemistry;
using EngineLayer;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public class ModList
    {

        #region Private Fields

        public readonly List<ModificationWithLocation> Mods;

        private string FullFileName;

        public ModList(string modFile, List<ModificationWithLocation> Mods)
        {
            this.FullFileName = modFile;
            this.Mods = Mods;
        }

        #endregion Private Fields

        #region Public Constructors


        #endregion Public Constructors

        #region Public Properties

        public string FileName
        {
            get
            {
                return Path.GetFileName(FullFileName);
            }
        }

        public string Description { get; private set; }

        #endregion Public Properties

    }
}