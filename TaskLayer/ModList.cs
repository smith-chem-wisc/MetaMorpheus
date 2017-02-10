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

        public ModList(string modFile)
        {
            this.FullFileName = modFile;
            this.Mods = UsefulProteomicsDatabases.PtmListLoader.ReadMods(modFile).ToList();
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