﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class GptmdParameters
    {
        #region Public Constructors

        public GptmdParameters()
        {
            ListOfModsGptmd = GlobalEngineLevelSettings.AllModsKnown.Where(b =>
                b.modificationType.Equals("Glycan") ||
                b.modificationType.Equals("Mod") ||
                b.modificationType.Equals("PeptideTermMod") ||
                b.modificationType.Equals("Metal") ||
                b.modificationType.Equals("ProteinTermMod")).Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();
        }

        #endregion Public Constructors

        #region Public Properties

        public List<Tuple<string, string>> ListOfModsGptmd { get; set; }

        #endregion Public Properties
    }
}