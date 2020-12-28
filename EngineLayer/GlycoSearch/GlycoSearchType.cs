using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.GlycoSearch
{
    //Localization of multiple glycans on one peptides can be divide into the following groups based on the quanlity of the localization. Similar to Proteomform Level.
    public enum LocalizationLevel
    {
        Level1,
        Level1b,
        Level2,
        Level3
    }

    public enum GlycoType
    {
        SinglePep,
        OGlycoPep,
        NGlycoPep,
        MixedGlycoPep
    }

    public enum GlycoSearchType
    {
        OGlycanSearch,
        NGlycanSearch,
        N_O_GlycanSearch

    }
}
