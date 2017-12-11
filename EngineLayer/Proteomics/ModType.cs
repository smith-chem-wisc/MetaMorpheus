using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Text;

namespace EngineLayer
{
    public enum ModType
    {
        [Description("Common Fixed")]
        CommonFixed,
        [Description("Common Variable")]
        CommonVariable,
        [Description("Detached")]
        Detached,
        [Description("Glycan")]
        Glycan,
        [Description("missing")]
        missing,
        [Description("Mod")]
        Mod,
        [Description("ProteinTermMod")]
        ProteinTermMod,
        [Description("PeptideTermMod")]
        PeptideTermMod,
        [Description("Metal")]
        Metal,
        [Description("TrypsinDigestedMod")]
        TrypsinDigestedMod,
        [Description("1 nucelotide substitution")]
        OneNucelotideSubstitution,
        [Description("2+ nucelotide substitution")]
        TwoPlusNucelotideSubstitution,
        [Description("Surfactant")]
        Surfactant,
        [Description("TandemMassTag")]
        TandemMassTag,
        [Description("Unimod")]
        Unimod,
        [Description("Uniprot")]
        Uniprot

    }

    

}