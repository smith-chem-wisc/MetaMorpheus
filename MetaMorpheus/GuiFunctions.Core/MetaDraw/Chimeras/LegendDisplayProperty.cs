namespace GuiFunctions.MetaDraw;

/// <summary>
/// Which property of a chimeric identification the legend shows. Extracted from
/// ChimeraLegendCanvas so MetaDrawSettings does not depend on a WPF Canvas subclass.
/// </summary>
public enum LegendDisplayProperty
{
    ProteinName,
    ProteinAccession,
    BaseSequence,
    FullSequence,
    Modifications
}
