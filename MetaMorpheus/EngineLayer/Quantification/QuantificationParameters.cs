using System.Collections.Generic;
using Omics.Digestion;
using Omics.Modifications;

namespace EngineLayer.Quantification;

public class QuantificationParameters
{
    // Base Properties
    public string SearchTaskId { get; set; }
    public List<string> CurrentRawFileList { get; set; }
    public string OutputFolder { get; set; }


    // Quantification Parameters
    public HashSet<IDigestionParams> ListOfDigestionParams { get; set; }
    public bool Normalize { get; set; }
    public double QuantifyPpmTol { get; set; }
    public bool UseSharedPeptides { get; set; }


    // LFQ Parameters
    public bool MatchBetweenRuns { get; set; }
    public double MbrFdrThreshold { get; set; }


    // Multiplex Parameters - Optional
    public Modification MultiplexMod { get; set; }

    // Silac Parameters - Optional
    public List<SilacLabel> SilacLabels { get; set; }
    public SilacLabel StartTurnoverLabel { get; set; } //used for SILAC turnover experiments
    public SilacLabel EndTurnoverLabel { get; set; } //used for SILAC turnover experiments
    public Dictionary<string, int> ModsToWriteSelection { get; set; }
}