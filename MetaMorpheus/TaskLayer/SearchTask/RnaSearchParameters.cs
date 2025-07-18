using System.Collections.Generic;

namespace TaskLayer;

/// <summary>
/// Inherits all params and properties from normal search parameters
/// </summary>
public class RnaSearchParameters : SearchParameters
{
    public RnaSearchParameters() : base()
    {
        // Have not been generalized/Optimized for RNA
        DoLocalizationAnalysis = false;
        DoLabelFreeQuantification = false;
        MinAllowedInternalFragmentLength = 0;

        SearchType = SearchType.Classic;

        // Output Options
        WritePepXml = false;
        WriteMzId = false;
        UpdateSpectralLibrary = false;
        WriteSpectralLibrary = false;
        ModsToWriteSelection = new Dictionary<string, int>
        {
            //Key is modification type.

            //Value is integer 0, 1, 2 and 3 interpreted as:
            //   0:   Do not Write
            //   1:   Write if in DB and Observed
            //   2:   Write if in DB
            //   3:   Write if Observed
            {"Biological", 3},
            {"Digestion Termini", 3},
            {"Metal", 3},
        };
    }
}