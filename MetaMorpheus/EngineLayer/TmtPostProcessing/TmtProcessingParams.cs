using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.TmtPostProcessing
{
    public class TmtProcessingParams
    {
        IsobaricMassTagType TmtTagType;
        TmtReferenceType ReferenceType;
        string? ReferenceChannelName;
        bool NormalizeLoading;
        bool TrimLowIntensityPsms;
        bool UseMedianPeptideForProtein;

        public TmtProcessingParams(IsobaricMassTagType tmtTagType, TmtReferenceType referenceType, string? referenceChannelName = null,
            bool normalizeLoading = true,  bool trimLowIntensityPsms = true, bool useMedianPeptideForProtein = true)
        {
            TmtTagType = tmtTagType;
            ReferenceType = referenceType;
            ReferenceChannelName = referenceChannelName;
            NormalizeLoading = normalizeLoading;
            TrimLowIntensityPsms = trimLowIntensityPsms;
            UseMedianPeptideForProtein = useMedianPeptideForProtein;
        }
    }
    public enum TmtReferenceType
    {
        VirtualReference,
        BridgeReference
    }
}
