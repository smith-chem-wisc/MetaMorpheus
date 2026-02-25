// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

namespace EngineLayer.FdrAnalysisDia
{
    /// <summary>
    /// Which score to use as the primary ranking criterion for DIA FDR calculation.
    /// </summary>
    public enum DiaFdrScoreType
    {
        /// <summary>Normalized dot product between library and extracted intensities</summary>
        DotProduct,

        /// <summary>Spectral angle: 1 - (2/π) * arccos(normalized dot product)</summary>
        SpectralAngle
    }
}