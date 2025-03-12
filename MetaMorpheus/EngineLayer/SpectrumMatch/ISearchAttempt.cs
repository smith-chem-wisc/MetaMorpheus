using System;

namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// An interface used during the search process by the SearchLog stored within the SpectralMatch. 
    /// </summary>
    public interface ISearchAttempt : IEquatable<ISearchAttempt>
    {
        double Score { get; }
        bool IsDecoy { get; }
        int Notch { get; }
    }
}
