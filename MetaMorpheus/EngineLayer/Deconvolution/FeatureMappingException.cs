using System;

namespace EngineLayer.Deconvolution;

public class FeatureMappingException : MetaMorpheusException
{
    public FeatureMappingException(string message, Exception innerException = null) : base(message, innerException)
    {
    }
}
