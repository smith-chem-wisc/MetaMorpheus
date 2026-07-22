using EngineLayer;
using System;

namespace EngineLayer.Deconvolution.FeatureFileMapping;

public class FeatureMappingException : MetaMorpheusException
{
    public FeatureMappingException(string message, Exception innerException = null) : base(message, innerException)
    {
    }
}
