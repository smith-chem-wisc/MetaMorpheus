using System;

namespace EngineLayer
{
    [Serializable]
    public class MetaMorpheusException : Exception
    {
        public MetaMorpheusException(string message) : base(message)
        {
        }
    }
}