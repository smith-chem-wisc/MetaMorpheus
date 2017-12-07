using System;

namespace EngineLayer
{
    [Serializable]//this should not be a problem anyway
    public class MetaMorpheusException : Exception
    {
        #region Public Constructors

        public MetaMorpheusException(string message) : base(message)
        {
        }

        #endregion Public Constructors
    }
}