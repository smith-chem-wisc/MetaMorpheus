using System;

namespace EngineLayer
{
    [Serializable]
    public class MetaMorpheusException : Exception
    {

        #region Public Constructors

        public MetaMorpheusException(string message) : base(message)
        {
        }

        #endregion Public Constructors

    }
}