using System;
using System.Runtime.Serialization;

namespace InternalLogicEngineLayer
{
    [Serializable]
    internal class EngineRunException : Exception
    {
        #region Public Constructors

        public EngineRunException()
        {
        }

        public EngineRunException(string message) : base(message)
        {
        }

        public EngineRunException(string message, Exception innerException) : base(message, innerException)
        {
        }

        #endregion Public Constructors

        #region Protected Constructors

        protected EngineRunException(SerializationInfo info, StreamingContext context) : base(info, context)
        {
        }

        #endregion Protected Constructors
    }
}