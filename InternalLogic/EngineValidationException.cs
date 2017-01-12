using System;
using System.Runtime.Serialization;

namespace InternalLogicEngineLayer
{
    [Serializable]
    public class EngineValidationException : Exception
    {
        #region Public Constructors

        public EngineValidationException()
        {
        }

        public EngineValidationException(string message) : base(message)
        {
        }

        public EngineValidationException(string message, Exception innerException) : base(message, innerException)
        {
        }

        #endregion Public Constructors

        #region Protected Constructors

        protected EngineValidationException(SerializationInfo info, StreamingContext context) : base(info, context)
        {
        }

        #endregion Protected Constructors
    }
}