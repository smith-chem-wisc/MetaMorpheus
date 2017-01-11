using System;
using System.Runtime.Serialization;

namespace InternalLogic
{
    [Serializable]
    public class EngineValidationException : Exception
    {
        public EngineValidationException()
        {
        }

        public EngineValidationException(string message) : base(message)
        {
        }

        public EngineValidationException(string message, Exception innerException) : base(message, innerException)
        {
        }

        protected EngineValidationException(SerializationInfo info, StreamingContext context) : base(info, context)
        {
        }
    }
}