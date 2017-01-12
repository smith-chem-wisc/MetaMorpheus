using System;
using System.Runtime.Serialization;

namespace InternalLogicEngineLayer
{
    [Serializable]
    internal class EngineRunException : Exception
    {
        public EngineRunException()
        {
        }

        public EngineRunException(string message) : base(message)
        {
        }

        public EngineRunException(string message, Exception innerException) : base(message, innerException)
        {
        }

        protected EngineRunException(SerializationInfo info, StreamingContext context) : base(info, context)
        {
        }
    }
}