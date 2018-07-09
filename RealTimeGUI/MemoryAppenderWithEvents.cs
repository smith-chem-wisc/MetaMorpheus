using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using log4net;
using log4net.Core;
using log4net.Appender;

namespace RealTimeGUI
{
    public class MemoryAppenderWithEvents : MemoryAppender
    {
        public event EventHandler Updated;

        protected override void Append(LoggingEvent loggingEvent)
        {
            base.Append(loggingEvent);

            Updated?.Invoke(this, new EventArgs());
        }
    }
}
