using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using log4net;
using log4net.Core;
using log4net.Appender;
using System.Windows.Controls;
using System.Windows;

namespace RealTimeGUI
{
    public class TextBoxAppender : AppenderSkeleton
    {
        public RichTextBox RichTextBox { get; set; }

        protected override void Append(LoggingEvent loggingEvent)
        {
            Action operation = () => { this.RichTextBox.AppendText(RenderLoggingEvent(loggingEvent)); };
            this.RichTextBox.Dispatcher.Invoke(operation);
        }
    }
}
