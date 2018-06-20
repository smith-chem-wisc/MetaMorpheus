using System;

namespace RealTimeGUI
{
    public class NotificationEventArgs : EventArgs
    {
        #region Public Constructors

        public NotificationEventArgs(string notification)
        {
            this.Notification = notification;
        }

        #endregion Public Constructors

        #region Public Properties

        public string Notification { get; private set; }

        #endregion Public Properties
    }
}