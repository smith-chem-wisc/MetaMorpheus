using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RealTimeGUI
{
    public class Notifications : INotifyPropertyChanged
    {
        private string notification;
        public string Notification
        {
            get { return notification; }
            set
            {
                    notification = value;
                    NotifyPropertyChanged("Notification");

            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        public void NotifyPropertyChanged(string propName)
        {
            if (PropertyChanged != null)
                PropertyChanged(this, new PropertyChangedEventArgs(propName));
        }
    }
}
