using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Threading;

namespace MetaMorpheusGUI
{
    class SearchTimer
    {
        public static DispatcherTimer Timer;

        public static void Initialize()
        {
            Timer = new DispatcherTimer();
            Timer.Interval = TimeSpan.FromMilliseconds(300);
        }

        // starts timer to keep track of user keystrokes
        public static void Set()
        {
            // Reset the timer
            Timer.Stop();
            Timer.Start();
        }
    }
}
