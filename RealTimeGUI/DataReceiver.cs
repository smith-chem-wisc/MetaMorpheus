using System;
using System.Threading;
using System.Collections.Generic;
using System.Windows;

using Thermo.Interfaces.ExactiveAccess_V1;
using Thermo.Interfaces.InstrumentAccess_V1.MsScanContainer;
using IMsScan = Thermo.Interfaces.InstrumentAccess_V2.MsScanContainer.IMsScan;

namespace RealTimeGUI
{
	/// <summary>
	/// Show incoming data packets and signals of acquisition start, acquisition stop and each scan.
	/// </summary>
	class DataReceiver
	{
		internal DataReceiver() { }

        string TextTime { get; set; }
        List<IMsScan> ListScan { get; set; }

        public static event EventHandler<NotificationEventArgs> DataReceiverNotificationEventHandler;

        internal void DoJob()
		{
            ListScan = new List<IMsScan>();
			using (IExactiveInstrumentAccess instrument = Connection.GetFirstInstrument())
			{
				IMsScanContainer orbitrap = instrument.GetMsScanContainer(0);
                //Console.WriteLine("Waiting 60 seconds for scans on detector " + orbitrap.DetectorClass + "...");
                string x = "Waiting 60 seconds for scans on detector " + orbitrap.DetectorClass + "...";
                DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(x + "/n"));
                MainWindow.mainWindow.TbNotifications.AppendText(x + "/n");

                orbitrap.AcquisitionStreamOpening += Orbitrap_AcquisitionStreamOpening;
				orbitrap.AcquisitionStreamClosing += Orbitrap_AcquisitionStreamClosing;
				orbitrap.MsScanArrived += Orbitrap_MsScanArrived;
				Thread.CurrentThread.Join(10000);
				orbitrap.MsScanArrived -= Orbitrap_MsScanArrived;
				orbitrap.AcquisitionStreamClosing -= Orbitrap_AcquisitionStreamClosing;
				orbitrap.AcquisitionStreamOpening -= Orbitrap_AcquisitionStreamOpening;
			}
		}

		private void Orbitrap_MsScanArrived(object sender, MsScanEventArgs e)
		{
			// If examining code takes longer, in particular for some scans, it is wise
			// to use a processing queue in order to get the system as responsive as possible.

			using (IMsScan scan = (IMsScan) e.GetScan())	// caution! You must dispose this, or you block shared memory!
			{
                //Console.WriteLine("\n{0:HH:mm:ss,fff} scan with {1} centroids arrived", DateTime.Now, scan.CentroidCount);
                ListScan.Add(scan);
                DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(ListScan.Count.ToString()));
            }
		}

		private void Orbitrap_AcquisitionStreamClosing(object sender, EventArgs e)
		{
			//Console.WriteLine("\n{0:HH:mm:ss,fff} {1}", DateTime.Now, "Acquisition stream closed (end of method)");
            TextTime = "\n{0:HH:mm:ss,fff} {1}" + DateTime.Now + "Acquisition stream closed (end of method)" + "/n";
            MainWindow.mainWindow.TbNotifications.AppendText(TextTime + "/n");
            DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(TextTime));
        }

		private void Orbitrap_AcquisitionStreamOpening(object sender, MsAcquisitionOpeningEventArgs e)
		{
            //Console.WriteLine("\n{0:HH:mm:ss,fff} {1}", DateTime.Now, "Acquisition stream opens (start of method)");
            TextTime = "\n{0:HH:mm:ss,fff} {1}" + DateTime.Now + "Acquisition stream opens (start of method)" + "/n";
            MainWindow.mainWindow.TbNotifications.AppendText(TextTime + "/n");
            DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(TextTime));
        }
	}
}
