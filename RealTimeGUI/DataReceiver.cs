using System;
using System.Threading;
using System.Collections.Generic;
using System.Windows;
using System.Runtime.Remoting.Messaging;

using Thermo.Interfaces.ExactiveAccess_V1;
using Thermo.Interfaces.InstrumentAccess_V1.MsScanContainer;
using IMsScan = Thermo.Interfaces.InstrumentAccess_V2.MsScanContainer.IMsScan;

//[assembly: log4net.Config.XmlConfigurator(Watch = true)]

namespace RealTimeGUI
{
	/// <summary>
	/// Show incoming data packets and signals of acquisition start, acquisition stop and each scan.
	/// </summary>
	public class DataReceiver
	{
        private static readonly log4net.ILog logD = log4net.LogManager.GetLogger(System.Reflection.MethodBase.GetCurrentMethod().DeclaringType);

        internal DataReceiver()
        {
            RTParameters = new RTParameters();
            ListScan = new List<IMsScan>();
            Notifications = new Notifications();
        }

        public IExactiveInstrumentAccess InstrumentAccess { get; set; }

        public IMsScanContainer ScanContainer { get; set; }

        public List<IMsScan> ListScan { get; set; }
        
        public RTParameters RTParameters { get; set; }

        public Notifications Notifications { get; set; }

        public static event EventHandler<NotificationEventArgs> DataReceiverNotificationEventHandler;

        internal void ReceiveData()
        {
            string x = "\n" + DateTime.Now + " Start receive scans on detector " + ScanContainer.DetectorClass + "..." + "CurrentThread:" + Thread.CurrentThread.ManagedThreadId.ToString();
            //string x = "Start receive scans on detector " + ScanContainer.DetectorClass + ".";
            //DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(x));
            //logD.Debug("Start receive scans on detector " + ScanContainer.DetectorClass);
            this.Notifications.Notification += x;

            ScanContainer.AcquisitionStreamOpening += Orbitrap_AcquisitionStreamOpening;
            ScanContainer.AcquisitionStreamClosing += Orbitrap_AcquisitionStreamClosing;
            ScanContainer.MsScanArrived += Orbitrap_MsScanArrived;
        }

        internal void StopReceiveData()
        {    
            ScanContainer.MsScanArrived -= Orbitrap_MsScanArrived;
            ScanContainer.AcquisitionStreamClosing -= Orbitrap_AcquisitionStreamClosing;
            ScanContainer.AcquisitionStreamOpening -= Orbitrap_AcquisitionStreamOpening;
            string x = "\n" + DateTime.Now + " Stop receive scans on detector " + ScanContainer.DetectorClass + "..." + "CurrentThread:" + Thread.CurrentThread.ManagedThreadId.ToString();
            //string x = "Stop receive scans on detector " + ScanContainer.DetectorClass + "...";
            //DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(x));
            //logD.Debug(" Stop receive scans on detector " + ScanContainer.DetectorClass);
            this.Notifications.Notification += x;
        }

        private void Orbitrap_MsScanArrived(object sender, MsScanEventArgs e)
		{
			// If examining code takes longer, in particular for some scans, it is wise
			// to use a processing queue in order to get the system as responsive as possible.

			using (IMsScan scan = (IMsScan) e.GetScan())	// caution! You must dispose this, or you block shared memory!
			{
                //Console.WriteLine("\n{0:HH:mm:ss,fff} scan with {1} centroids arrived", DateTime.Now, scan.CentroidCount);
                string x = "\n" + DateTime.Now + " " + Thread.CurrentThread.Name;
                //DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(x));
                //DataReceiverNotificationEventHandler?.BeginInvoke(this, new NotificationEventArgs(x), callBack, null);
                //ListScan.Add(scan);
                //Dump("Specific", scan.SpecificInformation);
                //logD.Debug("Get one scan");
                this.Notifications.Notification += x;

            }
		}

		private void Orbitrap_AcquisitionStreamClosing(object sender, EventArgs e)
		{
            string x = "\n" + DateTime.Now + "Acquisition stream closed (end of method)" + "\n";
            DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(x));
        }

		private void Orbitrap_AcquisitionStreamOpening(object sender, MsAcquisitionOpeningEventArgs e)
		{
            string x = "\n" + DateTime.Now + "Acquisition stream opens (start of method)" + "\n";
            DataReceiverNotificationEventHandler?.Invoke(this, new NotificationEventArgs(x));
        }

        private void Dump(string title, IInfoContainer container)
        {
            //Console.WriteLine(title);
            foreach (string key in container.Names)
            {
                string value;
                // the source has to be "Unknown" to match all sources. Otherwise, the source has to match.
                // Sources can address values appearing in Tune files, general settings, but also values from
                // status log or additional values to a scan.
                MsScanInformationSource source = MsScanInformationSource.Unknown;   // show everything
                try
                {
                    if (container.TryGetValue(key, out value, ref source))
                    {
                        string descr = source.ToString();
                        descr = descr.Substring(0, Math.Min(11, descr.Length));
                        //Console.WriteLine("   {0,-11} {1,-35} = {2}", descr, key, value);
                        logD.Debug(descr);
                    }
                }
                catch { /* up to and including 2.8SP1 there is a bug displaying items which are null and if Foundation 3.1SP4 is used with CommonCore */ }
            }
        }

        private void callBack(IAsyncResult asyncResult)
        {
            var syncResult = (AsyncResult)asyncResult;
            var invokedMethod = (EventHandler<NotificationEventArgs>)syncResult.AsyncDelegate;

            try
            {
                invokedMethod.EndInvoke(asyncResult);
            }
            catch (Exception)
            {

                throw;
            }
        }

        public void TestLog()
        {
            logD.Debug("Test log in DataReceiver");
        }


	}
}
