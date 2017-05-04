using IO.MzML;
using MassSpectrometry;
using System;
using System.Linq;

namespace Test
{
    public class FakeMsDataFile : MsDataFile<IMzmlScan>, IMsStaticDataFile<IMzmlScan>
    {

        #region Public Constructors

        public FakeMsDataFile(IMzmlScan[] FakeScans) : base(FakeScans)
        {
            this.Scans = FakeScans;
        }

        #endregion Public Constructors

        #region Public Methods

        public override int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            int ok = Array.BinarySearch(Scans.Select(b => b.RetentionTime).ToArray(), retentionTime);
            if (ok < 0)
                ok = ~ok;
            return ok + 1;
        }

        public override IMzmlScan GetOneBasedScan(int scanNumber)
        {
            return Scans[scanNumber - 1];
        }

        #endregion Public Methods

    }
}