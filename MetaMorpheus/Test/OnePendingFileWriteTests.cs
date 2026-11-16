using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Calibration and averaging write each exported spectra file in the background while they work on the next.
    /// These tests pin what must not change: at most one write pending, files reported only once complete and in the
    /// order they were started, failures surfaced, and the calibration task's files written and listed as before.
    /// </summary>
    [TestFixture]
    public sealed class OnePendingFileWriteTests
    {
        [Test]
        public static void ASecondWriteWaitsForTheFirstAndReportsRunInOrder()
        {
            var writes = new OnePendingFileWrite();
            var firstMayFinish = new ManualResetEventSlim(false);
            var log = new List<string>();

            writes.Start(() => { firstMayFinish.Wait(); lock (log) log.Add("write 1"); }, () => log.Add("report 1"));
            Assert.That(writes.IsWritePending, Is.True);

            Task second = Task.Run(() => writes.Start(() => { lock (log) log.Add("write 2"); }, () => log.Add("report 2")));
            Assert.That(second.Wait(200), Is.False, "the second write started while the first was still pending");

            firstMayFinish.Set();
            second.Wait();
            writes.Finish();

            Assert.That(log, Is.EqualTo(new[] { "write 1", "report 1", "write 2", "report 2" }));
            Assert.That(writes.IsWritePending, Is.False);
        }

        [Test]
        public static void AFileIsReportedOnlyOnceItsWriteHasFinished()
        {
            var writes = new OnePendingFileWrite();
            int written = 0;
            int writtenWhenReported = -1;

            writes.Start(() => { Thread.Sleep(100); Interlocked.Exchange(ref written, 1); },
                () => writtenWhenReported = Volatile.Read(ref written));
            Assert.That(writtenWhenReported, Is.EqualTo(-1), "reported before Finish");

            writes.Finish();
            Assert.That(writtenWhenReported, Is.EqualTo(1));
        }

        [Test]
        public static void WithNoWriteTheReportRunsAtOnceAfterThePendingOne()
        {
            var writes = new OnePendingFileWrite();
            var log = new List<string>();

            writes.Start(() => Thread.Sleep(50), () => log.Add("report 1"));
            writes.Start(null, () => log.Add("report 2"));

            Assert.That(log, Is.EqualTo(new[] { "report 1", "report 2" }));
            Assert.That(writes.IsWritePending, Is.False);
        }

        [Test]
        public static void AFailedWriteIsRethrownAndItsFileNotReported()
        {
            var writes = new OnePendingFileWrite();
            bool reported = false;

            writes.Start(() => throw new IOException("disk full"), () => reported = true);

            Assert.That(() => writes.Finish(), Throws.TypeOf<IOException>().With.Message.EqualTo("disk full"));
            Assert.That(reported, Is.False);
            Assert.That(() => writes.Finish(), Throws.Nothing, "the failure was surfaced twice");
        }

        [Test]
        public static void AFailedWriteIsSurfacedByTheNextStart()
        {
            var writes = new OnePendingFileWrite();
            bool secondStarted = false;

            writes.Start(() => throw new IOException("disk full"), () => { });

            Assert.That(() => writes.Start(() => secondStarted = true, () => { }), Throws.TypeOf<IOException>());
            Assert.That(secondStarted, Is.False);
        }

        [Test]
        public static void AHandledFailureIsHandedOverAndTheFileStillReported()
        {
            var writes = new OnePendingFileWrite();
            Exception handled = null;
            bool reported = false;

            writes.Start(() => throw new IOException("disk full"), () => reported = true, e => handled = e);
            writes.Finish();

            Assert.That(handled, Is.TypeOf<IOException>());
            Assert.That(reported, Is.True);
        }

        [Test]
        public static void DisposeWaitsForThePendingWriteWithoutReportingOrRethrowing()
        {
            var writes = new OnePendingFileWrite();
            int finished = 0;
            bool reported = false;

            writes.Start(() => { Thread.Sleep(100); Interlocked.Exchange(ref finished, 1); throw new IOException("late"); },
                () => reported = true);
            Assert.That(() => writes.Dispose(), Throws.Nothing);

            Assert.That(Volatile.Read(ref finished), Is.EqualTo(1), "Dispose returned before the write finished");
            Assert.That(reported, Is.False);
        }

        /// <summary>
        /// Calibrates two files and checks, at the moment each calibrated mzML is reported written, that the file is
        /// complete - it reads back with every scan - and that the files are listed for the next task in input order.
        /// </summary>
        [Test]
        public static void CalibrationReportsEachCalibratedFileCompleteAndInOrder()
        {
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "OnePendingFileWriteCalibration");
            string outputFolder = Path.Combine(unitTestFolder, "output");
            if (Directory.Exists(unitTestFolder))
            {
                Directory.Delete(unitTestFolder, true);
            }
            Directory.CreateDirectory(outputFolder);

            var inputs = new List<string>();
            foreach (string name in new[] { "first.mzML", "second.mzML" })
            {
                string path = Path.Combine(unitTestFolder, name);
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML"), path, true);
                inputs.Add(path);
            }
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");

            var reportedMzmls = new List<string>();
            var incompleteWhenReported = new List<string>();
            EventHandler<SingleFileEventArgs> onWritten = (_, e) =>
            {
                if (!e.WrittenFile.EndsWith(".mzML", StringComparison.OrdinalIgnoreCase))
                {
                    return;
                }
                reportedMzmls.Add(e.WrittenFile);
                try
                {
                    MsDataFile written = MsDataFileReader.GetDataFile(e.WrittenFile).LoadAllStaticData();
                    if (written.NumSpectra == 0)
                    {
                        incompleteWhenReported.Add(e.WrittenFile);
                    }
                }
                catch (Exception)
                {
                    incompleteWhenReported.Add(e.WrittenFile);
                }
            };

            MetaMorpheusTask.FinishedWritingFileHandler += onWritten;
            try
            {
                var results = new CalibrationTask().RunTask(outputFolder, new List<DbForTask> { new DbForTask(database, false) }, inputs, "test");

                var expected = inputs.Select(p => Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(p) + CalibrationTask.CalibSuffix + ".mzML")).ToList();
                Assert.That(results.NewSpectra, Is.EqualTo(expected));
                Assert.That(reportedMzmls, Is.EqualTo(expected));
                Assert.That(incompleteWhenReported, Is.Empty);
                Assert.That(expected.Select(p => File.ReadAllText(p).TrimEnd()).All(t => t.EndsWith("</indexedmzML>")), Is.True);
            }
            finally
            {
                MetaMorpheusTask.FinishedWritingFileHandler -= onWritten;
                Directory.Delete(unitTestFolder, true);
            }
        }
    }
}
