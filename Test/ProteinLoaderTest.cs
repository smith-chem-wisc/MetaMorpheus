using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class ProteinLoaderTest
    {
        [Test]
        public void TestProteinLoad()
        {
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta"));
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fa"));
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta.gz"));
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fa.gz"));
        }

        public class ProteinLoaderTask : MetaMorpheusTask
        {
            public ProteinLoaderTask(string x)
                : this()
            { }

            protected ProteinLoaderTask()
                : base(MyTask.Search)
            { }

            public void Run(string dbPath)
            {
                RunSpecific("", new List<DbForTask> { new DbForTask(dbPath, false) }, null, "", null);
            }

            protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
            {
                LoadProteins("", dbFilenameList, true, DecoyType.None, new List<string>(), new CommonParameters());
                return null;
            }
        }
    }
}