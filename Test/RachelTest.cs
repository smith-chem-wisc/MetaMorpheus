using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class RachelTest
    {
        
        // want a psm whose base sequence is not ambigous but full sequence is (ptm is not localized): make sure this does not make it in DB
       
       [Test]
        public static void TestPrunedGeneration()
        {
            //Create GPTMD Task
            //Create Search Task
            GptmdTask task1 = new GptmdTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false
                },
            };

            SearchTask task2 = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false
                },
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchTarget = true,
                    WritePrunedDatabase = true,
                    SearchType = SearchType.Classic
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
                ("task1", task1),
                ("task2", task2),
            };
            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, Environment.CurrentDirectory);
            engine.Run();
            string outputFolderInThisTest = MySetUpClass.outputFolder;
            string final = Path.Combine(MySetUpClass.outputFolder, "task2","DbForPrunedDbGPTMDproteinPruned.xml");
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var ok);
            //ensures that protein out put contins the correct number of proteins to match the folowing conditions. 
                // all proteins in DB have baseSequence!=null (not ambiguous)
                // all proteins that belong to a protein group are written to DB
            Assert.AreEqual(proteins.Count(),20);
            int totalNumberOfMods = 0;
            foreach (Protein p in proteins)
            {
                int numberOfMods = p.OneBasedPossibleLocalizedModifications.Count();
                totalNumberOfMods=totalNumberOfMods + numberOfMods;
            }
            //tests that modifications are being done correctly
            Assert.AreEqual(totalNumberOfMods, 0);

        }
    }
}     

