using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using TaskLayer;
using EngineLayer.ClassicSearch;

namespace Test
{
    internal class IsotopeAnalysisTest
    {

        /// <summary>
        /// Test ensures peptide FDR is calculated and that it doesn't output PSM FDR results
        /// </summary>
        [Test]
        public static void PeptideFDRTest()
        {
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra2.mzml");
            if (!File.Exists(myFile2)) { File.Copy(myFile, myFile2); }
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPeptideFDR");
            DbForTask db = new DbForTask(myDatabase, true);
            Directory.CreateDirectory(folderPath);

            SearchTask searchTask = new SearchTask();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Exact, searchTask.SearchParameters.CustomMdac);

            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);
            MsDataFile myMsDataFile2 = myFileManager.LoadFile(myFile2, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("VIHDNFGIVEGLMTTVHAITATQK", "P16858") };

            string targetSpectralLibrary = Path.Combine(folderPath, @"P16858_target.msp");
            string decoySpectralLibrary = Path.Combine(folderPath, @"P16858_decoy.msp");

            List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            //SpectralLibrary sl = new(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, null, null, new List<string>(), writeSpectralLibrary).Run();

            IsotopeAnalysis testIsotopeAnalysis = new(new List<string> { myFile, myFile2 }, allPsmsArray.ToList(), massDiffAcceptor);


            // search something with multiple hits of the same peptide to see if peptide FDR is calculated at the end
            new SearchTask().RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile, myFile2 }, "normal");
            List<string> columns = null;
            int cumDecoys = 0;
            int cumTargets = 0;
            double finalQValue = 0;
            foreach (string line in File.ReadAllLines(Path.Combine(folderPath, @"AllPeptides.psmtsv")))
            {
                string[] lineline = line.Split('\t');
                if (line.StartsWith("File Name")) // header
                {
                    columns = lineline.ToList();
                }

                // since each PSM has a duplicate, these counts will be 1,3,5,7, etc. if peptide FDR isn't calculated
                // if peptide FDR is calculated, they will be 1,2,3,4, etc. as expected
                else if (lineline[columns.IndexOf("Decoy/Contaminant/Target")] == "D")
                {
                    Assert.AreEqual(++cumDecoys, int.Parse(lineline[columns.IndexOf("Cumulative Decoy")]));
                    finalQValue = double.Parse(lineline[columns.IndexOf("QValue")], CultureInfo.InvariantCulture);
                }
                else
                {
                    Assert.AreEqual(++cumTargets, int.Parse(lineline[columns.IndexOf("Cumulative Target")]));
                    finalQValue = double.Parse(lineline[columns.IndexOf("QValue")], CultureInfo.InvariantCulture);
                }
            }

            // test that the final q-value follows the (target / decoy) formula
            // intermediate q-values no longer always follow this formula, so I'm not testing them here
            Assert.AreEqual(cumDecoys / (double)cumTargets, finalQValue, 0.0001);
            Directory.Delete(folderPath, true);



            //var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            //string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            //MyFileManager myFileManager = new MyFileManager(true);
            //CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            //MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            //var variableModifications = new List<Modification>();
            //var fixedModifications = new List<Modification>();
            //var proteinList = new List<Protein> { new Protein("VIHDNFGIVEGLMTTVHAITATQK", "P16858") };

            //string targetSpectralLibrary = Path.Combine(testDir, @"P16858_target.msp");
            //string decoySpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp");

            //List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            //SpectralLibrary sl = new(specLibs);

            //var searchModes = new SinglePpmAroundZeroSearchMode(5);

            //Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            //bool writeSpectralLibrary = false;
            //new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
            //    proteinList, searchModes, commonParameters, null, sl, new List<string>(), writeSpectralLibrary).Run();
        }
    }
}
