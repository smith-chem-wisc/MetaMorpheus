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
using UsefulProteomicsDatabases;

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
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams());
            SearchParameters SearchParameters = new SearchParameters();
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("TaGe_SA_HeLa_04_subset_longestSeq.mzML", CommonParameters));

            var myMsDataFile = myFileManager.LoadFile(origDataFile, CommonParameters);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, DecoyType.Reverse, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, -1);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML", CommonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, CommonParameters, fsp, null, new List<string>(), SearchParameters.WriteSpectralLibrary).Run();
            SearchTask searchTask = new();
            Tolerance tolerance = new PpmTolerance(10.0); // I'm not sure if we actually need a wider tolerance, but
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(tolerance, MassDiffAcceptorType.Exact, searchTask.SearchParameters.CustomMdac);

            List<PeptideWithSetModifications> peptideWithSetModificationList = allPsmsArray.Where(p => p != null).SelectMany(p => p.BestMatchingPeptides).Select(p => p.Peptide).ToList();
            foreach (var pwsm in peptideWithSetModificationList)
            {
                if (pwsm.FullChemicalFormula.AtomCount == 0)
                {
                    int placeholder2 = 3;
                }
            }

            IsotopeAnalysis testIsotopeAnalysis = new(new List<string> { origDataFile }, allPsmsArray.ToList(), massDiffAcceptor);
            int placeholder = 0;

            
        }
    }
}
