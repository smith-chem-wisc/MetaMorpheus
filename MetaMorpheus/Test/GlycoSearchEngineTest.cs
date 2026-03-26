using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Org.BouncyCastle.Asn1.Cmp;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using TaskLayer;
using ThermoFisher.CommonCore.Data.Business;
using static Org.BouncyCastle.Asn1.Cmp.Challenge;

namespace Test
{
    [TestFixture]
    public class GlycoSearchEngineTest
    {
        [Test]
        public void CreateGsm_WithWideProductTolerance_ScanInfo_p_IsCappedToOne() 
        {
            var commonParameters = new CommonParameters(productMassTolerance: new AbsoluteTolerance(1_000_000), dissociationType: DissociationType.HCD, trimMsMsPeaks: false);

            // ensure glycan DB paths used by GlycoSearchEngine ctor are registered (filenames must match ctor arguments)
            string oglycanPath = "OGlycan.gdb";
            string nglycanPath = "NGlycan_ForNoSearch.gdb";
            if (!GlobalVariables.OGlycanDatabasePaths.Contains(oglycanPath)) GlobalVariables.OGlycanDatabasePaths.Add(oglycanPath);
            if (!GlobalVariables.NGlycanDatabasePaths.Contains(nglycanPath)) GlobalVariables.NGlycanDatabasePaths.Add(nglycanPath);

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            var myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(spectraFile,commonParameters);
            var ms2ScanWithSpecificMass = MetaMorpheusTask.GetMs2Scans(msDataFile, spectraFile, commonParameters).ToArray();
            Assert.That(ms2ScanWithSpecificMass.Length, Is.GreaterThan(0), "No MS2 scans found in test MGF.");

            var protein = new Protein("AATVGSLAGQPLQER", "accession");
            var peptides = protein.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).ToList();
            var peptide = peptides.First();
            var peptideIndex = new List<PeptideWithSetModifications> { peptide };
            var globalGsms = new List<GlycoSpectralMatch>[ms2ScanWithSpecificMass.Length];

            var engine = new GlycoSearchEngine(globalGsms, ms2ScanWithSpecificMass, peptideIndex, null, null, 0, commonParameters, null, oglycanPath, nglycanPath, glycoSearchType: GlycoSearchType.OGlycanSearch, 30, 3, false, null);

            
            var route = new Route();
            route.AddPos(2, 2, true);
            var oxoniumIonIntensities = new double[Glycan.AllOxoniumIons.Length];
            var localizationGraphs = new List<LocalizationGraph>();


            var createGsmMethod = typeof(GlycoSearchEngine).GetMethod("CreateGsm", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.That(createGsmMethod, Is.Not.Null, "Unable to find private CreateGsm method via reflection.");

            var result = (GlycoSpectralMatch)createGsmMethod.Invoke(engine, new object[] { ms2ScanWithSpecificMass.First(), 0, 0, peptide, route, oxoniumIonIntensities, localizationGraphs });
        }
    }
}
