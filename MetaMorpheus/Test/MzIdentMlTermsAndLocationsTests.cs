using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Linq;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// What a repository reads from MetaMorpheus's mzIdentML: each controlled-vocabulary term tagged with its own
    /// vocabulary, a file format that matches the file, and locations that name the deposited files rather than the
    /// searcher's own folders.
    /// </summary>
    [TestFixture]
    public static class MzIdentMlTermsAndLocationsTests
    {
        // An accession's prefix names its vocabulary, and the cvRef must say the same.
        private static readonly Dictionary<string, string> CvRefOfPrefix = new()
        {
            { "MS", "PSI-MS" },
            { "UNIMOD", "UNIMOD" },
            { "MOD", "PSI-MOD" },
            { "UO", "UO" },
        };

        [Test]
        [TestCase("Unimod", "35", "UNIMOD:35", "UNIMOD", TestName = "Unimod reference")]
        [TestCase("PSI-MOD", "MOD:00046", "MOD:00046", "PSI-MOD", TestName = "PSI-MOD reference as UniProt writes it")]
        [TestCase("PSI-MOD", "00046", "MOD:00046", "PSI-MOD", TestName = "PSI-MOD reference without its prefix")]
        public static void ModificationTermCarriesItsOwnVocabulary(string database, string reference, string accession, string cvRef)
        {
            var mod = Mod("Oxidation", new Dictionary<string, IList<string>> { { database, new List<string> { reference } } });

            XDocument doc = WriteAndLoad(mod, "terms.mzID");

            // the peptide's own modification and the search's ModificationParams entry
            var terms = ModificationTerms(doc);
            Assert.That(terms, Has.Count.EqualTo(2));
            foreach (var term in terms)
            {
                Assert.That((string)term.Attribute("accession"), Is.EqualTo(accession));
                Assert.That((string)term.Attribute("cvRef"), Is.EqualTo(cvRef));
                Assert.That((string)term.Attribute("name"), Is.EqualTo("Oxidation"));
            }
        }

        [Test]
        public static void ResidReferenceGivesWayToThePsiModOne()
        {
            // UniProt's entries that have no Unimod reference carry RESID and PSI-MOD together, and mzIdentML takes its
            // modification terms from UNIMOD or PSI-MOD only. These are that UniProt entry's own references.
            var mod = Mod("(E)-2,3-didehydrobutyrine", new Dictionary<string, IList<string>>
            {
                { "RESID", new List<string> { "AA0547" } },
                { "PSI-MOD", new List<string> { "MOD:01470" } },
            });

            XDocument doc = WriteAndLoad(mod, "resid.mzID");

            var terms = ModificationTerms(doc);
            Assert.That(terms, Has.Count.EqualTo(2));
            foreach (var term in terms)
            {
                Assert.That((string)term.Attribute("accession"), Is.EqualTo("MOD:01470"));
                Assert.That((string)term.Attribute("cvRef"), Is.EqualTo("PSI-MOD"));
            }
            Assert.That(doc.ToString(), Does.Not.Contain("RESID:"));
        }

        [Test]
        [TestCase(true, TestName = "RESID reference only")]
        [TestCase(false, TestName = "no reference")]
        public static void ModificationWithNoUnimodOrPsiModTermIsAnUnknownModification(bool hasResid)
        {
            var references = hasResid
                ? new Dictionary<string, IList<string>> { { "RESID", new List<string> { "AA0547" } } }
                : null;
            var mod = Mod("Custom", references);

            XDocument doc = WriteAndLoad(mod, "unknown.mzID");

            var terms = ModificationTerms(doc);
            Assert.That(terms, Has.Count.EqualTo(2));
            foreach (var term in terms)
            {
                Assert.That((string)term.Attribute("accession"), Is.EqualTo("MS:1001460"));
                Assert.That((string)term.Attribute("name"), Is.EqualTo("unknown modification"));
                Assert.That((string)term.Attribute("cvRef"), Is.EqualTo("PSI-MS"));
                Assert.That((string)term.Attribute("value"), Is.EqualTo("Custom"));
            }
        }

        [Test]
        public static void EveryTermIsTaggedWithTheVocabularyItsAccessionNames()
        {
            // The whole document, not only the modifications: the units of the tolerances (UO), the software and the
            // scores (PSI-MS), and every modification branch at once.
            var mods = new List<Modification>
            {
                Mod("Oxidation", new Dictionary<string, IList<string>> { { "Unimod", new List<string> { "35" } } }),
                Mod("PsiModOnly", new Dictionary<string, IList<string>> { { "PSI-MOD", new List<string> { "MOD:00046" } } }),
                Mod("Custom", null),
            };
            XDocument doc = Write(mods, Psm(mods[0], "run1.mzML"), "every_term.mzID");

            var declared = doc.Descendants().Where(e => e.Name.LocalName == "cv")
                .Select(e => (string)e.Attribute("id")).ToHashSet();

            var tagged = doc.Descendants().Where(e => e.Name.LocalName == "cvParam")
                .SelectMany(e => new[] { (e.Attribute("accession"), e.Attribute("cvRef")), (e.Attribute("unitAccession"), e.Attribute("unitCvRef")) })
                .Where(pair => pair.Item1 != null)
                .ToList();
            Assert.That(tagged.Select(p => p.Item2?.Value), Does.Contain("UNIMOD").And.Contain("PSI-MOD").And.Contain("UO"),
                "the document must exercise every vocabulary, or this test proves nothing about them");

            foreach (var (accession, cvRef) in tagged)
            {
                string prefix = accession.Value.Split(':')[0];
                Assert.That(CvRefOfPrefix.ContainsKey(prefix), $"'{accession.Value}' is from no vocabulary mzIdentML uses");
                Assert.That(cvRef?.Value, Is.EqualTo(CvRefOfPrefix[prefix]), $"'{accession.Value}' is tagged with the wrong vocabulary");
                Assert.That(declared, Does.Contain(cvRef.Value), $"cvRef '{cvRef.Value}' is not declared in the cvList");
            }
        }

        [Test]
        public static void CvListNamesEachVocabularyAsItself()
        {
            XDocument doc = WriteAndLoad(Mod("Custom", null), "cvlist.mzID");
            var cvs = doc.Descendants().Where(e => e.Name.LocalName == "cv")
                .ToDictionary(e => (string)e.Attribute("id"));

            Assert.That((string)cvs["UNIMOD"].Attribute("fullName"), Is.EqualTo("UNIMOD"));
            Assert.That((string)cvs["UNIMOD"].Attribute("uri"), Does.Contain("unimod"));
            Assert.That((string)cvs["UO"].Attribute("fullName"), Is.EqualTo("UNIT-ONTOLOGY"));
            Assert.That((string)cvs["UO"].Attribute("uri"), Is.EqualTo("http://purl.obolibrary.org/obo/uo.obo"));
        }

        [Test]
        [TestCase("run1.raw", "MS:1000563", "MS:1000768")]
        [TestCase("run1.RAW", "MS:1000563", "MS:1000768")]
        [TestCase("run1.Raw", "MS:1000563", "MS:1000768")]
        [TestCase("run1.mzML", "MS:1000584", "MS:1001530")]
        public static void FileFormatFollowsTheExtensionInAnyCase(string fileName, string fileFormat, string spectrumIdFormat)
        {
            var mod = Mod("Custom", null);
            XDocument doc = Write(new List<Modification> { mod }, Psm(mod, fileName), "format.mzID");

            var spectraData = doc.Descendants().Single(e => e.Name.LocalName == "SpectraData");
            string AccessionUnder(string element) => (string)spectraData.Descendants()
                .Single(e => e.Name.LocalName == element).Descendants().Single(e => e.Name.LocalName == "cvParam").Attribute("accession");

            Assert.That(AccessionUnder("FileFormat"), Is.EqualTo(fileFormat));
            Assert.That(AccessionUnder("SpectrumIDFormat"), Is.EqualTo(spectrumIdFormat));
        }

        [Test]
        public static void LocationsNameTheFilesNotTheLocalFolders()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "MzIdentMlLocations", "someUser");
            string database = Path.Combine(folder, "human.fasta");
            string spectra = Path.Combine(folder, "run1-calib.mzML");

            var mod = Mod("Custom", null);
            XDocument doc = Write(new List<Modification> { mod }, Psm(mod, spectra, database), "locations.mzID");

            var spectraData = doc.Descendants().Single(e => e.Name.LocalName == "SpectraData");
            var searchDatabase = doc.Descendants().Single(e => e.Name.LocalName == "SearchDatabase");
            Assert.That((string)spectraData.Attribute("location"), Is.EqualTo("run1-calib.mzML"));
            Assert.That((string)spectraData.Attribute("name"), Is.EqualTo("run1-calib"));
            Assert.That((string)searchDatabase.Attribute("location"), Is.EqualTo("human.fasta"));

            // no attribute anywhere in the file carries the local folder
            Assert.That(doc.ToString(), Does.Not.Contain(folder).And.Not.Contain("someUser"));

            // the identifications still point at that spectra file and that database
            Assert.That(doc.Descendants().Where(e => e.Name.LocalName == "SpectrumIdentificationResult")
                .Select(e => (string)e.Attribute("spectraData_ref")), Is.EqualTo(new[] { (string)spectraData.Attribute("id") }));
            Assert.That(doc.Descendants().Where(e => e.Name.LocalName == "DBSequence")
                .Select(e => (string)e.Attribute("searchDatabase_ref")), Is.EqualTo(new[] { (string)searchDatabase.Attribute("id") }));
        }

        private static Modification Mod(string originalId, Dictionary<string, IList<string>> databaseReference)
        {
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            return new Modification(_originalId: originalId, _modificationType: "test", _target: motif,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915, _databaseReference: databaseReference);
        }

        // PEPTIDE with the modification on its T, as a fixed modification so every candidate carries it
        private static SpectralMatch Psm(Modification mod, string spectraFile, string databaseFile = "temp")
        {
            var protein = new Protein("PEPTIDE", "ACC1", databaseFilePath: databaseFile);
            var peptide = protein.Digest(new DigestionParams(), new List<Modification> { mod }, new List<Modification>()).First();
            Assert.That(peptide.NumMods, Is.EqualTo(1), "the test peptide must carry the modification");

            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, 2, true,
                Polarity.Positive, 1.0, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "controllerType=0 controllerNumber=1 scan=1",
                double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, spectraFile, new CommonParameters());

            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            return psm;
        }

        private static XDocument WriteAndLoad(Modification mod, string fileName) =>
            Write(new List<Modification> { mod }, Psm(mod, "run1.raw"), fileName);

        private static XDocument Write(List<Modification> fixedMods, SpectralMatch psm, string fileName)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, fileName);
            try
            {
                MzIdentMLWriter.WriteMzIdentMl(new List<SpectralMatch> { psm }, new List<ProteinGroup>(), new List<Modification>(),
                    fixedMods, null, new List<DigestionAgent>(), new PpmTolerance(20), new PpmTolerance(5), 2, path, false);
                return XDocument.Load(path);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        private static List<XElement> ModificationTerms(XDocument doc) => doc.Descendants()
            .Where(e => e.Name.LocalName == "Modification" || e.Name.LocalName == "SearchModification")
            .SelectMany(e => e.Elements().Where(c => c.Name.LocalName == "cvParam"))
            .ToList();
    }
}
