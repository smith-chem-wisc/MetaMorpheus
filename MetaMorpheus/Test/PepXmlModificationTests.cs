using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// pepXML carries a peptide's modifications, as the format defines them: each modified residue's position and the
    /// mass of the modified residue, and the modified N- and C-terminus masses as attributes of their own.
    /// </summary>
    /// <remarks>
    /// The writer built the modifications but attached them only to peptides that had none (the test on the list's
    /// count was inverted), so no pepXML file MetaMorpheus wrote ever contained a modification, and every unmodified
    /// peptide carried an empty modification_info element.
    /// </remarks>
    [TestFixture]
    public static class PepXmlModificationTests
    {
        private static readonly double HydrogenMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double HydroxylMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass + HydrogenMass;

        private static Modification Mod(string id, string residue, double mass)
        {
            ModificationMotif.TryGetMotif(residue, out var motif);
            return new Modification(_originalId: id, _modificationType: "test", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: mass);
        }

        private static PeptideWithSetModifications Peptide(Dictionary<int, Modification> mods) =>
            new(new Protein("PEPTIDE", "ACCESSION", databaseFilePath: "temp"), new DigestionParams(), 1, 7, CleavageSpecificity.Full, "", 0, mods, 0);

        private static XElement SearchHit(PeptideWithSetModifications peptide, string name)
        {
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null,
                null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());
            var psm = new PeptideSpectralMatch(peptide, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(PepXmlModificationTests)}_{name}.pep.XML");
            try
            {
                PepXMLWriter.WritePepXml(new List<SpectralMatch> { psm }, new List<DbForTask> { new DbForTask("temp", false) },
                    new List<Modification>(), new List<Modification>(), new CommonParameters(), path);
                return XDocument.Load(path).Descendants().Single(e => e.Name.LocalName == "search_hit");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public static void PepXmlWritesResidueAndTerminalModifications()
        {
            // keys are one-based with the N terminus as 1: E2 -> 3, T4 -> 5, the C terminus of a 7-residue peptide -> 9.
            // Filled out of position order on purpose.
            var mods = new Dictionary<int, Modification>
            {
                { 9, Mod("C-term", "E", -0.984016) },
                { 5, Mod("on T", "T", 79.966331) },
                { 1, Mod("N-term", "P", 42.010565) },
                { 3, Mod("on E", "E", 14.01565) },
            };

            XElement modificationInfo = SearchHit(Peptide(mods), "modified").Elements().Single(e => e.Name.LocalName == "modification_info");

            var residues = modificationInfo.Elements().Where(e => e.Name.LocalName == "mod_aminoacid_mass")
                .Select(e => (Position: int.Parse(e.Attribute("position").Value), Mass: double.Parse(e.Attribute("mass").Value))).ToList();
            Assert.That(residues.Select(r => r.Position), Is.EqualTo(new[] { 2, 4 }), "residue modifications in position order");
            Assert.That(residues[0].Mass, Is.EqualTo(Residue.GetResidue('E').MonoisotopicMass + 14.01565).Within(1e-9));
            Assert.That(residues[1].Mass, Is.EqualTo(Residue.GetResidue('T').MonoisotopicMass + 79.966331).Within(1e-9));

            Assert.That(double.Parse(modificationInfo.Attribute("mod_nterm_mass").Value), Is.EqualTo(HydrogenMass + 42.010565).Within(1e-9));
            Assert.That(double.Parse(modificationInfo.Attribute("mod_cterm_mass").Value), Is.EqualTo(HydroxylMass - 0.984016).Within(1e-9));
        }

        [Test]
        public static void PepXmlWritesNoModificationInfoForAnUnmodifiedPeptide()
        {
            XElement searchHit = SearchHit(Peptide(new Dictionary<int, Modification>()), "unmodified");
            Assert.That(searchHit.Elements().Any(e => e.Name.LocalName == "modification_info"), Is.False);
        }

        [Test]
        public static void PepXmlWritesOnlyTheTerminusAttributesThePeptideHas()
        {
            var mods = new Dictionary<int, Modification> { { 5, Mod("on T", "T", 79.966331) } };
            XElement modificationInfo = SearchHit(Peptide(mods), "residueonly").Elements().Single(e => e.Name.LocalName == "modification_info");
            Assert.That(modificationInfo.Attribute("mod_nterm_mass"), Is.Null);
            Assert.That(modificationInfo.Attribute("mod_cterm_mass"), Is.Null);
            Assert.That(modificationInfo.Elements().Count(e => e.Name.LocalName == "mod_aminoacid_mass"), Is.EqualTo(1));
        }
    }
}
