using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// A peptide's modifications reach mzIdentML in position order, whatever order its modification dictionary was
    /// filled in.
    /// </summary>
    /// <remarks>
    /// That order depends on how the peptide was built. Digestion inserts modifications in the order its patterns
    /// place them; a peptide read back from the Modern search's peptide-index cache has them in position order. A
    /// Modern search uses the peptides it built for its first spectra file and reads the cache back for every later
    /// one, so before this the same search wrote its mzIdentML Modification elements in a different order depending on
    /// whether the cache already existed.
    /// </remarks>
    [TestFixture]
    public static class ModificationOutputOrderTests
    {
        // AllModsOneIsNterminus keys are one-based with the N-terminus as 1, so residue i of PEPTIDE has key i + 1
        private static readonly int[] ResiduePositions = { 2, 4, 7 };

        private static (PeptideWithSetModifications InPositionOrder, PeptideWithSetModifications InReverseOrder) SamePeptideFilledTwoWays()
        {
            var protein = new Protein("PEPTIDE", "ACCESSION", databaseFilePath: "temp");
            ModificationMotif.TryGetMotif("T", out var tMotif);
            ModificationMotif.TryGetMotif("E", out var eMotif);
            var onT = new Modification(_originalId: "ORDER_T", _modificationType: "test", _target: tMotif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            var onE = new Modification(_originalId: "ORDER_E", _modificationType: "test", _target: eMotif, _locationRestriction: "Anywhere.", _monoisotopicMass: 20);
            var mods = new List<(int Key, Modification Mod)> { (3, onE), (5, onT), (8, onE) };

            PeptideWithSetModifications Build(IEnumerable<(int Key, Modification Mod)> fillOrder) =>
                new(protein, new DigestionParams(), 1, 7, CleavageSpecificity.Full, "", 0, fillOrder.ToDictionary(p => p.Key, p => p.Mod), 0);

            var reversed = Build(Enumerable.Reverse(mods));
            Assert.That(reversed.AllModsOneIsNterminus.Keys, Is.EqualTo(new[] { 8, 5, 3 }), "the second peptide's dictionary really is filled out of position order");
            return (Build(mods), reversed);
        }

        private static SpectralMatch Psm(PeptideWithSetModifications peptide)
        {
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null,
                null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());
            var psm = new PeptideSpectralMatch(peptide, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            return psm;
        }

        private static XElement MzIdentMlPeptide(PeptideWithSetModifications peptide, string name)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(ModificationOutputOrderTests)}_{name}.mzID");
            try
            {
                MzIdentMLWriter.WriteMzIdentMl(new List<SpectralMatch> { Psm(peptide) }, new List<ProteinGroup>(), new List<Modification>(),
                    new List<Modification>(), new List<SilacLabel>(), new List<DigestionAgent>(), new PpmTolerance(20), new PpmTolerance(20),
                    0, path, true);
                return XDocument.Load(path).Descendants().Single(e => e.Name.LocalName == "Peptide");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public static void MzIdentMlWritesModificationsInPositionOrder()
        {
            var (inPositionOrder, inReverseOrder) = SamePeptideFilledTwoWays();

            XElement fromPositionOrder = MzIdentMlPeptide(inPositionOrder, "position");
            XElement fromReverseOrder = MzIdentMlPeptide(inReverseOrder, "reverse");

            List<int> locations = fromReverseOrder.Elements().Where(e => e.Name.LocalName == "Modification")
                .Select(e => int.Parse(e.Attribute("location").Value)).ToList();
            Assert.That(locations, Is.EqualTo(ResiduePositions));
            Assert.That(XNode.DeepEquals(fromPositionOrder, fromReverseOrder), "the same peptide gives the same mzIdentML Peptide element");
        }
    }
}
