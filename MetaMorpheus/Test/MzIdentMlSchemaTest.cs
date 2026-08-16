using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// mzIdentML 1.1 requires at least one PeptideHypothesis inside every
    /// ProteinDetectionHypothesis, and at least one ProteinDetectionHypothesis inside every
    /// ProteinAmbiguityGroup (both are minOccurs="1" in the HUPO-PSI schema). A file that breaks
    /// either is rejected outright by validating consumers such as ProteomeXchange, which is the
    /// mzid half of mzLib issue 259.
    ///
    /// The existing MzIdentMLWriter tests all pass an empty protein group list, so nothing has ever
    /// exercised the ProteinDetectionList branch.
    /// </summary>
    [TestFixture]
    public static class MzIdentMlSchemaTest
    {
        [Test]
        public static void ProteinGroupSupportedOnlyByAnAmbiguousPsmStillWritesValidMzid()
        {
            // A PSM whose sequence could not be resolved has FullSequence == null. The writer skips
            // those when building PeptideEvidence, but the protein groups they support are still
            // written, so the group's ProteinDetectionHypothesis has nothing to point at.
            var psm = AmbiguousPsm(out Protein protein);
            Assert.That(psm.FullSequence, Is.Null, "the PSM under test is meant to be ambiguous");

            var group = new ProteinGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods>(psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer)),
                new HashSet<IBioPolymerWithSetMods>());
            group.SequenceCoverageFraction.Add(0.5);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "ambiguous_only_group.mzID");
            Write(new List<SpectralMatch> { psm }, new List<ProteinGroup> { group }, path);

            try
            {
                AssertMinimumChildCounts(path);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public static void ProteinInAGroupWithNoPeptidesFromThisFileStillWritesValidMzid()
        {
            // The per-file shape ProteinGroup.ConstructSubsetProteinGroup produces: every protein of
            // the group is kept, but AllPeptides is narrowed to the peptides seen in this file. In a
            // multi-file search a protein whose evidence came from another file therefore has no
            // peptide to point at, and its ProteinDetectionHypothesis comes out empty.
            Protein withEvidence = new Protein("PEPTIDEK", "ACC1", databaseFilePath: "temp");
            Protein withoutEvidence = new Protein("ANOTHERPEPTIDEK", "ACC2", databaseFilePath: "temp");

            var peptide = withEvidence.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, Scan(), new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            var group = new ProteinGroup(
                new HashSet<IBioPolymer> { withEvidence, withoutEvidence },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods>());
            group.SequenceCoverageFraction.Add(0.5);
            group.SequenceCoverageFraction.Add(0.0);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "subset_group.mzID");
            Write(new List<SpectralMatch> { psm }, new List<ProteinGroup> { group }, path);

            try
            {
                AssertMinimumChildCounts(path);

                // ACC1 keeps its hypothesis; only the evidence-free ACC2 is left out.
                Assert.That(CountElements(path, "ProteinAmbiguityGroup"), Is.EqualTo(1));
                Assert.That(CountElements(path, "ProteinDetectionHypothesis"), Is.EqualTo(1));
                Assert.That(File.ReadAllText(path), Does.Contain("DBS_ACC1").And.Not.Contain("DBS_ACC2\""));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public static void ProteinGroupSupportedByAResolvedPsmWritesValidMzid()
        {
            // The control: the same shape of output, but with a PSM the writer does emit
            // PeptideEvidence for. Guards against a fix that drops every group indiscriminately.
            Protein protein = new Protein("PEPTIDEK", "ACC1", databaseFilePath: "temp");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, Scan(), new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            Assert.That(psm.FullSequence, Is.Not.Null);

            var group = new ProteinGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });
            group.SequenceCoverageFraction.Add(0.5);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "resolved_group.mzID");
            Write(new List<SpectralMatch> { psm }, new List<ProteinGroup> { group }, path);

            try
            {
                AssertMinimumChildCounts(path);

                // The group must actually be present, not silently dropped.
                Assert.That(CountElements(path, "ProteinDetectionHypothesis"), Is.EqualTo(1));
                Assert.That(CountElements(path, "ProteinAmbiguityGroup"), Is.EqualTo(1));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public static void SilacSearchWithAnAmbiguousPsmDoesNotThrow()
        {
            // The SILAC pre-filter reads p.FullSequence.Contains("|"), but an ambiguous PSM has a
            // null FullSequence -- the neighbouring p.BaseSequence != null check does not cover it.
            // Only reachable when SilacLabels is non-null, i.e. an actual SILAC search.
            var psm = AmbiguousPsm(out Protein protein);
            var labels = new List<SilacLabel> { new SilacLabel('K', 'a', "C{13}6", 6.0201290768) };

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "silac_ambiguous.mzID");

            try
            {
                Assert.That(() => Write(new List<SpectralMatch> { psm }, new List<ProteinGroup>(), path, labels),
                    Throws.Nothing);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        private static SpectralMatch AmbiguousPsm(out Protein protein)
        {
            protein = new Protein("PEPTIDEK", "ACC1", databaseFilePath: "temp");

            // A variable mod with two candidate sites on the same peptide gives two peptides of
            // identical mass and different FullSequence, which is what makes a PSM ambiguous.
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            var variableMod = new Modification(_originalId: "Oxidation", _modificationType: "Common Variable",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915);

            var candidates = protein
                .Digest(new DigestionParams(), new List<Modification>(), new List<Modification> { variableMod })
                .Where(p => p.AllModsOneIsNterminus.Count == 1)
                .Take(2)
                .ToList();
            Assert.That(candidates.Count, Is.EqualTo(2), "digestion did not produce two singly-modified candidates");

            var psm = new PeptideSpectralMatch(candidates[0], 0, 10, 0, Scan(), new CommonParameters(), new List<MatchedFragmentIon>());
            psm.AddOrReplace(candidates[1], 10, 0, true, new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            return psm;
        }

        private static Ms2ScanWithSpecificMass Scan()
        {
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true,
                Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1",
                double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            return new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());
        }

        // SearchParameters.SilacLabels is left null for a non-SILAC search, which is what the
        // default production path passes.
        private static void Write(List<SpectralMatch> psms, List<ProteinGroup> groups, string path,
            List<SilacLabel> silacLabels = null)
        {
            MzIdentMLWriter.WriteMzIdentMl(psms, groups, new List<Modification>(), new List<Modification>(),
                silacLabels, new List<DigestionAgent>(), new PpmTolerance(20), new PpmTolerance(20),
                0, path, true);
        }

        private static int CountElements(string path, string localName)
        {
            int count = 0;
            using var reader = XmlReader.Create(path);
            while (reader.Read())
            {
                if (reader.NodeType == XmlNodeType.Element && reader.LocalName == localName) count++;
            }
            return count;
        }

        /// <summary>
        /// Walks the ProteinDetectionList and asserts the two minOccurs="1" rules the writer can
        /// break. Checking the shape directly rather than validating against mzIdentML1.1.0.xsd
        /// keeps the schema out of the test data.
        /// </summary>
        private static void AssertMinimumChildCounts(string path)
        {
            var settings = new XmlReaderSettings { IgnoreWhitespace = true };
            using var reader = XmlReader.Create(path, settings);

            string currentPag = null, currentPdh = null;
            int pdhInPag = 0, peptideHypothesisInPdh = 0;

            while (reader.Read())
            {
                if (reader.NodeType == XmlNodeType.Element)
                {
                    switch (reader.LocalName)
                    {
                        case "ProteinAmbiguityGroup":
                            currentPag = reader.GetAttribute("id");
                            pdhInPag = 0;
                            break;
                        case "ProteinDetectionHypothesis":
                            currentPdh = reader.GetAttribute("id");
                            pdhInPag++;
                            peptideHypothesisInPdh = 0;
                            break;
                        case "PeptideHypothesis":
                            peptideHypothesisInPdh++;
                            break;
                    }
                    // A childless element is serialized self-closing, so it never yields an
                    // EndElement node for the checks below to fire on.
                    if (reader.IsEmptyElement)
                    {
                        if (reader.LocalName == "ProteinDetectionHypothesis")
                        {
                            Assert.Fail($"ProteinDetectionHypothesis '{currentPdh}' is empty; the schema requires at least one PeptideHypothesis");
                        }
                        if (reader.LocalName == "ProteinAmbiguityGroup")
                        {
                            Assert.Fail($"ProteinAmbiguityGroup '{currentPag}' is empty; the schema requires at least one ProteinDetectionHypothesis");
                        }
                    }
                }
                else if (reader.NodeType == XmlNodeType.EndElement)
                {
                    if (reader.LocalName == "ProteinDetectionHypothesis")
                    {
                        Assert.That(peptideHypothesisInPdh, Is.GreaterThanOrEqualTo(1),
                            $"ProteinDetectionHypothesis '{currentPdh}' has no PeptideHypothesis child, which mzIdentML 1.1 forbids");
                    }
                    else if (reader.LocalName == "ProteinAmbiguityGroup")
                    {
                        Assert.That(pdhInPag, Is.GreaterThanOrEqualTo(1),
                            $"ProteinAmbiguityGroup '{currentPag}' has no ProteinDetectionHypothesis child, which mzIdentML 1.1 forbids");
                    }
                }
            }
        }
    }
}
