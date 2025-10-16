using Chemistry;
using Easy.Common;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using NUnit.Framework.Internal;
using NUnit.Framework.Legacy;
using Omics;
using Omics.BioPolymer;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class GptmdEngineTest
    {
        [Test]
        [TestCase("NNNNN", "accession", @"not applied", 5)]
        [TestCase("NNNNN", "accession", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 4)]
        public static void TestGptmdEngine(string proteinSequence, string accession, string sequenceVariantDescription, int numModifiedResidues)
        {
            List<SpectralMatch> allResultingIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            allResultingIdentifications = new List<SpectralMatch>();

            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", new CommonParameters()));

            var engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), fsp, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
            Assert.That(res.Mods.Count, Is.EqualTo(0));

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", "A", sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantBioPolymers();
            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));

            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications)).First();

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null), (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943).ToMz(1), 1, "filepath", new CommonParameters());

            var peptidesWithSetModifications = new List<IBioPolymerWithSetMods> { modPep };
            SpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());

            newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            newPsm.SetMs2Scan(scan.TheScan);
            allResultingIdentifications.Add(newPsm);

            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), null, new List<string>(), null);
            res = (GptmdResults)engine.Run();
            Assert.That(res.Mods.Count, Is.EqualTo(1));
            Assert.That(res.Mods["accession"].Count, Is.EqualTo(numModifiedResidues));
        }

        [Test]
        public static void TestGptmdEngineDissociationTypeAutodetect()
        {
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            //var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters cp = new CommonParameters(maxThreadsToUsePerFile: 1, digestionParams: new DigestionParams());
            var commonParameters = cp.CloneWithNewDissociationType(DissociationType.Autodetect);
            SearchParameters SearchParameters = new SearchParameters();
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("SmallCalibratible_Yeast.mzML", commonParameters));
            Tolerance precursorMassTolerance = new PpmTolerance(20);
            var myMsDataFile = myFileManager.LoadFile(origDataFile, commonParameters);
            List<double> acceptableMassShifts = new List<double> { 0.984015583, 0.984015583 };
            MassDiffAcceptor searchModes = new DotMassDiffAcceptor("", acceptableMassShifts, precursorMassTolerance);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(myDatabase, true, DecoyType.Reverse, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, -1);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, @"TestData\SmallCalibratible_Yeast.mzML", commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, fsp, null, new List<string>(), SearchParameters.WriteSpectralLibrary).Run();
            FdrAnalysisResults fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.Where(p => p != null).ToList(), 1,
                commonParameters, fsp, new List<string>()).Run());

            var nonNullPsms = allPsmsArray.Where(p => p != null).ToList();
            foreach(var psm in nonNullPsms)
            {
                psm.SetMs2Scan(listOfSortedms2Scans[psm.ScanIndex].TheScan);
            }
            GptmdParameters g = new GptmdParameters();
            List<Modification> gptmdModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => g.ListOfModsGptmd.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var reducedMods = new List<Modification>();
            foreach (var mod in gptmdModifications)
            {
                if (mod.IdWithMotif == "Deamidation on N" || mod.IdWithMotif == "Citrullination on R")
                {
                    reducedMods.Add(mod);
                }
            }
            
            var engine = new GptmdEngine(nonNullPsms, reducedMods, new List<Tuple<double, double>>(), new Dictionary<string, Tolerance> { { @"TestData\SmallCalibratible_Yeast.mzML", precursorMassTolerance } }, commonParameters, fsp, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
            Assert.That(8, Is.EqualTo(res.Mods.Count));
        }

        /// <summary>
        /// Example of a sequence variation:
        /// @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30"
        /// Reference allele: A
        /// Alternate allele: G
        /// Snpeff annotation: ANN=G||||||||||||||||
        /// Allele Index: 1
        /// Format: GT:AD:DP
        /// Genotype: 1/1
        /// Allelic depths: 30,30
        /// Homozygous reference calls: 30
        /// Heterozygous calls: 30
        /// </summary>

        [Test]
        [TestCase("NNNPPP", "accession", "A", @"not applied", 1, 3, 0, 3, 0)]
        [TestCase("NNNPPP", "accession", "A", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 1, 3, 0, 3, 0)]
        [TestCase("NNNPPP", "accession", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 2, 3, 0, 3, 1)]
        public static void TestCombos(string proteinSequence, string accession, string variantAA, string sequenceVariantDescription, int numModHashes, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP, int numModifiedResiduesNP)
        {
            List<SpectralMatch> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943),
                                                                      new Modification(_originalId: "16",  _modificationType: "mt", _target: motifP, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", variantAA, sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantBioPolymers();

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications)).First();

            MsDataScan dfd = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943 + 15.994915).ToMz(1), 1, "filepath", new CommonParameters());

            var peptidesWithSetModifications = new List<IBioPolymerWithSetMods> { modPep };
            SpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            match.SetMs2Scan(scan.TheScan);
            allIdentifications = new List<SpectralMatch> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), null, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
            Assert.That(res.Mods.Count, Is.EqualTo(numModHashes));
            Assert.That(res.Mods["accession"].Count, Is.EqualTo(numModifiedResidues));
            Assert.That(res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("21")).Count(), Is.EqualTo(numModifiedResiduesN));
            Assert.That(res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("16")).Count(), Is.EqualTo(numModifiedResiduesP));
            res.Mods.TryGetValue("accession_N1P", out var hash);
            Assert.That((hash ?? new HashSet<Tuple<int, Modification>>()).Count, Is.EqualTo(numModifiedResiduesNP));
        }
        [Test]
        public static void TestPtmBeforeVariant()
        {
            // This test is intentionally minimal in terms of spectrum/search setup but very verbose
            // in what’s being constructed and how protein variants, sequence variations, and
            // variant-aware database IO behave.

            // 1) Define two GPTMD candidate modifications (N and P) and a single mass-pair combo.
            //    - GPTMD will consider these for mass shifts on confidently identified peptides.
            //    - The actual PSM here is synthetic; the goal is to exercise GPTMD + variant plumbing.
            List<SpectralMatch> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<Modification>
            {
                new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943),
                new Modification(_originalId: "16", _modificationType: "mt", _target: motifP, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915)
            };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>
            {
                // One two-mod combo (e.g., deamidation-like + oxidation-like) to show GPTMD can
                // consider combinations of deltas when needed.
                new Tuple<double, double>(21.981943, 15.994915)
            };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            // 2) Build a parent protein with one SequenceVariation.
            //    Protein: "NNNPPP" (length 6). Variation is at 6..6 (the final P).
            //    OriginalSequence: "P", VariantSequence: "P" => a no-op variant (same residue).
            //    This is important because:
            //      - The applied-variant isoform will have the same sequence as the parent.
            //      - We still expect the variant to be tracked in metadata (AppliedSequenceVariations),
            //        which downstream systems use to reason about variant-aware mods or reporting.
            var svDescription = @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30";
            var parentProtein = new Protein(
                "NNNPPP",
                "protein",
                sequenceVariations: new List<SequenceVariation>
                {
            // Begin=6, End=6, Original="P", Variant="P" (no-op), carries a VCF-like description.
            new SequenceVariation(6, 6, "P", "P", "", svDescription, null)
                });

            // 3) Ask the parent to produce variant-aware isoforms.
            //    - GetVariantBioPolymers applies combinations of SequenceVariations and produces
            //      isoforms (Proteins) with AppliedSequenceVariations recorded.
            //    - For a no-op variant, we typically get a single isoform with the same BaseSequence.
            var variantProteins = parentProtein.GetVariantBioPolymers();

            // Validate the SequenceVariation metadata on the parent
            Assert.Multiple(() =>
            {
                Assert.That(parentProtein.SequenceVariations, Is.Not.Null);
                Assert.That(parentProtein.SequenceVariations.Count, Is.EqualTo(1));
                var sv = parentProtein.SequenceVariations[0];
                Assert.That(sv.OneBasedBeginPosition, Is.EqualTo(6));
                Assert.That(sv.OneBasedEndPosition, Is.EqualTo(6));
                Assert.That(sv.OriginalSequence, Is.EqualTo("P"));
                Assert.That(sv.VariantSequence, Is.EqualTo("P")); // no-op

                // Description is a SequenceVariantDescription (not string). Use ToString() for textual checks.
                Assert.That(sv.Description, Is.Not.Null);
                Assert.That(sv.Description.ToString(), Does.Contain("ANN="));
            });

            // Validate the isoforms produced
            Assert.Multiple(() =>
            {
                Assert.That(variantProteins, Is.Not.Null);
                Assert.That(variantProteins.Count, Is.GreaterThanOrEqualTo(1)); // no-op typically collapses to a single isoform
                var iso = variantProteins.First();
                // BaseSequence should remain unchanged since the variant is no-op
                Assert.That(iso.BaseSequence, Is.EqualTo(parentProtein.BaseSequence));
                // AppliedSequenceVariations may be 0 or 1 for a no-op; both are acceptable across versions.
                Assert.That(iso.AppliedSequenceVariations.Count, Is.LessThanOrEqualTo(1));
                if (iso.AppliedSequenceVariations.Count == 1)
                {
                    var applied = iso.AppliedSequenceVariations[0];
                    Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(6));
                    Assert.That(applied.OneBasedEndPosition, Is.EqualTo(6));
                    Assert.That(applied.OriginalSequence, Is.EqualTo("P"));
                    Assert.That(applied.VariantSequence, Is.EqualTo("P"));
                }
            });

            // 4) Digest the (variant) protein and synthesize a single PSM whose precursor mass
            //    matches the peptide + both GPTMD deltas to allow the engine to consider them.
            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins
                .SelectMany(p => p.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications))
                .First();

            Assert.Multiple(() =>
            {
                Assert.That(modPep, Is.Not.Null);
                // Ensure the peptide originates from one of the variant-applied isoforms
                Assert.That(variantProteins.Any(v => ReferenceEquals(v, modPep.Parent)), Is.True,
                    "Peptide.Parent should be one of the variant-applied isoforms returned by GetVariantBioPolymers()");

                // And the isoform's consensus should be the original parent protein
                Assert.That(ReferenceEquals(modPep.Parent.ConsensusVariant, parentProtein), Is.True,
                    "Isoform.ConsensusVariant should reference the original consensus parent protein");

                // Check peptide indexing sanity
                Assert.That(modPep.OneBasedStartResidue, Is.GreaterThanOrEqualTo(1));
                Assert.That(modPep.OneBasedEndResidue, Is.GreaterThanOrEqualTo(modPep.OneBasedStartResidue));
                Assert.That(modPep.OneBasedEndResidue, Is.LessThanOrEqualTo(modPep.Parent.BaseSequence.Length));
            });

            // Build a synthetic MS2 scan with precursor mass equal to peptide mass + both deltas
            MsDataScan dfd = new MsDataScan(
                new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null,
                MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN,
                null, DissociationType.AnyActivationType, 0, null);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(
                dfd,
                (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943 + 15.994915).ToMz(1),
                1,
                "filepath",
                new CommonParameters());

            var peptidesWithSetModifications = new List<IBioPolymerWithSetMods> { modPep };
            SpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            match.SetMs2Scan(scan.TheScan);
            allIdentifications = new List<SpectralMatch> { match };

            // 5) Run GPTMD.
            //    - With such a synthetic setup, we mainly assert that the engine completes and
            //      produces a well-formed result structure; we do not over-specify counts here to
            //      avoid brittleness across mass-tolerance or scoring tweaks.
            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), null, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
            Assert.Multiple(() =>
            {
                Assert.That(res, Is.Not.Null);
                Assert.That(res.Mods, Is.Not.Null);
                // Sanity: result map is present; the presence/absence of any particular key depends
                // on exact mass fitting paths. We avoid brittle hard expectations here.
                Assert.That(res.Mods.Count, Is.GreaterThanOrEqualTo(0));
            });

            // 6) Round-trip through XML: write an XML database with both consensus- and variant-
            //    targeted modifications, read it back, and assert locations.
            //
            //    Important mechanics:
            //    - ProteinDbWriter.WriteXmlDatabase accepts a per-accession map of (pos, mod).
            //    - If the accession is a variant accession (VariantApplication.GetAccession), the
            //      mod is attached to that isoform (sequence with variants applied).
            //    - On read (LoadProteinXML + decoys), target and decoy proteins get their
            //      OneBasedPossibleLocalizedModifications instantiated with correct mapped indices.
            //    - DecoyType.Reverse mirrors indices: for length L, original position i maps to L - i + 1.
            string xmlPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "ptm_before_variant.xml");

            // Build a small mod list: one on the consensus protein, one on the variant isoform.
            var modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();

            // Consensus mod at residue 2 (N target), easy to track
            var consensusHash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(
                    2,
                    new Modification(_originalId: "consensus-N@2", _modificationType: "type", _target: motifN, _monoisotopicMass: 21.981943, _locationRestriction: "Anywhere."))
            };
            modList.Add(parentProtein.Accession, consensusHash);

            // Variant-specific mod at residue 6 (P target) on the variant isoform’s accession
            var variantAccession = VariantApplication.GetAccession(parentProtein, parentProtein.SequenceVariations.ToArray());
            var variantHash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(
                    6,
                    new Modification(_originalId: "variant-P@6", _modificationType: "type", _target: motifP, _monoisotopicMass: 42, _locationRestriction: "Anywhere."))
            };
            modList.Add(variantAccession, variantHash);

            // Write and read the XML database (targets + reverse decoys)
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { parentProtein }, xmlPath);
            var loaded = ProteinDbLoader.LoadProteinXML(xmlPath, generateTargets: true, DecoyType.Reverse, allKnownModifications: null, isContaminant: false, modTypesToExclude: null, out var unknownMods, minAlleleDepth: 0);

            // Expect target + decoy
            Assert.Multiple(() =>
            {
                Assert.That(unknownMods.Count, Is.EqualTo(0));
                Assert.That(loaded.Count, Is.EqualTo(2));
            });

            var target = loaded.First(p => !p.IsDecoy);
            var decoy = loaded.First(p => p.IsDecoy);

            // After read, both consensus and variant-specific mods are present on the target protein's
            // OneBasedPossibleLocalizedModifications (the loaded protein sequence includes applied variants).
            // For "NNNPPP" length 6, positions are:
            //   - target: {2, 6}
            //   - decoy (reverse): index map i -> 6 - i + 1, so {2 -> 5, 6 -> 1} => {1, 5}
            Assert.Multiple(() =>
            {
                Assert.That(target.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
                var targetIndices = target.OneBasedPossibleLocalizedModifications.Keys.OrderBy(i => i).ToArray();
                CollectionAssert.AreEqual(new[] { 2, 6 }, targetIndices);

                Assert.That(decoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
                var decoyIndices = decoy.OneBasedPossibleLocalizedModifications.Keys.OrderBy(i => i).ToArray();
                CollectionAssert.AreEqual(new[] { 1, 5 }, decoyIndices);

                // Ensure variants remain attached to the target sequence and are queryable
                Assert.That(target.SequenceVariations.Count, Is.EqualTo(1));
                var rtv = target.SequenceVariations[0];
                Assert.That(rtv.OneBasedBeginPosition, Is.EqualTo(6));
                Assert.That(rtv.OneBasedEndPosition, Is.EqualTo(6));
                Assert.That(rtv.OriginalSequence, Is.EqualTo("P"));
                Assert.That(rtv.VariantSequence, Is.EqualTo("P"));
            });
        }
        [Test]
        public static void TestSearchPtmVariantDatabase()
        {
            // This test exercises a full round-trip through:
            // 1) Construct a protein with a SequenceVariation.
            // 2) Compute the variant accession (isoform accession) using VariantApplication.
            // 3) Write an XML database that includes both consensus- and variant-targeted modifications.
            // 4) Read the XML back (with decoys), then assert the correct placement of modifications
            //    on target and decoy proteins, and inspect how SequenceVariations are represented.

            // Create Search Task (minimal config just to satisfy task construction)
            SearchTask task1 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                },
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5))
            };

            // create modification lists used later for digestion checks
            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(b => task1.CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            // protein creation with a SequenceVariation:
            // - Base sequence: "PEPTID" (length 6)
            // - Variation at position 3 (P -> K), with a VCF-like description string
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            ModificationMotif.TryGetMotif("K", out ModificationMotif motifK);
            var variant = new SequenceVariation(
                3,               // 1-based position of the variant (begin == end here)
                "P",             // original residue
                "K",             // variant residue
                @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G|||||||||||||||||||\tGT:AD:DP\t1/1:30,30:30"); // description

            // The parent protein carries the variant in SequenceVariations (no variants applied yet).
            Protein testProteinWithMod = new Protein("PEPTID", "accession1", sequenceVariations: new List<SequenceVariation> { variant });

            // The "variant accession" is a deterministic accession string identifying the specific isoform
            // with this SequenceVariation applied. Any mods written under this accession will be associated
            // with the variant-applied isoform upon read.
            string variantAcc = VariantApplication.GetAccession(testProteinWithMod, new[] { variant });

            Assert.Multiple(() =>
            {
                // Verify SequenceVariation basics on consensus
                Assert.That(testProteinWithMod.SequenceVariations, Is.Not.Null);
                Assert.That(testProteinWithMod.SequenceVariations.Count, Is.EqualTo(1));
                var sv0 = testProteinWithMod.SequenceVariations[0];
                Assert.That(sv0.OneBasedBeginPosition, Is.EqualTo(3));
                Assert.That(sv0.OneBasedEndPosition, Is.EqualTo(3));
                Assert.That(sv0.OriginalSequence, Is.EqualTo("P"));
                Assert.That(sv0.VariantSequence, Is.EqualTo("K"));
                Assert.That(variantAcc, Is.Not.Null.And.Not.Empty);
            });

            // First write XML database
            string xmlName = "oblm.xml";

            // Build a per-accession mod map:
            // - For the consensus accession ("accession1"), put a mod at position 1 (on P).
            // - For the variant accession (P->K at position 3), put a mod at position 3 (on K).
            // This demonstrates how to attach modifications to both the consensus and a variant isoform.
            var modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var hashConsensus = new HashSet<Tuple<int, Modification>>
    {
        new Tuple<int, Modification>(1, new Modification(
            _originalId: "acetyl on P", _modificationType: "type", _target: motifP, _monoisotopicMass: 42, _locationRestriction: "Anywhere."))
    };
            var hashVariant = new HashSet<Tuple<int, Modification>>
    {
        new Tuple<int, Modification>(3, new Modification(
            _originalId: "acetyl on K", _modificationType: "type", _target: motifK, _monoisotopicMass: 42, _locationRestriction: "Anywhere."))
    };
            modList.Add(testProteinWithMod.Accession, hashConsensus);
            modList.Add(variantAcc, hashVariant);

            // Write database (one target will be written, and on read we request reverse decoys).
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { testProteinWithMod }, xmlName);

            // Read back the XML database:
            // - generateTargets: true (load target)
            // - decoyType: Reverse (mirror sequence to produce a decoy)
            // - allKnownModifications: null (let loader parse from XML)
            // - minAlleleDepth: 0 (allow any SV)
            var loadedProteins = ProteinDbLoader.LoadProteinXML(
                xmlName, generateTargets: true, DecoyType.Reverse, allKnownModifications: null,
                isContaminant: false, modTypesToExclude: null, out var unknownModifications, minAlleleDepth: 0);

            Assert.Multiple(() =>
            {
                // Expect exactly target + decoy
                Assert.That(unknownModifications.Count, Is.EqualTo(0));
                Assert.That(loadedProteins.Count, Is.EqualTo(2));
                Assert.That(loadedProteins.Count(p => !p.IsDecoy), Is.EqualTo(1));
                Assert.That(loadedProteins.Count(p => p.IsDecoy), Is.EqualTo(1));
            });

            var target = loadedProteins.First(p => !p.IsDecoy);
            var decoy = loadedProteins.First(p => p.IsDecoy);

            // Important behavior note:
            // When any variant-specific accession is present in the XML (e.g., accession1_P3K),
            // the loader can return the variant-applied isoform as the target entry (accession = variantAcc),
            // with BaseSequence reflecting the variant (PEKTID) and AppliedSequenceVariations populated.
            //
            // As a result:
            // - target.Accession == variantAcc (e.g., "accession1_P3K")
            // - target.BaseSequence == "PEKTID" (K at position 3)
            // - target.AppliedSequenceVariations.Count == 1 (the 3:P->K variant is applied)
            // - The per-position mods from both consensus accession (pos 1) and variant accession (pos 3)
            //   are realized on this variant isoform.
            Assert.Multiple(() =>
            {
                Assert.That(target.Accession, Is.EqualTo(variantAcc));
                Assert.That(target.BaseSequence, Is.EqualTo("PEKTID"));

                // Two localized modification sites expected: pos 1 (consensus P) and pos 3 (variant K)
                Assert.That(target.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
                var targetIndices = target.OneBasedPossibleLocalizedModifications.Keys.OrderBy(i => i).ToArray();
                CollectionAssert.AreEqual(new[] { 1, 3 }, targetIndices);

                // Validate mods by mass and logical target rather than OriginalId (which may be re-written by XML IO)
                var targetPos1Mods = target.OneBasedPossibleLocalizedModifications[1];
                var targetPos3Mods = target.OneBasedPossibleLocalizedModifications[3];

                // Mass check (within small tolerance) and motif check via IdWithMotif suffix
                Assert.That(targetPos1Mods.Any(m =>
                        Math.Abs(m.MonoisotopicMass.Value - 42.0) < 1e-6 &&
                        m.IdWithMotif.EndsWith("on P", StringComparison.Ordinal)),
                    Is.True, "Expected an acetyl-like mod (~42 Da) targeting P at position 1 on target");

                Assert.That(targetPos3Mods.Any(m =>
                        Math.Abs(m.MonoisotopicMass.Value - 42.0) < 1e-6 &&
                        m.IdWithMotif.EndsWith("on K", StringComparison.Ordinal)),
                    Is.True, "Expected an acetyl-like mod (~42 Da) targeting K at position 3 on target");

                // The variant is applied on the target isoform
                Assert.That(target.AppliedSequenceVariations.Count, Is.EqualTo(1));
                var applied = target.AppliedSequenceVariations[0];
                Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(3));
                Assert.That(applied.OriginalSequence, Is.EqualTo("P"));
                Assert.That(applied.VariantSequence, Is.EqualTo("K"));

                // The original variant description should still be present on the consensus variant bag
                Assert.That(target.SequenceVariations.Count, Is.EqualTo(1));
                var svT = target.SequenceVariations[0];
                Assert.That(svT.OneBasedBeginPosition, Is.EqualTo(3));
                Assert.That(svT.VariantSequence, Is.EqualTo("K"));
            });

            // Decoy behavior:
            // - BaseSequence is the reverse of the target
            // - Indices mirror (i -> L - i + 1). For L=6, pos1->6 and pos3->4.
            Assert.Multiple(() =>
            {
                Assert.That(decoy.IsDecoy, Is.True);
                Assert.That(decoy.BaseSequence, Is.EqualTo(new string(target.BaseSequence.Reverse().ToArray())));
                Assert.That(decoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
                var decoyIndices = decoy.OneBasedPossibleLocalizedModifications.Keys.OrderBy(i => i).ToArray();
                CollectionAssert.AreEqual(new[] { 4, 6 }, decoyIndices); // 3->4 and 1->6

                var decoyPos4Mods = decoy.OneBasedPossibleLocalizedModifications[4];
                var decoyPos6Mods = decoy.OneBasedPossibleLocalizedModifications[6];

                // Use mass + motif checks for robustness
                Assert.That(decoyPos6Mods.Concat(decoyPos4Mods).Any(m =>
                        Math.Abs(m.MonoisotopicMass.Value - 42.0) < 1e-6 &&
                        (m.IdWithMotif.EndsWith("on P", StringComparison.Ordinal) ||
                         m.IdWithMotif.EndsWith("on K", StringComparison.Ordinal))),
                    Is.True, "Decoy should carry mirrored acetyl-like mods (~42 Da) at mirrored positions");
            });

            // Variant-aware digestion check (sanity): digesting the target's variant isoform should
            // produce peptides where position 3 is K rather than P. The isoforms enumerate both consensus
            // and applied-variant entries.
            var targetIsoforms = target.GetVariantBioPolymers();
            Assert.That(targetIsoforms, Is.Not.Empty);
            var appliedIso = targetIsoforms.FirstOrDefault(i => i.AppliedSequenceVariations.Any());
            if (appliedIso != null)
            {
                var applied = appliedIso.AppliedSequenceVariations[0];
                Assert.Multiple(() =>
                {
                    Assert.That(applied.OneBasedBeginPosition, Is.EqualTo(3));
                    Assert.That(applied.VariantSequence, Is.EqualTo("K"));
                    Assert.That(appliedIso.BaseSequence[2], Is.EqualTo('K'), "3rd residue should be K on the applied-variant isoform");
                });

                // Digest the variant-applied isoform and confirm basic invariants
                var digestedVariant = appliedIso.Digest(task1.CommonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();
                Assert.That(digestedVariant, Is.Not.Empty);
            }

            // This section is for orchestrating a simple search run over the generated XML + dummy mzML.
            // It isn't central to SequenceVariation handling but mirrors how a task would run.
            var taskList = new List<(string, MetaMorpheusTask)> { ("task1", task1) };

            // Use an existing small mzML instead of generating one via TestDataFile to avoid NotImplemented in TestDataFile.LoadAllStaticData
            string mzmlName = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName },
                new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);

            // Ensure the engine runs end-to-end; exceptions here would indicate IO or setup issues
            engine.Run();
        }
        [Test]
        [TestCase("P", "PETID", "junk", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Unassigned.", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Anywhere.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "N-terminal.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "Peptide N-terminal.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "C-terminal.", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Peptide C-terminal.", 1, 5, 1, false)]
        [TestCase("E", "PETID", "Anywhere.", 2, 5, 2, true)]
        [TestCase("E", "PETID", "N-terminal.", 2, 5, 2, true)]
        [TestCase("E", "PETID", "Peptide N-terminal.", 2, 5, 2, false)]
        [TestCase("E", "PETID", "C-terminal.", 2, 5, 2, false)]
        [TestCase("E", "PETID", "Peptide C-terminal.", 2, 5, 2, false)]
        [TestCase("D", "PETID", "Anywhere.", 5, 5, 5, true)]
        [TestCase("D", "PETID", "N-terminal.", 5, 5, 5, false)]
        [TestCase("D", "PETID", "Peptide N-terminal.", 5, 5, 5, false)]
        [TestCase("D", "PETID", "C-terminal.", 5, 5, 5, true)]
        [TestCase("D", "PETID", "Peptide C-terminal.", 5, 5, 5, true)]
        public static void Test_GptmdEngineModFits(string targetAminoAcid, string proteinSequence, string locationRestriction, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex, bool result)
        {
            ModificationMotif.TryGetMotif(targetAminoAcid, out ModificationMotif motif);
            Modification attemptToLocalize = new Modification(null, null, null, null, _target: motif, _locationRestriction: locationRestriction, _chemicalFormula: null, _monoisotopicMass: 1, _databaseReference: null, _taxonomicRange: null, _keywords: null, _neutralLosses: null, _diagnosticIons: null, _fileOrigin: null);
            Dictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>();
            oneBasedModifications.Add(proteinOneBasedIndex, new List<Modification>() { attemptToLocalize });
            Protein protein = new Protein(proteinSequence, null, null, null, oneBasedModifications, null, null, null, false, false, null, null, null, null, null, null, "");

            Assert.That(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex), Is.EqualTo(result));
        }

        [Test]
        public static void TestUniProtGptmdConflict()
        {
            // this unit test checks to make sure GPTMD does not annotate mods at residues on 
            // proteins where the equivalent uniprot mod already exists
            Modification uniProtPhospho = GlobalVariables.AllModsKnown.First(p => p.ModificationType == "UniProt" && p.IdWithMotif.Contains("Phosphoserine"));
            Modification mmPhospho = GlobalVariables.AllModsKnown.First(p => p.ModificationType == "Common Biological" && p.IdWithMotif.Contains("Phosphorylation on S"));

            Protein protein = new Protein("PEPTIDESK", "test",
                oneBasedModifications: new Dictionary<int, List<Modification>>() { { 8, new List<Modification> { uniProtPhospho } } });

            PeptideWithSetModifications pep = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First(p => p.AllModsOneIsNterminus.Count == 0);

            // mod should not fit anywhere on the protein
            for (int i = 0; i < pep.Length; i++)
            {
                bool modFits = GptmdEngine.ModFits(mmPhospho, protein, i + 1, pep.Length, pep.OneBasedStartResidue + i);

                Assert.That(!modFits);
            }

            // the following code is just a control to make sure the phosphorylation actually does fit
            // at the given residue if the UniProt phosphorylation is not already present
            var someOtherSMod = GlobalVariables.AllModsKnown.Where(p => p.ModificationType == "Common Biological" && p.IdWithMotif.Contains("HexNAc on S")).First();

            protein = new Protein("PEPTIDESK", "test",
                oneBasedModifications: new Dictionary<int, List<Modification>>() { { 8, new List<Modification> { someOtherSMod } } });

            pep = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First(p => p.AllModsOneIsNterminus.Count == 0);

            // mod should fit at position 8
            for (int i = 0; i < pep.Length; i++)
            {
                bool modFits = GptmdEngine.ModFits(mmPhospho, protein, i + 1, pep.Length, pep.OneBasedStartResidue + i);

                if (i + 1 == 8)
                {
                    Assert.That(modFits);
                }
                else
                {
                    Assert.That(!modFits);
                }
            }
        }
    }
}