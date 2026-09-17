using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Readers;
using TaskLayer;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    /// <summary>
    /// Non-specific search over nucleic acids: the rnase singleN or singleC gives seed oligos fixed at one end, and the
    /// non-specific search engine trims the other end to the precursor mass, as it does for peptides.
    /// </summary>
    /// <remarks>
    /// The spectrum is the intact GUACUG sixmer (5'-OH, 3'-OH) that Classic and Modern search already find. Embedded in a
    /// longer molecule, it can only be found by trimming, and only if the trimmed end carries the right terminus: singleC
    /// leaves 5'-OH at the cut (matches), singleN leaves a 3'-phosphate (does not). Both facts are tested.
    /// </remarks>
    [TestFixture]
    public class TestRnaNonSpecificSearch
    {
        private static readonly string SixmerFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "GUACUG_NegativeMode_Sliced.mzML");
        private static readonly IHasChemicalFormula FivePrimeHydroxyl = ChemicalFormula.ParseFormula("O-3P-1");
        private static readonly IHasChemicalFormula ThreePrimePhosphate = ChemicalFormula.ParseFormula("H2O4P");
        private static readonly IHasChemicalFormula ThreePrimeHydroxyl = ChemicalFormula.ParseFormula("OH");

        private static CommonParameters RnaCommonParameters(IDigestionParams digestionParams, bool addCompIons = false) => new(
            dissociationType: DissociationType.CID,
            deconvolutionMaxAssumedChargeState: -20,
            deconvolutionIntensityRatio: 3,
            deconvolutionMassTolerance: new PpmTolerance(20),
            precursorMassTolerance: new PpmTolerance(10),
            productMassTolerance: new PpmTolerance(20),
            scoreCutoff: 5,
            totalPartitions: 1,
            maxThreadsToUsePerFile: 1,
            doPrecursorDeconvolution: true,
            useProvidedPrecursorInfo: false,
            addCompIons: addCompIons,
            digestionParams: digestionParams);

        private static RnaDigestionParams Seeds(string rnase, FragmentationTerminus terminus) =>
            new(rnase, minLength: 1, fragmentationTerminus: terminus);

        private static readonly MassDiffAcceptor TenPpm = new SinglePpmAroundZeroSearchMode(10);

        private static Ms2ScanWithSpecificMass[] SixmerScans(CommonParameters commonParameters) =>
            MetaMorpheusTask.GetMs2Scans(MsDataFileReader.GetDataFile(SixmerFilePath), SixmerFilePath, commonParameters)
                .OrderBy(b => b.PrecursorMass)
                .ToArray();

        /// <summary>Runs the non-specific engine over these molecules and returns its matches, by FDR category.</summary>
        private static SpectralMatch[][] RunNonSpecific(CommonParameters commonParameters, List<IBioPolymer> targets, List<Modification> variableMods = null)
        {
            variableMods ??= [];
            var scans = SixmerScans(commonParameters);
            var indexResults = (IndexingResults)new IndexingEngine(targets, variableMods, [], null, null, null, 0, DecoyType.None,
                commonParameters, null, 30000, true, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            int categories = Enum.GetValues<FdrCategory>().Length;
            var matches = new SpectralMatch[categories][];
            for (int i = 0; i < categories; i++)
            {
                matches[i] = new SpectralMatch[scans.Length];
            }
            var coisolationIndex = Enumerable.Range(0, scans.Length).Select(i => new List<int> { i }).ToArray();

            new NonSpecificEnzymeSearchEngine(matches, scans, coisolationIndex, indexResults.PeptideIndex, indexResults.FragmentIndex,
                indexResults.PrecursorIndex, 0, commonParameters, null, variableMods, TenPpm, 0, new List<string>()).Run();
            return matches;
        }

        /// <summary>The best score Modern search gives the intact sixmer, which a trimmed GUACUG must equal.</summary>
        private static double ModernSearchScoreOfTheSixmer()
        {
            var commonParameters = RnaCommonParameters(new RnaDigestionParams("top-down"));
            var scans = SixmerScans(commonParameters);
            var indexResults = (IndexingResults)new IndexingEngine([new RNA("GUACUG")], [], [], null, null, null, 0, DecoyType.None,
                commonParameters, null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();
            var matches = new SpectralMatch[scans.Length];
            new ModernSearchEngine(matches, scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0,
                commonParameters, null, TenPpm, 0, ["reference"]).Run();
            return matches.Where(m => m != null).Max(m => m.Score);
        }

        private static SpectralMatch BestNonSpecificMatch(SpectralMatch[][] matches)
        {
            Assert.That(matches[(int)FdrCategory.FullySpecific].All(m => m == null), "a nucleic acid candidate has no specific rnase to be fully specific to");
            Assert.That(matches[(int)FdrCategory.SemiSpecific].All(m => m == null));
            // the engine leaves ambiguity to ResolveFdrCategorySpecificPsms, so sequences are unset until resolved
            var found = matches[(int)FdrCategory.NonSpecific].Where(m => m != null).ToList();
            found.ForEach(m => m.ResolveAllAmbiguities());
            return found.OrderByDescending(m => m.Score).FirstOrDefault();
        }

        #region Trimming arithmetic

        /// <summary>
        /// Every sub-oligo mass the search tests equals the mass mzLib gives the same oligo built independently: the
        /// kept nucleotides, the fixed end's terminus, and the rnase's terminus at the cut.
        /// </summary>
        [Test]
        [TestCase("singleN")]
        [TestCase("singleC")]
        public static void SubOligoMasses_EqualTheMassesOfTheSameOligosBuiltDirectly(string rnase)
        {
            const string sequence = "GUACUGCCAU";
            var digestionParams = Seeds(rnase, FragmentationTerminus.Both);
            var seeds = new RNA(sequence).Digest(digestionParams, [], []).Cast<OligoWithSetMods>().ToList();
            Assert.That(seeds, Has.Count.EqualTo(sequence.Length), "premise: one seed per position");

            FragmentationTerminus fixedTerminus = OligoSeedTrimming.FixedTerminus(digestionParams.Rnase).Value;
            IHasChemicalFormula openEnd = OligoSeedTrimming.OpenEndTermini(digestionParams.Rnase).Single();
            int checkedSubOligos = 0;
            foreach (var seed in seeds)
            {
                double[] masses = OligoSeedTrimming.MassesWithoutOpenEndTerminus(seed, fixedTerminus);
                for (int n = 1; n < seed.BaseSequence.Length; n++)
                {
                    bool fivePrimeFixed = fixedTerminus == FragmentationTerminus.FivePrime;
                    string kept = fivePrimeFixed ? seed.BaseSequence[..n] : seed.BaseSequence[^n..];
                    var direct = new OligoWithSetMods(kept,
                        fivePrimeTerminus: fivePrimeFixed ? seed.FivePrimeTerminus : openEnd,
                        threePrimeTerminus: fivePrimeFixed ? openEnd : seed.ThreePrimeTerminus);

                    Assert.That(masses[n] + openEnd.MonoisotopicMass, Is.EqualTo(direct.MonoisotopicMass).Within(1e-9), $"{seed.BaseSequence} kept {n}");

                    var trimmed = OligoSeedTrimming.Trim(seed, fixedTerminus, n, openEnd);
                    Assert.That(trimmed.BaseSequence, Is.EqualTo(kept));
                    Assert.That(trimmed.MonoisotopicMass, Is.EqualTo(direct.MonoisotopicMass).Within(1e-9));
                    Assert.That(trimmed.CleavageSpecificityForFdrCategory, Is.EqualTo(CleavageSpecificity.None));
                    checkedSubOligos++;
                }
            }
            Assert.That(checkedSubOligos, Is.EqualTo(sequence.Length * (sequence.Length - 1) / 2));
        }

        /// <summary>
        /// Modifications move with their nucleotides: a trimmed oligo keeps the ones on the nucleotides it keeps (at
        /// shifted keys for a 3'-fixed seed), drops the rest, and gains the open-end modification at the new terminus.
        /// </summary>
        [Test]
        [TestCase("singleN")]
        [TestCase("singleC")]
        public static void Trim_KeepsTheModificationsOnKeptNucleotides_AndAddsTheOpenEndModification(string rnase)
        {
            ModificationMotif.TryGetMotif("U", out var uridine);
            var methyl = new Modification("Methyl", _modificationType: "Test", _target: uridine, _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("CH2"));
            var terminal = new Modification("TerminalTest", _modificationType: "Test", _target: uridine, _locationRestriction: "3'-terminal.",
                _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O"));

            var digestionParams = Seeds(rnase, FragmentationTerminus.Both);
            FragmentationTerminus fixedTerminus = OligoSeedTrimming.FixedTerminus(digestionParams.Rnase).Value;
            bool fivePrimeFixed = fixedTerminus == FragmentationTerminus.FivePrime;
            var rna = new RNA("GUACUG");
            var seed = rna.Digest(digestionParams, [], []).Cast<OligoWithSetMods>().First(s => s.BaseSequence == "GUACUG");

            // methyl on U2 (key 3) and on U5 (key 6)
            var modified = new OligoWithSetMods(rna, digestionParams, 1, 6, 0, seed.CleavageSpecificityForFdrCategory,
                new Dictionary<int, Modification> { [3] = methyl, [6] = methyl }, 0, seed.FivePrimeTerminus, seed.ThreePrimeTerminus);
            IHasChemicalFormula openEnd = OligoSeedTrimming.OpenEndTermini(digestionParams.Rnase).Single();

            var trimmed = OligoSeedTrimming.Trim(modified, fixedTerminus, 3, openEnd, terminal);

            // 5' fixed keeps GUA (methyl on U2 stays at key 3); 3' fixed keeps CUG (methyl on U5 moves to key 3)
            Assert.That(trimmed.BaseSequence, Is.EqualTo(fivePrimeFixed ? "GUA" : "CUG"));
            Assert.That(trimmed.AllModsOneIsNterminus[3], Is.SameAs(methyl));
            Assert.That(trimmed.AllModsOneIsNterminus[fivePrimeFixed ? 5 : 1], Is.SameAs(terminal), "the open-end modification sits on the new terminus");
            Assert.That(trimmed.AllModsOneIsNterminus, Has.Count.EqualTo(2));

            double[] masses = OligoSeedTrimming.MassesWithoutOpenEndTerminus(modified, fixedTerminus);
            Assert.That(masses[3] + openEnd.MonoisotopicMass + terminal.MonoisotopicMass.Value,
                Is.EqualTo(trimmed.MonoisotopicMass).Within(1e-9), "the mass tested before trimming is the mass of what is reported");
        }

        /// <summary>
        /// The open end's modifications are found by the nucleic acid location restrictions. Matching "C-terminal" or
        /// "N-terminal", as the peptide search does, finds none of them.
        /// </summary>
        [Test]
        public static void VariableTerminalMods_ForNucleicAcids_AreTheOnesOnTheOpenEnd()
        {
            ModificationMotif.TryGetMotif("X", out var any);
            Modification Mod(string restriction) => new(restriction, _modificationType: "Test", _target: any,
                _locationRestriction: restriction, _chemicalFormula: ChemicalFormula.ParseFormula("H"));
            var mods = new List<Modification> { Mod("5'-terminal."), Mod("Oligo 5'-terminal."), Mod("3'-terminal."), Mod("Oligo 3'-terminal."), Mod("Anywhere.") };

            var singleN = NonSpecificEnzymeSearchEngine.GetVariableTerminalMods(
                NonSpecificEnzymeSearchEngine.SeedTerminus(Seeds("singleN", FragmentationTerminus.Both)), mods);
            var singleC = NonSpecificEnzymeSearchEngine.GetVariableTerminalMods(
                NonSpecificEnzymeSearchEngine.SeedTerminus(Seeds("singleC", FragmentationTerminus.Both)), mods);

            Assert.That(singleN.Select(m => m.IdWithMotif), Is.EquivalentTo(new[] { "3'-terminal. on X", "Oligo 3'-terminal. on X" }));
            Assert.That(singleC.Select(m => m.IdWithMotif), Is.EquivalentTo(new[] { "5'-terminal. on X", "Oligo 5'-terminal. on X" }));
        }

        /// <summary>For nucleic acids the rnase decides which end is fixed, whatever FragmentationTerminus says.</summary>
        [Test]
        [TestCase("singleN", FragmentationTerminus.Both, FragmentationTerminus.FivePrime)]
        [TestCase("singleN", FragmentationTerminus.N, FragmentationTerminus.FivePrime)]
        [TestCase("singleC", FragmentationTerminus.Both, FragmentationTerminus.ThreePrime)]
        [TestCase("RNase T1", FragmentationTerminus.Both, FragmentationTerminus.Both)]
        public static void SeedTerminus_ForNucleicAcids_ComesFromTheRnase(string rnase, FragmentationTerminus terminus, FragmentationTerminus expected)
        {
            Assert.That(NonSpecificEnzymeSearchEngine.SeedTerminus(new RnaDigestionParams(rnase, fragmentationTerminus: terminus)), Is.EqualTo(expected));
        }

        #endregion

        #region Engine

        /// <summary>
        /// singleC over AAGUACUG gives the seed AAGUACUG and trims its 5' end to GUACUG, which carries 5'-OH, exactly
        /// the intact sixmer. It is reported as an oligo match in the non-specific category, with the score Modern
        /// search gives the sixmer itself.
        /// </summary>
        [Test]
        public static void SingleC_TrimsTheFivePrimeEnd_AndFindsTheSixmerWithTheModernSearchScore()
        {
            var commonParameters = RnaCommonParameters(Seeds("singleC", FragmentationTerminus.ThreePrime));

            var best = BestNonSpecificMatch(RunNonSpecific(commonParameters, [new RNA("AAGUACUG")]));

            Assert.That(best, Is.Not.Null, "the trimmed sixmer was not found");
            Assert.That(best, Is.TypeOf<OligoSpectralMatch>());
            Assert.That(best.BaseSequence, Is.EqualTo("GUACUG"));
            Assert.That(best.OneBasedStartResidue, Is.EqualTo(3));
            Assert.That(best.Score, Is.EqualTo(ModernSearchScoreOfTheSixmer()).Within(1e-9));
        }

        /// <summary>
        /// singleN over GUACUGAA trims the 3' end, and the cut leaves the rnase's 3'-phosphate, which the sixmer does
        /// not have. So nothing is found: the engine uses the terminus, not a peptide-style water. A custom singleN
        /// whose cut leaves 3'-OH finds it -- proving the terminus comes from the rnase.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void SingleN_TrimsTheThreePrimeEnd_WithTheTerminusTheRnaseLeaves()
        {
            List<IBioPolymer> targets = [new RNA("GUACUGAA")];

            var withPhosphate = BestNonSpecificMatch(RunNonSpecific(RnaCommonParameters(Seeds("singleN", FragmentationTerminus.FivePrime)), targets));
            Assert.That(withPhosphate?.BaseSequence, Is.Not.EqualTo("GUACUG"), "a 3'-phosphate GUACUG must not match the 3'-OH sixmer");

            const string customName = "singleN leaving 3'-OH (test)";
            RnaseDictionary.Dictionary[customName] = new Rnase(customName, CleavageSpecificity.SingleN, [],
                threePrimeTerminusRemainder: [ThreePrimeHydroxyl]);
            try
            {
                var withHydroxyl = BestNonSpecificMatch(RunNonSpecific(RnaCommonParameters(Seeds(customName, FragmentationTerminus.FivePrime)), targets));

                Assert.That(withHydroxyl, Is.Not.Null, "the trimmed sixmer was not found");
                Assert.That(withHydroxyl.BaseSequence, Is.EqualTo("GUACUG"));
                Assert.That(withHydroxyl.OneBasedStartResidue, Is.EqualTo(1));
                Assert.That(withHydroxyl.Score, Is.EqualTo(ModernSearchScoreOfTheSixmer()).Within(1e-9));
            }
            finally
            {
                RnaseDictionary.Dictionary.Remove(customName);
            }
        }

        /// <summary>A seed that is the whole molecule is found through the precursor index, untrimmed.</summary>
        [Test]
        [TestCase("singleN", FragmentationTerminus.FivePrime)]
        [TestCase("singleC", FragmentationTerminus.ThreePrime)]
        public static void AWholeMoleculeSeed_IsFoundUntrimmed(string rnase, FragmentationTerminus terminus)
        {
            var best = BestNonSpecificMatch(RunNonSpecific(RnaCommonParameters(Seeds(rnase, terminus)), [new RNA("GUACUG")]));

            Assert.That(best, Is.Not.Null);
            Assert.That(best.BaseSequence, Is.EqualTo("GUACUG"));
            Assert.That(best.Score, Is.EqualTo(ModernSearchScoreOfTheSixmer()).Within(1e-9));
        }

        /// <summary>
        /// A 3'-terminal variable modification is tried on the trimmed 3' end of a singleN seed, and reported on the
        /// match. singleN over GUACUGAA gives GUACUG with a 3'-phosphate, HPO3 heavier than the sixmer; a made-up
        /// 3'-terminal modification of minus HPO3 on G turns it back into the sixmer, so only a search that applies the
        /// modification to the trimmed end finds it.
        /// </summary>
        [Test]
        public static void SingleN_TriesOpenEndModifications_OnTheTrimmedEnd()
        {
            ModificationMotif.TryGetMotif("G", out var guanosine);
            var removesPhosphate = new Modification("RemovesPhosphate", _modificationType: "Test", _target: guanosine,
                _locationRestriction: "3'-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H-1O-3P-1"));
            var commonParameters = RnaCommonParameters(Seeds("singleN", FragmentationTerminus.FivePrime));

            var best = BestNonSpecificMatch(RunNonSpecific(commonParameters, [new RNA("GUACUGAA")], [removesPhosphate]));

            Assert.That(best, Is.Not.Null, "the modified, trimmed sixmer was not found");
            Assert.That(best.BaseSequence, Is.EqualTo("GUACUG"));
            Assert.That(best.FullSequence, Does.Contain("RemovesPhosphate"), "the open-end modification must be on the reported oligo");
            Assert.That(best.Score, Is.EqualTo(ModernSearchScoreOfTheSixmer()).Within(1e-9));
        }

        /// <summary>
        /// A seed whose interior nucleotide can carry an open-end modification is also filed in the precursor index at the
        /// mass of the modified sub-oligo ending there, so the search considers it for that precursor. The peptide
        /// version of this threw for oligos.
        /// </summary>
        /// <remarks>
        /// Tested on the index directly: in <see cref="SingleN_TriesOpenEndModifications_OnTheTrimmedEnd"/> a c ion
        /// less HPO3 weighs exactly what an a ion weighs, so the fragment route finds that seed without these bins.
        /// </remarks>
        [Test]
        public static void PrecursorIndex_FilesSeedsAtTheMassOfTheirModifiedSubOligos()
        {
            ModificationMotif.TryGetMotif("C", out var cytidine);
            var terminal = new Modification("TerminalOnC", _modificationType: "Test", _target: cytidine,
                _locationRestriction: "3'-terminal.", _monoisotopicMass: 100);
            var commonParameters = RnaCommonParameters(Seeds("singleN", FragmentationTerminus.FivePrime));

            var indexResults = (IndexingResults)new IndexingEngine([new RNA("GUACUGAA")], [terminal], [], null, null, null, 0, DecoyType.None,
                commonParameters, null, 30000, true, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            int seedId = indexResults.PeptideIndex.FindIndex(p => p.BaseSequence == "GUACUGAA");
            double subOligoMass = new OligoWithSetMods("GUAC", threePrimeTerminus: ThreePrimePhosphate).MonoisotopicMass + 100;
            int bin = (int)Math.Round(subOligoMass * 1000); // FragmentBinsPerDalton

            Assert.That(seedId, Is.GreaterThanOrEqualTo(0), "premise: the seed is indexed");
            Assert.That(indexResults.PrecursorIndex[bin], Does.Contain(seedId), "GUAC + terminal mod is not in the precursor index");
        }

        [Test]
        public static void Engine_RefusesAnRnaseThatGivesNoSeeds()
        {
            var commonParameters = RnaCommonParameters(new RnaDigestionParams("RNase T1"));
            var thrown = Assert.Throws<MetaMorpheusException>(() =>
                new NonSpecificEnzymeSearchEngine([], [], [], [], null, null, 0, commonParameters, null, [], TenPpm, 0, []));
            Assert.That(thrown.Message, Does.Contain("singleN or singleC"));
        }

        [Test]
        public static void Engine_RefusesComplementaryIonsForNucleicAcids()
        {
            var commonParameters = RnaCommonParameters(Seeds("singleN", FragmentationTerminus.FivePrime), addCompIons: true);
            var thrown = Assert.Throws<MetaMorpheusException>(() =>
                new NonSpecificEnzymeSearchEngine([], [], [], [], null, null, 0, commonParameters, null, [], TenPpm, 0, []));
            Assert.That(thrown.Message, Does.Contain("Complementary ions"));
        }

        #endregion

        #region Task

        /// <summary>
        /// End to end through SearchTask: a non-specific search over a nucleic acid database, which used to be refused,
        /// runs, sets the terminus from the rnase (the settings say Both), and writes the trimmed sixmer as an OSM.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void SearchTask_NonSpecificOverANucleicAcidDatabase_ReportsTheTrimmedOligo()
        {
            var original = GlobalVariables.AnalyteType;
            string outputDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "RnaNonSpecificSearchTask");
            if (Directory.Exists(outputDir)) Directory.Delete(outputDir, true);
            Directory.CreateDirectory(outputDir);

            try
            {
                GlobalVariables.AnalyteType = AnalyteType.Oligo;
                string database = Path.Combine(outputDir, "extended6mer.fasta");
                File.WriteAllLines(database, [
                    ">id:1|Name:extended6mer|SOterm:6mer|Type:tRNA|Subtype:Ala|Feature:VGC|Cellular_Localization:freezer|Species:standard",
                    "AAGUACUG"]);

                var task = new SearchTask
                {
                    CommonParameters = RnaCommonParameters(new RnaDigestionParams("singleC", minLength: 3)),
                    SearchParameters = new RnaSearchParameters
                    {
                        SearchType = SearchType.NonSpecific,
                        DecoyType = DecoyType.Reverse,
                        MassDiffAcceptorType = MassDiffAcceptorType.Custom,
                        CustomMdac = "Custom interval [-5,5]",
                        DisposeOfFileWhenDone = true
                    }
                };
                Assert.That(task.GetSeedDigestionRefusal(), Is.Null, "premise: the run is not refused");

                task.RunTask(outputDir, [new DbForTask(database, false)], [SixmerFilePath], "NonSpecificRna");

                string osms = Path.Combine(outputDir, "AllOSMs.osmtsv");
                Assert.That(File.Exists(osms), "no OSM file was written");
                var lines = File.ReadAllLines(osms);
                var header = lines[0].Split('\t');
                int baseSequence = Array.IndexOf(header, SpectrumMatchFromTsvHeader.BaseSequence);
                int start = Array.IndexOf(header, SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence);
                var rows = lines.Skip(1).Select(l => l.Split('\t')).ToList();
                Assert.That(rows.Any(r => r[baseSequence] == "GUACUG" && r[start].Contains("[3 to 8]")),
                    "the trimmed sixmer is not in the OSM file");
            }
            finally
            {
                GlobalVariables.AnalyteType = original;
                if (Directory.Exists(outputDir)) Directory.Delete(outputDir, true);
            }
        }

        #endregion
    }
}
