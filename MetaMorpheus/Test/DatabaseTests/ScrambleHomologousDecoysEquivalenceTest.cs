using EngineLayer;
using EngineLayer.DatabaseLoading;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    /// <summary>
    /// ScrambleHomologousDecoys looks for decoy digestion products whose base sequence is also a target's, and scrambles
    /// those decoys. It digests without enumerating modified forms, since only base sequences are compared. These tests
    /// hold it to the previous implementation, kept verbatim below: the same decoys replaced, by equal objects - same
    /// sequence, accession and modification dictionaries in the same key order holding the same modification objects -
    /// and every other entry left as the very same object, on one thread or several.
    /// </summary>
    [TestFixture]
    internal static class ScrambleHomologousDecoysEquivalenceTest
    {
        // K and R twice, so random proteins digest into tryptic peptides of ordinary length
        private const string AminoAcids = "ACDEFGHIKLMNPQRSTVWYKR";
        private const string Nucleotides = "ACGU";

        [Test]
        [TestCase("trypsin", 2, 7, CleavageSpecificity.Full, InitiatorMethionineBehavior.Variable, 1)]
        [TestCase("trypsin", 0, 5, CleavageSpecificity.Full, InitiatorMethionineBehavior.Retain, 2)]
        [TestCase("trypsin", 1, 6, CleavageSpecificity.Semi, InitiatorMethionineBehavior.Cleave, 3)]
        [TestCase("Arg-C", 2, 5, CleavageSpecificity.Full, InitiatorMethionineBehavior.Variable, 4)]
        [TestCase("top-down", 0, 7, CleavageSpecificity.Full, InitiatorMethionineBehavior.Variable, 5)]
        public static void ProteinDecoysAreScrambledAsBefore(string protease, int missedCleavages, int minLength,
            CleavageSpecificity searchMode, InitiatorMethionineBehavior initiatorMethionine, int seed)
        {
            var random = new Random(seed);
            var digestionParams = new DigestionParams(protease, missedCleavages, minLength, maxModificationIsoforms: 1024,
                initiatorMethionineBehavior: initiatorMethionine, maxModsForPeptides: 2, searchModeType: searchMode);
            List<IBioPolymer> database = RandomProteinDatabase(random, digestionParams, 300);

            int scrambled = AssertScrambledAsBefore(database, digestionParams);
            Assert.That(scrambled, Is.GreaterThan(20), "the planted collisions must exercise the scrambling path");
        }

        [Test]
        [TestCase("RNase T1", 0, 3, 11)]
        [TestCase("RNase T1", 2, 4, 12)]
        [TestCase("RNase A", 1, 3, 13)]
        [TestCase("top-down", 0, 3, 14)]
        public static void RnaDecoysAreScrambledAsBefore(string rnase, int missedCleavages, int minLength, int seed)
        {
            var random = new Random(seed);
            var digestionParams = new RnaDigestionParams(rnase, missedCleavages, minLength, maxModificationIsoforms: 1024, maxMods: 2);
            List<IBioPolymer> database = RandomRnaDatabase(random, digestionParams, 200);

            int scrambled = AssertScrambledAsBefore(database, digestionParams);
            Assert.That(scrambled, Is.GreaterThan(10), "the planted collisions must exercise the scrambling path");
        }

        /// <summary>
        /// A real UniProt XML, with its localized modifications and the reverse decoys the loader makes. Few of those
        /// collide on their own, so a second copy of the database adds decoys that do.
        /// </summary>
        [Test]
        [TestCase(CleavageSpecificity.Full, 2, 7)]
        [TestCase(CleavageSpecificity.Semi, 1, 5)]
        public static void LoadedDatabaseIsScrambledAsBefore(CleavageSpecificity searchMode, int missedCleavages, int minLength)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "LowResSnip_B6_mouse_11700_117500.xml");
            var digestionParams = new DigestionParams("trypsin", missedCleavages, minLength, searchModeType: searchMode);
            var commonParameters = new CommonParameters(digestionParams: digestionParams);
            List<IBioPolymer> database = DatabaseLoadingEngine.LoadBioPolymers("test", new List<DbForTask> { new(path, false) }, true,
                DecoyType.Reverse, GlobalVariables.AllModTypesKnown.ToList(), commonParameters, out _);
            Assert.That(database.Any(p => p.OneBasedPossibleLocalizedModifications.Count > 0), "the database must carry localized modifications");

            // decoys that splice target peptides into reversed targets, keeping the real modifications
            var random = new Random(7);
            var targets = database.Where(p => !p.IsDecoy).Cast<Protein>().ToList();
            foreach (Protein target in targets)
            {
                Protein donor = targets[random.Next(targets.Count)];
                string decoySequence = SpliceTargetProducts(random, new string(target.BaseSequence.Reverse().ToArray()), donor, digestionParams);
                database.Add(new Protein(decoySequence, "PLANTED_" + target.Accession, isDecoy: true,
                    oneBasedModifications: ShiftedModifications(target.OneBasedPossibleLocalizedModifications, decoySequence.Length)));
            }

            int scrambled = AssertScrambledAsBefore(database, digestionParams);
            Assert.That(scrambled, Is.GreaterThan(10));
        }

        [Test]
        public static void CollisionDigestionDropsModifiedFormsOnly()
        {
            var digestionParams = new DigestionParams("trypsin", 2, 7, 64, 1024, InitiatorMethionineBehavior.Variable, 3,
                CleavageSpecificity.Semi, Omics.Fragmentation.FragmentationTerminus.Both, false, true, true);
            var collision = (DigestionParams)DatabaseLoadingEngine.GetCollisionDigestionParams(digestionParams);
            Assert.That(collision, Is.Not.SameAs(digestionParams), "the caller's parameters must not be changed");
            Assert.That(digestionParams.MaxMods, Is.EqualTo(3));
            Assert.That(collision.MaxMods, Is.EqualTo(0));
            Assert.That(collision.Protease, Is.SameAs(digestionParams.Protease));
            Assert.That(collision.SpecificProtease, Is.SameAs(digestionParams.SpecificProtease));
            Assert.That((collision.MaxMissedCleavages, collision.MinLength, collision.MaxLength, collision.MaxModificationIsoforms,
                    collision.InitiatorMethionineBehavior, collision.SearchModeType, collision.FragmentationTerminus,
                    collision.GeneratehUnlabeledProteinsForSilac, collision.KeepNGlycopeptide, collision.KeepOGlycopeptide),
                Is.EqualTo((digestionParams.MaxMissedCleavages, digestionParams.MinLength, digestionParams.MaxLength, digestionParams.MaxModificationIsoforms,
                    digestionParams.InitiatorMethionineBehavior, digestionParams.SearchModeType, digestionParams.FragmentationTerminus,
                    digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParams.KeepNGlycopeptide, digestionParams.KeepOGlycopeptide)));

            var rnaParams = new RnaDigestionParams("RNase T1", 1, 4, 50, 128, 2);
            var rnaCollision = (RnaDigestionParams)DatabaseLoadingEngine.GetCollisionDigestionParams(rnaParams);
            Assert.That(rnaCollision, Is.Not.SameAs(rnaParams));
            Assert.That((rnaCollision.MaxMods, rnaParams.MaxMods), Is.EqualTo((0, 2)));
            Assert.That(rnaCollision.Rnase, Is.SameAs(rnaParams.Rnase));
        }

        /// <summary>
        /// With mzLib's opt-in RespectCleavageBlockingModifications on, a product can exist only in a modified form, so the
        /// collision check must keep digesting with the parameters it was given.
        /// </summary>
        [Test]
        public static void CollisionDigestionKeepsModifiedFormsWhenModificationsBlockCleavage()
        {
            var blocking = new DigestionParams("trypsin", 0, 5, 64, 1024, InitiatorMethionineBehavior.Retain, 1,
                respectCleavageBlockingModifications: true);
            Assert.That(DatabaseLoadingEngine.GetCollisionDigestionParams(blocking), Is.SameAs(blocking));

            var notBlocking = new DigestionParams("trypsin", 0, 5, 64, 1024, InitiatorMethionineBehavior.Retain, 1);
            var collision = DatabaseLoadingEngine.GetCollisionDigestionParams(notBlocking);
            Assert.That(collision.MaxMods, Is.EqualTo(0));

            // why: with acetylated K blocking trypsin, the read-through product across the acetylated K is yielded only
            // in its modified form, so digesting with MaxMods = 0 would not see its base sequence at all
            var protein = new Protein("PEPTIDEKAAAAAAKGGGGGR", "blocking");
            var acetylK = GlobalVariables.AllModsKnown.Where(m => m.IdWithMotif == "Acetylation on K").ToList();
            Assert.That(acetylK.Any(m => m.BlocksCleavage), Is.True);
            HashSet<string> BaseSequences(IDigestionParams p) =>
                protein.Digest(p, new List<Modification>(), acetylK).Select(pep => pep.BaseSequence).ToHashSet();
            HashSet<string> searched = BaseSequences(blocking);
            HashSet<string> withoutModifiedForms = BaseSequences(new DigestionParams("trypsin", 0, 5, 64, 1024,
                InitiatorMethionineBehavior.Retain, 0, respectCleavageBlockingModifications: true));
            Assert.That(searched, Does.Contain("PEPTIDEKAAAAAAK"));
            Assert.That(withoutModifiedForms, Does.Not.Contain("PEPTIDEKAAAAAAK"));

            // and with the flag off, MaxMods = 0 yields exactly the base sequences of the full digestion
            Assert.That(BaseSequences(collision), Is.EquivalentTo(BaseSequences(notBlocking)));
        }

        /// <summary>
        /// Likewise with mzLib's opt-in RespectCleavagePromotingModifications on: StcE cuts RPPITQ|SSL only when Thr5, at P2,
        /// carries an O-glycan, so RPPITQ is yielded only glycosylated.
        /// </summary>
        [Test]
        public static void CollisionDigestionKeepsModifiedFormsWhenModificationsPromoteCleavage()
        {
            var promoting = new DigestionParams("StcE", 0, 1, 64, 1024, InitiatorMethionineBehavior.Retain, 1,
                respectCleavagePromotingModifications: true);
            Assert.That(DatabaseLoadingEngine.GetCollisionDigestionParams(promoting), Is.SameAs(promoting));

            ModificationMotif.TryGetMotif("T", out ModificationMotif threonine);
            var oGlycan = new Modification(_originalId: "HexNAc", _modificationType: "O-linked glycosylation", _target: threonine,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373);
            var protein = new Protein("RPPITQSSL", "promoting",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 5, new List<Modification> { oGlycan } } });
            HashSet<string> BaseSequences(IDigestionParams p) =>
                protein.Digest(p, new List<Modification>(), new List<Modification>()).Select(pep => pep.BaseSequence).ToHashSet();

            // why: digesting with MaxMods = 0 would not see RPPITQ's base sequence at all
            HashSet<string> withoutModifiedForms = BaseSequences(new DigestionParams("StcE", 0, 1, 64, 1024,
                InitiatorMethionineBehavior.Retain, 0, respectCleavagePromotingModifications: true));
            Assert.That(BaseSequences(promoting), Does.Contain("RPPITQ"));
            Assert.That(withoutModifiedForms, Does.Not.Contain("RPPITQ"));
        }

        #region comparison

        /// <summary>
        /// Runs the previous implementation and the current one on separate copies of the same list and asserts the
        /// outcomes are the same. Returns how many decoys were scrambled.
        /// </summary>
        private static int AssertScrambledAsBefore(List<IBioPolymer> database, IDigestionParams digestionParams)
        {
            var expected = new List<IBioPolymer>(database);
            int expectedScrambled = ScrambleHomologousDecoysBefore(expected, digestionParams);

            // one thread, and several: the outcome must not depend on how the work is divided
            foreach (int threads in new[] { 1, 3, 8 })
            {
                var actual = new List<IBioPolymer>(database);
                int actualScrambled = DatabaseLoadingEngine.ScrambleHomologousDecoys(actual, digestionParams, threads);

                Assert.That(actualScrambled, Is.EqualTo(expectedScrambled), $"{threads} threads");
                AssertSameDatabase(database, expected, actual);
            }
            return expectedScrambled;
        }

        private static void AssertSameDatabase(List<IBioPolymer> input, List<IBioPolymer> expected, List<IBioPolymer> actual)
        {
            Assert.That(actual.Count, Is.EqualTo(expected.Count));
            for (int i = 0; i < expected.Count; i++)
            {
                if (ReferenceEquals(expected[i], input[i]))
                {
                    Assert.That(actual[i], Is.SameAs(input[i]), $"entry {i} ({input[i].Accession}) must be left as it was");
                    continue;
                }

                Assert.That(actual[i], Is.Not.SameAs(input[i]), $"entry {i} ({input[i].Accession}) must be scrambled");
                Assert.That(actual[i].GetType(), Is.EqualTo(expected[i].GetType()));
                Assert.That(actual[i].BaseSequence, Is.EqualTo(expected[i].BaseSequence), $"sequence of entry {i}");
                Assert.That((actual[i].Accession, actual[i].Name, actual[i].IsDecoy, actual[i].IsContaminant),
                    Is.EqualTo((expected[i].Accession, expected[i].Name, expected[i].IsDecoy, expected[i].IsContaminant)));
                AssertSameModifications(expected[i].OneBasedPossibleLocalizedModifications, actual[i].OneBasedPossibleLocalizedModifications, i);
                AssertSameModifications(expected[i].OriginalNonVariantModifications, actual[i].OriginalNonVariantModifications, i);
            }
        }

        // same keys in the same order, each holding the same modification objects in the same order
        private static void AssertSameModifications(IDictionary<int, List<Modification>> expected, IDictionary<int, List<Modification>> actual, int entry)
        {
            Assert.That(actual.Keys.ToList(), Is.EqualTo(expected.Keys.ToList()), $"modification positions of entry {entry}");
            foreach (int key in expected.Keys)
            {
                Assert.That(actual[key].Count, Is.EqualTo(expected[key].Count));
                for (int m = 0; m < expected[key].Count; m++)
                {
                    Assert.That(actual[key][m], Is.SameAs(expected[key][m]), $"modification {m} at {key} of entry {entry}");
                }
            }
        }

        #endregion

        #region the previous implementation, verbatim

        private static readonly ListPool<string> SequencePool = new();

        public static int ScrambleHomologousDecoysBefore<TBioPolymer>(List<TBioPolymer> bioPolymers, IDigestionParams digestionParams) where TBioPolymer : IBioPolymer
        {
            // Sanitize the decoys
            // TODO: Fix this so that it accounts for multi-protease searches. Currently, we only consider the first protease
            // when looking for target/decoy collisions
            HashSet<string> targetPeptideSequences = new();
            foreach (var bioPolymer in bioPolymers.Where(p => !p.IsDecoy))
            {
                // When thinking about decoy collisions, we can ignore modifications
                foreach (var peptide in bioPolymer.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    targetPeptideSequences.Add(peptide.BaseSequence);
                }
            }

            var variable = new List<Modification>();
            var fixedM = new List<Modification>();
            int scrambled = 0;

            // Now, we iterate through the decoys and scramble the sequences that correspond to target peptides
            for (int i = 0; i < bioPolymers.Count; i++)
            {
                var bioPolymer = bioPolymers[i];
                if (!bioPolymer.IsDecoy)
                    continue;

                // Use a pooled list to avoid repeated allocations in tight loops
                List<string> peptidesToReplace = null;
                foreach (var peptide in bioPolymer.Digest(digestionParams, variable, fixedM))
                {
                    if (targetPeptideSequences.Contains(peptide.BaseSequence))
                    {
                        peptidesToReplace ??= SequencePool.Get();
                        peptidesToReplace.Add(peptide.BaseSequence);
                    }
                }

                if (peptidesToReplace is { Count: > 0 })
                {
                    bioPolymers[i] = DecoySequenceValidator.ScrambleDecoyBioPolymer(
                        bioPolymer,
                        digestionParams,
                        forbiddenSequences: targetPeptideSequences,
                        peptidesToReplace);
                    SequencePool.Return(peptidesToReplace);
                    scrambled++;
                }
            }
            return scrambled;
        }

        #endregion

        #region random databases

        /// <summary>
        /// Targets with random localized modifications, and for each one a decoy: its reverse, with digestion products of
        /// random targets spliced in (so it collides), also modified. One decoy in ten is an exact copy of a target, which
        /// is the only way a top-down decoy collides.
        /// </summary>
        private static List<IBioPolymer> RandomProteinDatabase(Random random, DigestionParams digestionParams, int targetCount)
        {
            List<Modification> mods = ProteinModifications();
            var targets = new List<Protein>();
            for (int t = 0; t < targetCount; t++)
            {
                int length = 20 + random.Next(280);
                char[] sequence = new char[length];
                for (int r = 0; r < length; r++)
                {
                    sequence[r] = AminoAcids[random.Next(AminoAcids.Length)];
                }
                if (random.Next(2) == 0)
                {
                    sequence[0] = 'M';
                }
                string s = new(sequence);
                targets.Add(new Protein(s, "P" + t, name: "target " + t, oneBasedModifications: RandomModifications(random, s, mods)));
            }

            var database = new List<IBioPolymer>(targets);
            for (int t = 0; t < targetCount; t++)
            {
                Protein target = targets[t];
                string decoySequence = random.Next(10) == 0
                    ? target.BaseSequence
                    : SpliceTargetProducts(random, new string(target.BaseSequence.Reverse().ToArray()), targets[random.Next(targetCount)], digestionParams);
                database.Add(new Protein(decoySequence, "DECOY_P" + t, name: "decoy " + t, isDecoy: true,
                    oneBasedModifications: RandomModifications(random, decoySequence, mods)));
            }
            return database;
        }

        private static List<IBioPolymer> RandomRnaDatabase(Random random, RnaDigestionParams digestionParams, int targetCount)
        {
            List<Modification> mods = RnaModifications();
            var targets = new List<RNA>();
            for (int t = 0; t < targetCount; t++)
            {
                int length = 12 + random.Next(100);
                char[] sequence = new char[length];
                for (int r = 0; r < length; r++)
                {
                    sequence[r] = Nucleotides[random.Next(Nucleotides.Length)];
                }
                string s = new(sequence);
                targets.Add(new RNA(s, "R" + t, RandomModifications(random, s, mods), name: "target " + t));
            }

            var database = new List<IBioPolymer>(targets);
            for (int t = 0; t < targetCount; t++)
            {
                RNA target = targets[t];
                string decoySequence = random.Next(10) == 0
                    ? target.BaseSequence
                    : SpliceTargetProducts(random, new string(target.BaseSequence.Reverse().ToArray()), targets[random.Next(targetCount)], digestionParams);
                database.Add(new RNA(decoySequence, "DECOY_R" + t, RandomModifications(random, decoySequence, mods), name: "decoy " + t, isDecoy: true));
            }
            return database;
        }

        // inserts one to three digestion products of the donor into the sequence, at random points
        private static string SpliceTargetProducts(Random random, string sequence, IBioPolymer donor, IDigestionParams digestionParams)
        {
            List<string> products = donor.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Select(p => p.BaseSequence).Distinct().ToList();
            if (products.Count == 0)
            {
                return sequence;
            }
            int splices = 1 + random.Next(3);
            for (int s = 0; s < splices; s++)
            {
                int at = random.Next(sequence.Length + 1);
                sequence = sequence.Insert(at, products[random.Next(products.Count)]);
            }
            return sequence;
        }

        // keeps the modifications whose position still lies inside the new sequence
        private static Dictionary<int, List<Modification>> ShiftedModifications(IDictionary<int, List<Modification>> mods, int length)
            => mods.Where(kv => kv.Key >= 1 && kv.Key <= length).ToDictionary(kv => kv.Key, kv => kv.Value.ToList());

        // about one residue in eight carries one or two modifications whose motif fits it
        private static Dictionary<int, List<Modification>> RandomModifications(Random random, string sequence, List<Modification> mods)
        {
            var result = new Dictionary<int, List<Modification>>();
            for (int position = 1; position <= sequence.Length; position++)
            {
                if (random.Next(8) != 0)
                {
                    continue;
                }
                char residue = sequence[position - 1];
                List<Modification> fitting = mods.Where(m => m.Target.ToString() == "X" || m.Target.ToString() == residue.ToString()).ToList();
                if (fitting.Count == 0)
                {
                    continue;
                }
                var here = new List<Modification> { fitting[random.Next(fitting.Count)] };
                if (random.Next(4) == 0)
                {
                    Modification second = fitting[random.Next(fitting.Count)];
                    if (!here.Contains(second))
                    {
                        here.Add(second);
                    }
                }
                result.Add(position, here);
            }
            return result;
        }

        private static List<Modification> ProteinModifications() => new()
        {
            NewModification("Phospho", "S", "Anywhere.", 79.966331),
            NewModification("Phospho", "T", "Anywhere.", 79.966331),
            NewModification("Phospho", "Y", "Anywhere.", 79.966331),
            NewModification("Acetyl", "K", "Anywhere.", 42.010565),
            NewModification("Methyl", "K", "Anywhere.", 14.01565),
            NewModification("Methyl", "R", "Anywhere.", 14.01565),
            NewModification("Oxidation", "M", "Anywhere.", 15.994915),
            NewModification("Deamidation", "N", "Anywhere.", 0.984016),
            NewModification("Acetyl", "X", "N-terminal.", 42.010565),
            NewModification("Amidation", "X", "C-terminal.", -0.984016),
            NewModification("Formyl", "X", "Peptide N-terminal.", 27.994915),
        };

        private static List<Modification> RnaModifications() => new()
        {
            NewModification("Methyl", "A", "Anywhere.", 14.01565),
            NewModification("Methyl", "C", "Anywhere.", 14.01565),
            NewModification("Methyl", "G", "Anywhere.", 14.01565),
            NewModification("Pseudouridine", "U", "Anywhere.", 0),
            NewModification("Thio", "U", "Anywhere.", 15.977156),
        };

        private static Modification NewModification(string id, string motif, string location, double mass)
        {
            ModificationMotif.TryGetMotif(motif, out ModificationMotif target);
            return new Modification(_originalId: id, _modificationType: "Test", _target: target, _locationRestriction: location,
                _monoisotopicMass: mass);
        }

        #endregion
    }
}
