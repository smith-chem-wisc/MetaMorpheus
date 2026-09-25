using EngineLayer;
using EngineLayer.Indexing;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// The peptide index is cached on disk and reused whenever <see cref="IndexingEngine.ToString"/>
    /// matches the params file written beside it (MetaMorpheusTask.SameSettings compares the two as
    /// literal text). These tests pin the property that makes reuse safe: anything that changes which
    /// peptides land in the index must change the key.
    ///
    /// The key used to carry modification *counts* only ("Number of variable mods: 1"). Two searches
    /// whose modification lists differed but happened to be the same length therefore shared a key, and
    /// the second silently reused an index built for the first one's modifications. Because the index
    /// stores decorated PeptideWithSetModifications, that index is wrong and not merely incomplete.
    /// </summary>
    [TestFixture]
    public static class IndexCacheKeyTest
    {
        private const char CarriageReturn = (char)13;
        private const char LineFeed = (char)10;

        private static Modification Mod(string id, string motifString, double mass, string locationRestriction = "Anywhere.")
        {
            ModificationMotif.TryGetMotif(motifString, out ModificationMotif motif);
            return new Modification(_originalId: id, _modificationType: "Test", _target: motif,
                _locationRestriction: locationRestriction, _monoisotopicMass: mass);
        }

        private static Modification Oxidation => Mod("Oxidation on M", "M", 15.994915);
        private static Modification Phospho => Mod("Phospho on S", "S", 79.966331);
        private static Modification Carbamidomethyl => Mod("Carbamidomethyl on C", "C", 57.021464);
        private static Modification Iodoacetamide => Mod("Iodoacetamide derivative on C", "C", 58.005479);

        private static string KeyFor(List<Modification> variableMods = null, List<Modification> fixedMods = null,
            List<SilacLabel> silacLabels = null, SilacLabel startLabel = null, SilacLabel endLabel = null)
        {
            CommonParameters commonParameters = new CommonParameters();
            var fileSpecificParameters = new List<(string, CommonParameters)> { ("", commonParameters) };

            var engine = new IndexingEngine(
                new List<Protein> { new Protein("MNNNKQQQSTC", "accession") },
                variableMods ?? new List<Modification>(),
                fixedMods ?? new List<Modification>(),
                silacLabels, startLabel, endLabel, 1, DecoyType.None, commonParameters, fileSpecificParameters,
                30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            return engine.ToString();
        }

        private static string LineStartingWith(string key, string prefix) =>
            key.Split(LineFeed).Select(line => line.TrimEnd(CarriageReturn)).Single(line => line.StartsWith(prefix));

        // ----- the defect this change exists to fix -----

        [Test]
        public static void SwappingOneVariableModForAnotherChangesTheKey()
        {
            string withOxidation = KeyFor(variableMods: new List<Modification> { Oxidation });
            string withPhospho = KeyFor(variableMods: new List<Modification> { Phospho });

            Assert.That(withPhospho, Is.Not.EqualTo(withOxidation),
                "Same number of variable mods, different mods. A shared key here reuses an index that " +
                "contains no phospho peptides at all.");
        }

        [Test]
        public static void SwappingOneFixedModForAnotherChangesTheKey()
        {
            string withCarbamidomethyl = KeyFor(fixedMods: new List<Modification> { Carbamidomethyl });
            string withIodoacetamide = KeyFor(fixedMods: new List<Modification> { Iodoacetamide });

            Assert.That(withIodoacetamide, Is.Not.EqualTo(withCarbamidomethyl),
                "Same number of fixed mods, different mods. A shared key here reuses an index whose " +
                "every peptide mass is wrong.");
        }

        [Test]
        public static void SwappingModsWithinAMultiModListChangesTheKey()
        {
            string first = KeyFor(variableMods: new List<Modification> { Oxidation, Phospho });
            string second = KeyFor(variableMods: new List<Modification> { Oxidation, Mod("Acetyl on K", "K", 42.010565) });

            Assert.That(second, Is.Not.EqualTo(first));
        }

        [Test]
        public static void TheKeyNoLongerCountsModifications()
        {
            string key = KeyFor(variableMods: new List<Modification> { Oxidation });

            Assert.That(key, Does.Not.Contain("Number of variable mods"));
            Assert.That(key, Does.Not.Contain("Number of fixed mods"));
        }

        // ----- reuse must still happen when nothing relevant changed -----

        [Test]
        public static void IdenticalModificationsGiveTheSameKey()
        {
            string first = KeyFor(variableMods: new List<Modification> { Oxidation, Phospho },
                fixedMods: new List<Modification> { Carbamidomethyl });
            string second = KeyFor(variableMods: new List<Modification> { Oxidation, Phospho },
                fixedMods: new List<Modification> { Carbamidomethyl });

            Assert.That(second, Is.EqualTo(first),
                "Equal inputs must reuse the cached index, or every search re-indexes for nothing.");
        }

        [Test]
        public static void TheKeyIsStableAcrossRepeatedCalls()
        {
            string key = KeyFor(variableMods: new List<Modification> { Oxidation });

            Assert.That(KeyFor(variableMods: new List<Modification> { Oxidation }), Is.EqualTo(key));
        }

        // ----- identity is the whole definition, not just the name -----

        [Test]
        public static void EditingAModificationsMassInPlaceChangesTheKey()
        {
            // A user editing a custom mods file, keeping the name and changing the mass. Keying on the
            // id alone would reuse an index built with the old mass.
            string original = KeyFor(variableMods: new List<Modification> { Mod("Custom on M", "M", 15.994915) });
            string edited = KeyFor(variableMods: new List<Modification> { Mod("Custom on M", "M", 16.994915) });

            Assert.That(edited, Is.Not.EqualTo(original));
        }

        [Test]
        public static void ChangingOnlyTheLocationRestrictionChangesTheKey()
        {
            string anywhere = KeyFor(variableMods: new List<Modification> { Mod("Acetyl on K", "K", 42.010565, "Anywhere.") });
            string nTerminal = KeyFor(variableMods: new List<Modification> { Mod("Acetyl on K", "K", 42.010565, "N-terminal.") });

            Assert.That(nTerminal, Is.Not.EqualTo(anywhere),
                "Same name and mass, different placement, so different peptides.");
        }

        [Test]
        public static void TheKeyNamesEachModificationAndQualifiesItWithAShortHash()
        {
            string line = LineStartingWith(KeyFor(variableMods: new List<Modification> { Oxidation, Phospho }), "Variable mods: ");

            Assert.That(line, Is.EqualTo("Variable mods: " +
                "Oxidation on M[" + IndexingEngine.ShortHash(Oxidation.ToString()) + "]," +
                "Phospho on S[" + IndexingEngine.ShortHash(Phospho.ToString()) + "]"),
                "Readable id first, so a human can diff two params files and see what changed.");
        }

        // ----- order is deliberately significant -----

        [Test]
        public static void ReorderingTheVariableModListChangesTheKey()
        {
            // Not cosmetic. GetVariableModificationPatterns is truncated at MaxModificationIsoforms by a
            // yield break, so which isoforms survive can depend on the order the mods arrive in. Sorting
            // the key would let two orders reuse each other's index.
            string forward = KeyFor(variableMods: new List<Modification> { Oxidation, Phospho });
            string reversed = KeyFor(variableMods: new List<Modification> { Phospho, Oxidation });

            Assert.That(reversed, Is.Not.EqualTo(forward));
        }

        // ----- SILAC labels also decide the peptides, and were absent from the key entirely -----

        [Test]
        public static void SwappingSilacLabelsChangesTheKey()
        {
            string lysine = KeyFor(silacLabels: new List<SilacLabel> { new SilacLabel('K', 'a', "C{13}6", 6.020129) });
            string arginine = KeyFor(silacLabels: new List<SilacLabel> { new SilacLabel('R', 'b', "C{13}6", 6.020129) });

            Assert.That(arginine, Is.Not.EqualTo(lysine));
        }

        [Test]
        public static void ChangingASilacLabelMassChangesTheKey()
        {
            string light = KeyFor(silacLabels: new List<SilacLabel> { new SilacLabel('K', 'a', "C{13}6", 6.020129) });
            string heavy = KeyFor(silacLabels: new List<SilacLabel> { new SilacLabel('K', 'a', "C{13}6", 8.014199) });

            Assert.That(heavy, Is.Not.EqualTo(light));
        }

        [Test]
        public static void AnAdditionalSilacLabelChangesTheKey()
        {
            SilacLabel plain = new SilacLabel('K', 'a', "C{13}6", 6.020129);

            SilacLabel compound = new SilacLabel('K', 'a', "C{13}6", 6.020129);
            compound.AddAdditionalSilacLabel(new SilacLabel('R', 'b', "C{13}6", 6.020129));

            Assert.That(KeyFor(silacLabels: new List<SilacLabel> { compound }),
                Is.Not.EqualTo(KeyFor(silacLabels: new List<SilacLabel> { plain })));
        }

        [Test]
        public static void TurnoverLabelsChangeTheKey()
        {
            SilacLabel start = new SilacLabel('K', 'a', "C{13}6", 6.020129);
            SilacLabel end = new SilacLabel('K', 'b', "N{15}2", 2.004);

            string none = KeyFor();
            string withTurnover = KeyFor(startLabel: start, endLabel: end);
            string swapped = KeyFor(startLabel: end, endLabel: start);

            Assert.That(withTurnover, Is.Not.EqualTo(none));
            Assert.That(swapped, Is.Not.EqualTo(withTurnover), "Start and end are not interchangeable.");
        }

        [Test]
        public static void MultipleSilacLabelsAreSeparated()
        {
            // Pins the separator. Without it two labels run together into one token, and a pair of lists
            // that differ only in where one label ends and the next begins would share a key.
            SilacLabel first = new SilacLabel('K', 'a', "C{13}6", 6.020129);
            SilacLabel second = new SilacLabel('R', 'b', "C{13}6", 6.020129);

            string line = LineStartingWith(KeyFor(silacLabels: new List<SilacLabel> { first, second }), "Silac labels: ");

            Assert.That(line, Is.EqualTo("Silac labels: " +
                IndexingEngine.DescribeSilacLabel(first) + "," + IndexingEngine.DescribeSilacLabel(second)));
        }

        [Test]
        public static void NoTurnoverLabelsRendersAsNoneRatherThanAPair()
        {
            // The constructor only builds TurnoverLabels when at least one label is given. If that guard
            // inverted, "no labels" would render as a pair of empties and stop being distinguishable from
            // a real pair that happens to be empty.
            string line = LineStartingWith(KeyFor(), "Turnover labels: ");

            Assert.That(line, Is.EqualTo("Turnover labels: none"));
        }

        [Test]
        public static void OneTurnoverLabelIsDistinguishableFromNone()
        {
            SilacLabel start = new SilacLabel('K', 'a', "C{13}6", 6.020129);

            Assert.That(KeyFor(startLabel: start), Is.Not.EqualTo(KeyFor()));
        }

        // ----- degenerate inputs must be described, not thrown on -----

        [Test]
        public static void NullAndEmptyModificationListsAreDescribedAndDistinct()
        {
            Assert.That(IndexingEngine.DescribeModifications(null), Is.EqualTo("none"));
            Assert.That(IndexingEngine.DescribeModifications(new List<Modification>()), Is.Empty);
        }

        [Test]
        public static void NullSilacLabelsAreDescribedAndDistinct()
        {
            Assert.That(IndexingEngine.DescribeSilacLabels(null), Is.EqualTo("none"));
            Assert.That(IndexingEngine.DescribeSilacLabels(new List<SilacLabel>()), Is.Empty);
            Assert.That(IndexingEngine.DescribeSilacLabel(null), Is.EqualTo("none"));
        }

        [Test]
        public static void AModificationWithNoIdIsStillDescribed()
        {
            // Modification.IdWithMotif is null when the mod was built without an id; ValidModification
            // rejects those, but the key must not throw on the way to finding that out.
            string described = IndexingEngine.DescribeModifications(new List<Modification> { new Modification() });

            Assert.That(described, Does.StartWith("unnamed["));
        }

        [Test]
        public static void ANullEntryInTheModificationListIsStillDescribed()
        {
            Assert.That(() => IndexingEngine.DescribeModifications(new List<Modification> { null }), Throws.Nothing);
        }

        // ----- the short hash itself -----

        [Test]
        public static void ShortHashIgnoresLineEndings()
        {
            // Modification.ToString builds its text with AppendLine, which emits the writing machine's
            // line ending. Spelled out as char codes so the test cannot be confused by its own source
            // file's line endings.
            string windowsStyle = "MM   15.994915" + CarriageReturn + LineFeed + "TG   M";
            string unixStyle = "MM   15.994915" + LineFeed + "TG   M";
            string oldMacStyle = "MM   15.994915" + CarriageReturn + "TG   M";

            Assert.That(IndexingEngine.ShortHash(unixStyle), Is.EqualTo(IndexingEngine.ShortHash(windowsStyle)));
            Assert.That(IndexingEngine.ShortHash(oldMacStyle), Is.EqualTo(IndexingEngine.ShortHash(windowsStyle)));
        }

        [Test]
        public static void ShortHashSeparatesDifferentDefinitions()
        {
            Assert.That(IndexingEngine.ShortHash("MM   15.994915"), Is.Not.EqualTo(IndexingEngine.ShortHash("MM   15.994916")));
        }

        [Test]
        public static void ShortHashIsSixteenLowercaseHexCharacters()
        {
            string shortHash = IndexingEngine.ShortHash(Oxidation.ToString());

            Assert.That(shortHash, Has.Length.EqualTo(16));
            Assert.That(shortHash.All(c => (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')), Is.True,
                "Expected lowercase hex, got: " + shortHash);
        }

        [Test]
        public static void ShortHashOfNothingIsNone()
        {
            Assert.That(IndexingEngine.ShortHash(null), Is.EqualTo("none"));
            Assert.That(IndexingEngine.ShortHash(string.Empty), Is.EqualTo("none"));
        }

        // ----- guard the lines that were already correct -----

        [Test]
        public static void AddingAModificationStillChangesTheKey()
        {
            string one = KeyFor(variableMods: new List<Modification> { Oxidation });
            string two = KeyFor(variableMods: new List<Modification> { Oxidation, Phospho });

            Assert.That(two, Is.Not.EqualTo(one));
        }
    }
}
