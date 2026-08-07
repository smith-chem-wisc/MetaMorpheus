using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System.Linq;
using Omics.Modifications;
using System.Collections.Generic;

namespace Test
{
    public class TestGlycan
    {
        /// <summary>
        /// Builds a small, valid N-glycan Kind[] (2 HexNAc, 3 Hex - a minimal trimannosyl-core-like
        /// composition) so RegenerateIons has real monosaccharides to work with.
        /// </summary>
        private static byte[] BuildTestKind()
        {
            byte[] kind = new byte[Glycan.KindCapacity];
            kind[0] = 3; // Hex
            kind[1] = 2; // HexNAc
            return kind;
        }

        [Test]
        public static void RegenerateIons_PopulatesIonsWhenNull()
        {
            var kind = BuildTestKind();
            var nGlycan = new Glycan(kind, "Nxs", GlycanType.N_glycan); // Ions == null (ToGenerateIons: false path)
            Assert.That(nGlycan.Ions, Is.Null); // Ions are only populated when RegenerateIons is called, or when the constructor is called with ToGenerateIons: true.
            Assert.That(nGlycan.DiagnosticIons, Is.Not.Null); // DiagnosticIons are always populated in the constructor, even when ToGenerateIons is false.

            nGlycan.RegenerateIons();

            Assert.That(nGlycan.Ions, Is.Not.Null);
            Assert.That(nGlycan.Ions.Count, Is.GreaterThan(0));
            Assert.That(nGlycan.Type, Is.EqualTo(GlycanType.N_glycan)); 

            var oGlycan = new Glycan(kind, "S", GlycanType.O_glycan);
            Assert.That(oGlycan.Ions, Is.Null);
            Assert.That(oGlycan.DiagnosticIons, Is.Not.Null);

            oGlycan.RegenerateIons();

            Assert.That(oGlycan.Ions, Is.Not.Null);
            Assert.That(oGlycan.Ions.Count, Is.GreaterThan(0));
            Assert.That(oGlycan.Type, Is.EqualTo(GlycanType.O_glycan));
        }

        [Test]
        public static void RegenerateIons_IsNoOpWhenAlreadyHydrated()
        {
            var kind = BuildTestKind();
            var glycan = new Glycan(kind, "Nxs", GlycanType.N_glycan);

            glycan.RegenerateIons();

            var originalIons = glycan.Ions;

            // Calling again should not rebuild/replace the already-populated Ions list, and return without doing anything.
            glycan.RegenerateIons();
            Assert.That(glycan.Ions, Is.SameAs(originalIons));
        }

        [Test]
        public static void RegenerateIons_DatabaseLoadedGlycanWithIons_ExpectedBehavior()
        {
            // Testing the behavior of RegenerateIons on a glycan that was loaded with ToGenerateIons: true.
            var nGlycanPath = GlobalVariables.NGlycanDatabasePaths.First(p => p.Contains("NGlycan.gdb"));
            var H1N2_withIon = GlycanDatabase.LoadGlycan(nGlycanPath, true, false).First(p => p.IdWithMotif == "H1N2 on Nxt");
            var H1N2_withoutIon = GlycanDatabase.LoadGlycan(nGlycanPath, false, false).First(p => p.IdWithMotif == "H1N2 on Nxt");

            // The glycan loaded with ToGenerateIons: false should have null Ions.
            Assert.That(H1N2_withoutIon, Is.Not.Null);
            Assert.That(H1N2_withoutIon.Ions, Is.Null);

            // The glycan loaded with ToGenerateIons: true should have non-null Ions and populated.
            Assert.That(H1N2_withIon, Is.Not.Null);
            Assert.That(H1N2_withIon.Ions, Is.Not.Null);
            var expectedGlycanIons = H1N2_withIon.Ions;

            H1N2_withoutIon.RegenerateIons();

            // After calling RegenerateIons on the glycan without Ions, it should now have non-null Ions, matching the behavior of the glycan that was loaded with ToGenerateIons: true.
            Assert.That(H1N2_withoutIon.Ions, Is.Not.Null);
            Assert.That(H1N2_withoutIon.Ions.Count, Is.GreaterThan(0));
            Assert.That(H1N2_withoutIon.Ions.Count, Is.EqualTo(expectedGlycanIons.Count));
            Assert.That(H1N2_withoutIon.Ions.Select(i => i.IonMass), Is.EqualTo(expectedGlycanIons.Select(i => i.IonMass)));
        }

        [Test]
        public static void NonGlycanNLinkedMod_IsSkippedByOfTypeGuard_NoCrash()
        {
            // Guards against a future protein-XML / GPTMD annotation carrying the
            // "N-linked glycosylation" ModificationType string on a plain Modification
            // (not a Glycan). At the consumption site such a mod must be filtered out by
            // OfType<Glycan>() rather than cast to null and NRE in GetGlycanYIons.
            ModificationMotif.TryGetMotif("Nxs", out var motif);
            var fakeNLinkedMod = new Modification(
                _originalId: "FakeN",
                _modificationType: "N-linked glycosylation",  // same string the consumption site keys on
                _target: motif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 203.079);                   // not a Glycan instance

            var allMods = new Dictionary<int, Modification> { { 1, fakeNLinkedMod } };

            // The consumption-site selection: OfType<Glycan>() must skip the plain Modification.
            var selectedGlycans = allMods.Values
                .OfType<Glycan>()
                .Where(g => g.ModificationType == "N-linked glycosylation")
                .ToList();

            Assert.That(selectedGlycans, Is.Empty,
                "A non-Glycan N-linked mod must be filtered out by OfType<Glycan>(), not cast to null.");

            // And it must not throw when the consumption loop runs over the (empty) selection.
            Assert.DoesNotThrow(() =>
            {
                foreach (var g in selectedGlycans)
                {
                    g.RegenerateIons();
                    _ = GlycoPeptides.GetGlycanYIons(1000.0, g);
                }
            });
        }

    }
}
