using System.Collections.Generic;
using System.Linq;
using Chemistry;
using EngineLayer.Truncation;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    /// <summary>
    /// Phase 3 tests: reverse-decoy generation for parents (01_Architecture.md decision #14), reusing
    /// the existing mzLib reversal so PTMs lock to residues, and not duplicating decoys already present.
    /// </summary>
    [TestFixture]
    public class TruncationDecoyTests
    {
        [Test]
        public void ReverseDecoy_Generated_PtmsLockToResidues()
        {
            // PEPT(phospho)IDE -> EDIT(phospho)PEP : the phospho stays with the T.
            var phospho = MakePhospho();
            var target = new TruncationParent(
                MakeProteoform("PEPTIDE", "ACC", new Dictionary<int, Modification> { { 5, phospho } }), // T is residue 4 -> key 5
                "ACC", "row", isDecoy: false);

            List<TruncationParent> combined = TruncationParentBuilder.AddReverseDecoys(new[] { target });

            TruncationParent decoy = combined.Single(p => p.IsDecoy);
            Assert.That(decoy.ProteinAccession, Is.EqualTo("DECOY_ACC"));
            Assert.That(decoy.Proteoform.BaseSequence, Is.EqualTo("EDITPEP"));
            var mod = decoy.Proteoform.AllModsOneIsNterminus.Single();
            Assert.That(mod.Key, Is.EqualTo(5)); // T sits at decoy residue 4 -> key 5
            Assert.That(mod.Value.OriginalId, Is.EqualTo("phospho"));
            Assert.That(decoy.Proteoform.BaseSequence[mod.Key - 2], Is.EqualTo('T')); // phospho still on a T
        }

        [Test]
        public void ReverseDecoy_NotDuplicated_WhenPass1Provided()
        {
            var target = new TruncationParent(MakeProteoform("PEPTIDE", "ACC", null), "ACC", "row", isDecoy: false);
            var existingDecoy = new TruncationParent(MakeProteoform("EDITPEP", "DECOY_ACC", null), "DECOY_ACC", "row", isDecoy: true);

            List<TruncationParent> combined = TruncationParentBuilder.AddReverseDecoys(new[] { target, existingDecoy });

            Assert.That(combined.Count, Is.EqualTo(2));                       // no fresh decoy generated
            Assert.That(combined.Count(p => p.IsDecoy), Is.EqualTo(1));
        }

        // ---------- helpers ----------

        private static Modification MakePhospho()
        {
            ModificationMotif.TryGetMotif("T", out var motif);
            return new Modification(_originalId: "phospho", _modificationType: "testMod", _target: motif,
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _locationRestriction: "Anywhere.");
        }

        private static PeptideWithSetModifications MakeProteoform(string sequence, string accession, Dictionary<int, Modification> mods)
        {
            var protein = new Protein(sequence, accession);
            var digestionParams = new DigestionParams(protease: "top-down", minPeptideLength: 1, maxPeptideLength: 100000);
            return new PeptideWithSetModifications(protein, digestionParams, 1, sequence.Length,
                CleavageSpecificity.Full, "top-down", 0, mods ?? new Dictionary<int, Modification>(), mods?.Count ?? 0);
        }
    }
}
