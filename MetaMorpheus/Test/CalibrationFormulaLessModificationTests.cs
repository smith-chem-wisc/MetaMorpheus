using Chemistry;
using EngineLayer;
using EngineLayer.Calibration;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Calibration and modifications with a mass but no chemical formula (issue #1875). The isobaric labels in Mods\TMT.txt
    /// are such modifications. Fixed on every peptide, they used to leave calibration only the few matches without a label.
    /// </summary>
    [TestFixture]
    internal class CalibrationFormulaLessModificationTests
    {
        private static Modification Tmt10(string residue, string location, bool withFormula)
        {
            ModificationMotif.TryGetMotif(residue, out var motif);
            // Unimod 737: TMT6/10/11, C8 13C4 H20 N 15N O2, 229.162932. Mods\TMT.txt gives the mass only.
            return withFormula
                ? new Modification(_originalId: "TMT10", _modificationType: "Formula Test", _target: motif, _locationRestriction: location,
                    _chemicalFormula: ChemicalFormula.ParseFormula("C8C{13}4H20NN{15}O2"))
                : new Modification(_originalId: "TMT10", _modificationType: "Formula Test", _target: motif, _locationRestriction: location,
                    _monoisotopicMass: 229.162932);
        }

        private static PeptideWithSetModifications LabelledPeptide(bool withFormula)
        {
            var nTerm = Tmt10("X", "Peptide N-terminal.", withFormula);
            var lys = Tmt10("K", "Anywhere.", withFormula);
            var mods = new Dictionary<string, Modification> { [nTerm.IdWithMotif] = nTerm, [lys.IdWithMotif] = lys };
            return new PeptideWithSetModifications("[Formula Test:TMT10 on X]LGEHNIDVLEGNEQFINAAK[Formula Test:TMT10 on K]", mods);
        }

        [Test]
        public static void AMatchWhoseModificationsAllHaveFormulasGetsTheEnvelopeItAlwaysGot()
        {
            var peptide = new PeptideWithSetModifications("LGEHNIDVLEGNEQFINAAK", GlobalVariables.AllModsKnownDictionary);
            var expected = IsotopicDistribution.GetDistribution(peptide.ThisChemicalFormula, 0.1, 0.001);

            var (masses, intensities) = DataPointAcquisitionEngine.TheoreticalIsotopicDistribution(peptide);

            Assert.That(masses, Is.EqualTo(expected.Masses.ToArray()));
            Assert.That(intensities, Is.EqualTo(expected.Intensities.ToArray()));
        }

        [Test]
        public static void AFormulaLessLabelGetsAnEnvelopeAnchoredAtItsExactMassAndCloseToTheTrueOne()
        {
            var formulaLess = LabelledPeptide(withFormula: false);
            var exact = LabelledPeptide(withFormula: true);
            Assert.That(formulaLess.ThisChemicalFormula, Is.Null, "the fixture must carry a formula-less modification");
            Assert.That(formulaLess.MonoisotopicMass, Is.EqualTo(exact.MonoisotopicMass).Within(1e-6));

            var (masses, intensities) = DataPointAcquisitionEngine.TheoreticalIsotopicDistribution(formulaLess);
            var truth = IsotopicDistribution.GetDistribution(exact.ThisChemicalFormula, 0.1, 0.001);

            Assert.That(masses, Is.Not.Empty);
            Assert.That(masses[0], Is.EqualTo(formulaLess.MonoisotopicMass).Within(1e-6), "anchored at the exact monoisotopic mass");
            // The three most abundant isotopologues, as calibration searches them: same peaks, within 2 mDa and 5% of the
            // tallest. Averagine stands in for two labels' atoms only; the residues keep their own formula.
            for (int i = 0; i < 3; i++)
            {
                Assert.That(masses[i], Is.EqualTo(truth.Masses[i]).Within(0.002), $"isotopologue {i} mass");
                Assert.That(intensities[i] / intensities.Max(), Is.EqualTo(truth.Intensities[i] / truth.Intensities.Max()).Within(0.05), $"isotopologue {i} intensity");
            }
        }

        [Test]
        public static void TheKnownFormulaLeavesOutOnlyTheFormulaLessModifications()
        {
            var formulaLess = LabelledPeptide(withFormula: false);
            var unlabelled = new PeptideWithSetModifications("LGEHNIDVLEGNEQFINAAK", GlobalVariables.AllModsKnownDictionary);

            Assert.That(DataPointAcquisitionEngine.KnownChemicalFormula(formulaLess), Is.EqualTo(unlabelled.ThisChemicalFormula));
        }

        /// <summary>
        /// A label fixed on every peptide N-terminus, with a mass (here 0 Da, so the search finds what it finds without it) and
        /// no formula. Before #1875 was fixed, every match carried it, so calibration skipped every match and failed with no MS1
        /// point at all. Now the envelope is built for it, and with no mass to fill it is the exact envelope, so calibration
        /// writes what it writes without the label, to the byte.
        /// </summary>
        [Test]
        public static void AFormulaLessLabelOnEveryPeptideCalibratesAsThePeptidesWithoutIt()
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            var label = new Modification(_originalId: "ZeroLabel", _modificationType: "Formula Test", _target: motif,
                _locationRestriction: "Peptide N-terminal.", _monoisotopicMass: 0);
            GlobalVariables.AddMods(new List<Modification> { label }, false);

            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "CalibrationFormulaLessLabel");
            if (Directory.Exists(root)) Directory.Delete(root, true);
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML");
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");

            string Calibrate(string name, bool labelled)
            {
                var fixedMods = new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C") };
                if (labelled) fixedMods.Add(("Formula Test", "ZeroLabel on X"));
                var task = new CalibrationTask { CommonParameters = new CommonParameters(listOfModsFixed: fixedMods) };
                string output = Path.Combine(root, name);
                Directory.CreateDirectory(output);
                task.RunTask(output, new List<DbForTask> { new DbForTask(database, false) }, new List<string> { spectra }, "test");
                return output;
            }

            try
            {
                string plain = Calibrate("plain", labelled: false);
                string labelled = Calibrate("labelled", labelled: true);

                string calibrated = Path.Combine(labelled, "SmallCalibratible_Yeast-calib.mzML");
                Assert.That(File.Exists(calibrated), "calibration must succeed with a formula-less label on every peptide");
                Assert.That(File.ReadAllText(Path.Combine(labelled, "SmallCalibratible_Yeast-calib.toml")),
                    Is.EqualTo(File.ReadAllText(Path.Combine(plain, "SmallCalibratible_Yeast-calib.toml"))));
            }
            finally
            {
                if (Directory.Exists(root)) Directory.Delete(root, true);
            }
        }
    }
}
