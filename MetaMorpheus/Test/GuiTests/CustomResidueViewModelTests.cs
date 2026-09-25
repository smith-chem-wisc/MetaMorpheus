using Chemistry;
using EngineLayer;
using GuiFunctions;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.IO;
using System.Linq;
using Transcriptomics;

namespace Test.GuiTests
{
    [TestFixture]
    public class CustomResidueViewModelTests
    {
        [Test]
        public void StartsWithRequiredFieldValidation()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel();

                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage, Is.EqualTo("A residue name is required."));
                Assert.That(viewModel.IsSymbolVisible, Is.False);
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void ProteinModeRejectsInvalidFormulaAndExistingLetter()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Test amino acid",
                    OneLetterCode = "q",
                    ChemicalFormula = "not a formula"
                };

                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage, Is.EqualTo("The chemical formula could not be parsed."));

                viewModel.ChemicalFormula = "C2H3NO";
                Assert.That(viewModel.CanSave, Is.True);

                viewModel.OneLetterCode = "A";
                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage, Does.Contain("already exists"));
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void RnaModeRequiresAndValidatesNucleotideSymbol()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = true;

            try
            {
                char unusedLetter = GetUnusedNucleotideLetter();
                using var viewModel = new CustomResidueViewModel
                {
                    Name = $"Test nucleotide {unusedLetter}",
                    OneLetterCode = unusedLetter.ToString(),
                    ChemicalFormula = "C5H5N2O2"
                };

                Assert.That(viewModel.IsRnaMode, Is.True);
                Assert.That(viewModel.IsSymbolVisible, Is.True);
                Assert.That(viewModel.ValidationMessage, Is.EqualTo("A three-letter nucleotide symbol is required."));

                viewModel.Symbol = $"T{unusedLetter}t";
                Assert.That(viewModel.CanSave, Is.True);

                viewModel.Symbol = "Ade";
                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage, Does.Contain("already exists"));
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void ModeChangesUpdateSymbolVisibility()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel();
                Assert.That(viewModel.IsSymbolVisible, Is.False);

                GuiGlobalParamsViewModel.Instance.IsRnaMode = true;

                Assert.That(viewModel.IsRnaMode, Is.True);
                Assert.That(viewModel.IsSymbolVisible, Is.True);
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void CancelCommandRequestsFailedClose()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel();
                CustomResidueDialogResultEventArgs result = null;
                viewModel.RequestClose += (_, args) => result = args;

                viewModel.CancelCommand.Execute(null);

                Assert.That(result, Is.Not.Null);
                Assert.That(result.Succeeded, Is.False);
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SaveCommandRegistersAndPersistsAminoAcid()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
            string path = Path.Combine(GlobalVariables.DataDir, "CustomAminoAcids", "CustomAminoAcids.txt");
            bool hadFile = File.Exists(path);
            string[] originalLines = hadFile ? File.ReadAllLines(path) : Array.Empty<string>();

            try
            {
                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Test saved amino acid",
                    OneLetterCode = GetUnusedAminoAcidLetter().ToString(),
                    ChemicalFormula = "C2H3NO"
                };
                CustomResidueDialogResultEventArgs result = null;
                viewModel.RequestClose += (_, args) => result = args;

                viewModel.SaveCommand.Execute(null);

                Assert.That(result?.Succeeded, Is.True);
                char letter = viewModel.OneLetterCode[0];
                Assert.That(Residue.TryGetResidue(letter, out Residue residue), Is.True);
                Assert.That(residue.ThisChemicalFormula.Formula, Is.EqualTo("C2H3NO"));
                Assert.That(File.ReadAllLines(path).Any(line => line.StartsWith($"Test saved amino acid\t{letter}\t")), Is.True);
            }
            finally
            {
                if (hadFile)
                    File.WriteAllLines(path, originalLines);
                else if (File.Exists(path))
                    File.Delete(path);
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SaveCommandRegistersAndPersistsNucleotide()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
            string path = Path.Combine(GlobalVariables.DataDir, "CustomNucleotides", "CustomNucleotides.txt");
            bool hadFile = File.Exists(path);
            string[] originalLines = hadFile ? File.ReadAllLines(path) : Array.Empty<string>();

            try
            {
                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Test saved nucleotide",
                    OneLetterCode = GetUnusedNucleotideLetter().ToString(),
                    ChemicalFormula = "C5H5N2O2"
                };
                viewModel.Symbol = $"T{viewModel.OneLetterCode}v";
                CustomResidueDialogResultEventArgs result = null;
                viewModel.RequestClose += (_, args) => result = args;

                viewModel.SaveCommand.Execute(null);

                Assert.That(result?.Succeeded, Is.True);
                char letter = viewModel.OneLetterCode[0];
                Assert.That(Nucleotide.TryGetResidue(letter, out Nucleotide nucleotide), Is.True);
                Assert.That(nucleotide.Symbol, Is.EqualTo(viewModel.Symbol));
                Assert.That(File.ReadAllLines(path).Any(line => line.StartsWith($"Test saved nucleotide\t{letter}\t{viewModel.Symbol}\t")), Is.True);
            }
            finally
            {
                if (hadFile)
                    File.WriteAllLines(path, originalLines);
                else if (File.Exists(path))
                    File.Delete(path);
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SaveCommandCreatesAminoAcidsFileWhenMissing()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
            string path = Path.Combine(GlobalVariables.DataDir, "CustomAminoAcids", "CustomAminoAcids.txt");
            bool hadFile = File.Exists(path);
            string[] originalLines = hadFile ? File.ReadAllLines(path) : Array.Empty<string>();

            try
            {
                if (File.Exists(path))
                    File.Delete(path);

                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Test missing amino-acid file",
                    OneLetterCode = GetUnusedAminoAcidLetter().ToString(),
                    ChemicalFormula = "C2H3NO"
                };

                viewModel.SaveCommand.Execute(null);

                Assert.That(File.Exists(path), Is.True);
                char letter = viewModel.OneLetterCode[0];
                Assert.That(File.ReadAllLines(path)
                    .Any(line => line.StartsWith($"Test missing amino-acid file\t{letter}\t")), Is.True);
            }
            finally
            {
                if (hadFile)
                    File.WriteAllLines(path, originalLines);
                else if (File.Exists(path))
                    File.Delete(path);

                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SaveCommandCreatesNucleotidesFileWhenMissing()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
            string path = Path.Combine(GlobalVariables.DataDir, "CustomNucleotides", "CustomNucleotides.txt");
            bool hadFile = File.Exists(path);
            string[] originalLines = hadFile ? File.ReadAllLines(path) : Array.Empty<string>();

            try
            {
                if (File.Exists(path))
                    File.Delete(path);

                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Test missing nucleotide file",
                    OneLetterCode = GetUnusedNucleotideLetter().ToString(),
                    ChemicalFormula = "C5H5N2O2"
                };
                viewModel.Symbol = $"T{viewModel.OneLetterCode}m";

                viewModel.SaveCommand.Execute(null);

                Assert.That(File.Exists(path), Is.True);
                char letter = viewModel.OneLetterCode[0];
                Assert.That(File.ReadAllLines(path)
                    .Any(line => line.StartsWith($"Test missing nucleotide file\t{letter}\t{viewModel.Symbol}\t")), Is.True);
            }
            finally
            {
                if (hadFile)
                    File.WriteAllLines(path, originalLines);
                else if (File.Exists(path))
                    File.Delete(path);

                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        // ── Identical-set guard tests ─────────────────────────────────────────

        [Test]
        public void SettingNameToSameValueDoesNotFirePropertyChanged()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel { Name = "Same" };
                int eventCount = 0;
                viewModel.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(viewModel.Name)) eventCount++; };

                viewModel.Name = "Same";

                Assert.That(eventCount, Is.Zero);
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SettingOneLetterCodeToSameValueDoesNotFirePropertyChanged()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel { OneLetterCode = "q" };
                int eventCount = 0;
                viewModel.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(viewModel.OneLetterCode)) eventCount++; };

                viewModel.OneLetterCode = "q";

                Assert.That(eventCount, Is.Zero);
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SettingSymbolToSameValueDoesNotFirePropertyChanged()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = true;

            try
            {
                using var viewModel = new CustomResidueViewModel { Symbol = "Xyz" };
                int eventCount = 0;
                viewModel.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(viewModel.Symbol)) eventCount++; };

                viewModel.Symbol = "Xyz";

                Assert.That(eventCount, Is.Zero);
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SettingChemicalFormulaToSameValueDoesNotFirePropertyChanged()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel { ChemicalFormula = "C2H3NO" };
                int eventCount = 0;
                viewModel.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(viewModel.ChemicalFormula)) eventCount++; };

                viewModel.ChemicalFormula = "C2H3NO";

                Assert.That(eventCount, Is.Zero);
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        // ── Remaining validation-message branch tests ─────────────────────────

        [Test]
        public void OneLetterCodeWrongLengthFailsValidation()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Test residue",
                    OneLetterCode = "ab"
                };

                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage,
                    Is.EqualTo("The one-letter abbreviation must contain exactly one character."));
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void ProteinModeRejectsReservedAminoAcidCharacter()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

            try
            {
                char reserved = GlobalVariables.InvalidAminoAcids.First();
                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Reserved test",
                    OneLetterCode = reserved.ToString(),
                    ChemicalFormula = "C2H3NO"
                };

                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage,
                    Is.EqualTo($"The amino acid character '{reserved}' is reserved and cannot be assigned."));
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void RnaModeRejectsDuplicateNucleotideLetter()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = true;

            try
            {
                // 'A' maps to Adenine in the default nucleotide table
                char existingLetter = 'A';
                Assert.That(Nucleotide.TryGetResidue(existingLetter, out _), Is.True,
                    "Precondition: 'A' must be a known nucleotide letter.");

                string unusedSymbol = $"X{GetUnusedNucleotideLetter()}x";
                using var viewModel = new CustomResidueViewModel
                {
                    Name = "Duplicate letter nucleotide",
                    OneLetterCode = existingLetter.ToString(),
                    Symbol = unusedSymbol,
                    ChemicalFormula = "C5H5N2O2"
                };

                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage,
                    Is.EqualTo($"The nucleotide letter '{existingLetter}' already exists."));
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void RnaModeRejectsDuplicateNucleotideName()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = true;

            try
            {
                // "Adenine" is a standard nucleotide name
                const string existingName = "Adenine";
                Assert.That(Nucleotide.TryGetResidue(existingName, out _), Is.True,
                    "Precondition: 'Adenine' must be a known nucleotide name.");

                char unusedLetter = GetUnusedNucleotideLetter();
                string unusedSymbol = $"X{unusedLetter}x";
                using var viewModel = new CustomResidueViewModel
                {
                    Name = existingName,
                    OneLetterCode = unusedLetter.ToString(),
                    Symbol = unusedSymbol,
                    ChemicalFormula = "C5H5N2O2"
                };

                Assert.That(viewModel.CanSave, Is.False);
                Assert.That(viewModel.ValidationMessage,
                    Is.EqualTo($"The nucleotide name '{existingName}' already exists."));
            }
            finally
            {
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        private static char GetUnusedAminoAcidLetter() => Enumerable.Range('a', 'z' - 'a' + 1)
            .Select(value => (char)value)
            .First(letter => !Residue.TryGetResidue(letter, out _));

        private static char GetUnusedNucleotideLetter() => Enumerable.Range('a', 'z' - 'a' + 1)
            .Select(value => (char)value)
            .First(letter => !Nucleotide.TryGetResidue(letter, out _));
    }
}
