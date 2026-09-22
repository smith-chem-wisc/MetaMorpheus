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
            string[] originalLines = File.ReadAllLines(path);

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
                File.WriteAllLines(path, originalLines);
                GuiGlobalParamsViewModel.Instance.IsRnaMode = originalMode;
            }
        }

        [Test]
        public void SaveCommandRegistersAndPersistsNucleotide()
        {
            bool originalMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
            string path = Path.Combine(GlobalVariables.DataDir, "CustomNucleotides", "CustomNucleotides.txt");
            string[] originalLines = File.ReadAllLines(path);

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
                File.WriteAllLines(path, originalLines);
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
