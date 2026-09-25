using Chemistry;
using EngineLayer;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Input;
using Transcriptomics;

namespace GuiFunctions
{
    public class CustomResidueViewModel : BaseViewModel, IDisposable
    {
        private string _name = "";
        private string _oneLetterCode = "";
        private string _symbol = "";
        private string _chemicalFormula = "";
        private string _validationMessage = "";
        private bool _isRnaMode;

        public CustomResidueViewModel()
        {
            SaveCommand = new RelayCommand(Save);
            CancelCommand = new RelayCommand(Cancel);

            _isRnaMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
            GuiGlobalParamsViewModel.Instance.PropertyChanged += GlobalParametersChanged;
            RefreshModeProperties();
            Validate();
        }

        public event EventHandler<CustomResidueDialogResultEventArgs> RequestClose;

        public ICommand SaveCommand { get; }
        public ICommand CancelCommand { get; }

        public string Name
        {
            get => _name;
            set
            {
                if (_name == value)
                    return;

                _name = value;
                OnPropertyChanged(nameof(Name));
                Validate();
            }
        }

        public string OneLetterCode
        {
            get => _oneLetterCode;
            set
            {
                if (_oneLetterCode == value)
                    return;

                _oneLetterCode = value;
                OnPropertyChanged(nameof(OneLetterCode));
                Validate();
            }
        }

        public string Symbol
        {
            get => _symbol;
            set
            {
                if (_symbol == value)
                    return;

                _symbol = value;
                OnPropertyChanged(nameof(Symbol));
                Validate();
            }
        }

        public string ChemicalFormula
        {
            get => _chemicalFormula;
            set
            {
                if (_chemicalFormula == value)
                    return;

                _chemicalFormula = value;
                OnPropertyChanged(nameof(ChemicalFormula));
                Validate();
            }
        }

        public bool IsRnaMode
        {
            get => _isRnaMode;
            private set
            {
                if (_isRnaMode == value)
                    return;

                _isRnaMode = value;
                OnPropertyChanged(nameof(IsRnaMode));
                OnPropertyChanged(nameof(IsSymbolVisible));
                Validate();
            }
        }

        public bool IsSymbolVisible => IsRnaMode;

        public bool CanSave => string.IsNullOrEmpty(ValidationMessage);

        public string ValidationMessage
        {
            get => _validationMessage;
            private set
            {
                if (_validationMessage == value)
                    return;

                _validationMessage = value;
                OnPropertyChanged(nameof(ValidationMessage));
                OnPropertyChanged(nameof(CanSave));
            }
        }

        private void GlobalParametersChanged(object sender, System.ComponentModel.PropertyChangedEventArgs e)
        {
            if (e.PropertyName == nameof(GuiGlobalParamsViewModel.IsRnaMode))
                IsRnaMode = GuiGlobalParamsViewModel.Instance.IsRnaMode;
        }

        private void RefreshModeProperties()
        {
            OnPropertyChanged(nameof(IsRnaMode));
            OnPropertyChanged(nameof(IsSymbolVisible));
        }

        private void Validate()
        {
            if (string.IsNullOrWhiteSpace(Name))
            {
                ValidationMessage = "A residue name is required.";
                return;
            }

            if (OneLetterCode.Length != 1)
            {
                ValidationMessage = "The one-letter abbreviation must contain exactly one character.";
                return;
            }

            if (IsRnaMode && string.IsNullOrWhiteSpace(Symbol))
            {
                ValidationMessage = "A three-letter nucleotide symbol is required.";
                return;
            }

            if (!TryParseFormula(out _))
            {
                ValidationMessage = "The chemical formula could not be parsed.";
                return;
            }

            char letter = OneLetterCode[0];
            if (IsRnaMode)
            {
                if (Nucleotide.TryGetResidue(letter, out _))
                {
                    ValidationMessage = $"The nucleotide letter '{letter}' already exists.";
                    return;
                }

                if (Nucleotide.TryGetResidue(Symbol, out _))
                {
                    ValidationMessage = $"The nucleotide symbol '{Symbol}' already exists.";
                    return;
                }

                if (Nucleotide.TryGetResidue(Name, out _))
                {
                    ValidationMessage = $"The nucleotide name '{Name}' already exists.";
                    return;
                }
            }
            else
            {
                if (GlobalVariables.InvalidAminoAcids.Contains(letter))
                {
                    ValidationMessage = $"The amino acid character '{letter}' is reserved and cannot be assigned.";
                    return;
                }

                if (Residue.TryGetResidue(letter, out _))
                {
                    ValidationMessage = $"The amino acid letter '{letter}' already exists.";
                    return;
                }
            }

            ValidationMessage = "";
        }

        private bool TryParseFormula(out ChemicalFormula formula)
        {
            try
            {
                formula = Chemistry.ChemicalFormula.ParseFormula(ChemicalFormula);
                return true;
            }
            catch
            {
                formula = null;
                return false;
            }
        }

        private void Save()
        {
            Validate();
            if (!CanSave || !TryParseFormula(out ChemicalFormula formula))
                return;

            try
            {
                if (IsRnaMode)
                    SaveNucleotide(formula);
                else
                    SaveAminoAcid(formula);

                RequestClose?.Invoke(this, new CustomResidueDialogResultEventArgs(true));
            }
            catch (Exception e)
            {
                ValidationMessage = $"Could not save the custom residue: {e.Message}";
            }
        }

        private void SaveAminoAcid(ChemicalFormula formula)
        {
            string directory = Path.Combine(GlobalVariables.DataDir, "CustomAminoAcids");
            string path = Path.Combine(directory, "CustomAminoAcids.txt");
            if (!File.Exists(path))
                GlobalVariables.WriteAminoAcidsFile();

            List<string> lines = File.ReadAllLines(path).ToList();
            lines.Add($"{Name}\t{OneLetterCode[0]}\t{formula.MonoisotopicMass}\t{formula.Formula}");
            File.WriteAllLines(path, lines);

            Residue.AddNewResiduesToDictionary(new List<Residue>
            {
                new Residue(Name, OneLetterCode[0], Name, formula, ModificationSites.Any)
            });
        }

        private void SaveNucleotide(ChemicalFormula formula)
        {
            string directory = Path.Combine(GlobalVariables.DataDir, "CustomNucleotides");
            string path = Path.Combine(directory, "CustomNucleotides.txt");
            if (!File.Exists(path))
                GlobalVariables.WriteNucleotidesFile();

            List<string> lines = File.ReadAllLines(path).ToList();
            lines.Add($"{Name}\t{OneLetterCode[0]}\t{Symbol}\t{formula.Formula}");
            File.WriteAllLines(path, lines);

            Nucleotide.AddResidue(Name, OneLetterCode[0], Symbol, formula);
        }

        private void Cancel()
        {
            RequestClose?.Invoke(this, new CustomResidueDialogResultEventArgs(false));
        }

        public void Dispose()
        {
            GuiGlobalParamsViewModel.Instance.PropertyChanged -= GlobalParametersChanged;
        }
    }

    public sealed class CustomResidueDialogResultEventArgs : EventArgs
    {
        public CustomResidueDialogResultEventArgs(bool succeeded)
        {
            Succeeded = succeeded;
        }

        public bool Succeeded { get; }
    }
}
