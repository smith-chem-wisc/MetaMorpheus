using Chemistry;
using Proteomics;
using System.Collections.Generic;

namespace MetaMorpheusGUI
{
    internal class SilacInfoForDataGrid
    {
        public readonly List<SilacLabel> SilacLabel;
        public readonly List<ChemicalFormula> LabelFormula;
        public string AminoAcidLabels { get; set; } //need the {get; set;} for the string to be displayed
        public SilacModificationWindow.ExperimentType LabelType { get; set; } //need the {get; set;} for the string to be displayed

        public SilacInfoForDataGrid(SilacLabel label, SilacModificationWindow.ExperimentType labelType)
        {
            SilacLabel = new List<SilacLabel> { label };
            LabelFormula = new List<ChemicalFormula> { ChemicalFormula.ParseFormula(label.LabelChemicalFormula) };
            AminoAcidLabels = label.OriginalAminoAcid + "(" + label.MassDifference + ")";
            LabelType = labelType;
        }

        public SilacInfoForDataGrid(SilacModificationWindow.ExperimentType labelType)
        {
            AminoAcidLabels = "Unlabeled";
            LabelType = labelType;
        }

        public void AddAdditionalLabel(SilacInfoForDataGrid label)
        {
            SilacLabel.AddRange(label.SilacLabel);
            LabelFormula.AddRange(label.LabelFormula);
            AminoAcidLabels += " & " + label.AminoAcidLabels;
        }
    }
}
