using Chemistry;
using Proteomics;
using System.Collections.Generic;

namespace MetaMorpheusGUI
{
    internal class SilacInfoForDataGrid
    {
        public readonly List<SilacLabel> SilacLabel;
        public readonly List<ChemicalFormula> LabelFormula;
        public string Label { get; set; } //need the {get; set;} for the string to be displayed

        public SilacInfoForDataGrid(SilacLabel label)
        {
            SilacLabel = new List<SilacLabel> { label };
            LabelFormula = new List<ChemicalFormula> { ChemicalFormula.ParseFormula(label.LabelChemicalFormula) };
            Label = label.OriginalAminoAcid + "(" + label.MassDifference + ")";
        }

        public void AddAdditionalLabel(SilacInfoForDataGrid label)
        {
            SilacLabel.AddRange(label.SilacLabel);
            LabelFormula.AddRange(label.LabelFormula);
            Label += " & " + label.Label;
        }
    }
}
