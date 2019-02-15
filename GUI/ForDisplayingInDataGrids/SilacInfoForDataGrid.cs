using Chemistry;
using Proteomics;

namespace MetaMorpheusGUI
{
    internal class SilacInfoForDataGrid
    {
        public readonly SilacLabel SilacLabel;
        public readonly ChemicalFormula LabelFormula;
        public string Label { get; set; }
        public SilacInfoForDataGrid(SilacLabel label)
        {
            SilacLabel = label;
            LabelFormula = ChemicalFormula.ParseFormula(label.LabelChemicalFormula);
            Label = label.OriginalAminoAcid + "(" + label.MassDifference + ")";
        }
    }
}
