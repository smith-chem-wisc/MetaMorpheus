using Chemistry;
using System;
using System.Globalization;

namespace EngineLayer
{
    public class MetaMorpheusModification
    {

        #region Public Constructors

        public MetaMorpheusModification(string nameInXml, ModificationType type, char aminoAcid, string database, char prevAA, double precursorMassShift, double fragmentMassShift, double observedMassShift, ChemicalFormula cf)
        {
            NameInXml = nameInXml;
            ThisModificationType = type;
            AminoAcid = aminoAcid;
            Database = database;
            PrevAminoAcid = prevAA;
            PrecursorMassShift = precursorMassShift;
            FragmentMassShift = fragmentMassShift;
            ObservedMassShift = observedMassShift;
            ChemicalFormula = cf;
        }

        public MetaMorpheusModification(string NameInXml)
        {
            this.NameInXml = NameInXml;
        }

        public MetaMorpheusModification(double v)
        {
            this.FragmentMassShift = v;
            ThisModificationType = ModificationType.AminoAcidResidue;
            PrevAminoAcid = '\0';
            this.NameInXml = "";
        }

        #endregion Public Constructors

        #region Public Properties

        public double FragmentMassShift { get; private set; }

        public string Description
        {
            get
            {
                return Database + ":" + NameInXml + (Math.Abs(PrecursorMassShift - FragmentMassShift) > 1e-3 ? ":fms" + FragmentMassShift.ToString("F3", CultureInfo.InvariantCulture) : "");
            }
        }

        public ModificationType ThisModificationType { get; private set; }
        public char AminoAcid { get; private set; }
        public double PrecursorMassShift { get; private set; }
        public string Database { get; private set; }
        public string NameInXml { get; private set; }
        public char PrevAminoAcid { get; private set; }
        public ChemicalFormula ChemicalFormula { get; private set; }
        public double ObservedMassShift { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return Description;
        }

        #endregion Public Methods

    }
}