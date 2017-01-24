using Chemistry;

namespace OldInternalLogic
{
    public class MorpheusModification
    {
        #region Public Constructors

        public MorpheusModification(string nameInXml, ModificationType type, char aminoAcid, string database, string databaseName, char prevAA, double alternativeMassShift, bool labile, ChemicalFormula cf)
        {
            this.NameInXml = nameInXml;
            ThisModificationType = type;
            AminoAcid = aminoAcid;
            MonoisotopicMassShift = cf.MonoisotopicMass;
            Database = database;
            DatabaseName = databaseName;
            PrevAminoAcid = prevAA;
            this.AlternativeMassShift = alternativeMassShift;
            this.Labile = labile;
            this.ChemicalFormula = cf;
        }

        public MorpheusModification(string NameInXml)
        {
            this.NameInXml = NameInXml;
        }

        public MorpheusModification(double v)
        {
            this.MonoisotopicMassShift = v;
            ThisModificationType = ModificationType.AminoAcidResidue;
            PrevAminoAcid = '\0';
            this.NameInXml = "";
        }

        #endregion Public Constructors

        #region Public Properties

        public string Description
        {
            get
            {
                return Database + (Labile ? ":labile" : "") + ":" + NameInXml;
            }
        }

        public bool Labile { get; private set; }
        public ModificationType ThisModificationType { get; private set; }
        public char AminoAcid { get; private set; }
        public double MonoisotopicMassShift { get; private set; }
        public string Database { get; private set; }
        public string DatabaseName { get; private set; }
        public string NameInXml { get; private set; }
        public char PrevAminoAcid { get; private set; }
        public double AlternativeMassShift { get; private set; }
        public ChemicalFormula ChemicalFormula { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return Description;
        }

        #endregion Public Methods
    }
}
