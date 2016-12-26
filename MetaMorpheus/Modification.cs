namespace MetaMorpheus
{
    public class MorpheusModification
    {
        public string Description
        {
            get
            {
                return Database + (labile ? ":labile" : "") + ":" + NameInXML;
            }
        }

        public bool labile { get; private set; }
        public ModificationType Type { get; private set; }
        public char AminoAcid { get; private set; }
        public double MonoisotopicMassShift { get; private set; }
        public bool DefaultFixed { get; private set; }
        public bool DefaultVariable { get; private set; }
        public string Database { get; private set; }
        public string DatabaseName { get; private set; }
        public string NameInXML { get; private set; }
        public char PrevAminoAcid { get; private set; }
        public double AlternativeMassShift { get; private set; }

        public MorpheusModification(string NameInXML, ModificationType type, char aminoAcid, double monoisotopicMassShift,
            string database, string databaseName, char prevAA, double AlternativeMassShift, bool labile)
        {
            this.NameInXML = NameInXML;
            Type = type;
            AminoAcid = aminoAcid;
            MonoisotopicMassShift = monoisotopicMassShift;
            Database = database;
            DatabaseName = databaseName;
            PrevAminoAcid = prevAA;
            this.AlternativeMassShift = AlternativeMassShift;
            this.labile = labile;
        }

        public MorpheusModification(string NameInXML)
        {
            this.NameInXML = NameInXML;
        }

        public override string ToString()
        {
            return Description;
        }
    }
}