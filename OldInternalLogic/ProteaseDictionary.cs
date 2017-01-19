using System;
using System.Collections.Generic;
using System.IO;

namespace OldInternalLogic
{
    public class ProteaseDictionary : Dictionary<string, Protease>
    {

        #region Private Fields

        private static readonly ProteaseDictionary instance = new ProteaseDictionary();

        #endregion Private Fields

        #region Private Constructors

        private ProteaseDictionary()
        {
			using (StreamReader proteases = new StreamReader(Path.Combine(Environment.CurrentDirectory, "proteases.tsv")))
            {
                proteases.ReadLine();

                while (proteases.Peek() != -1)
                {
                    string line = proteases.ReadLine();
                    string[] fields = line.Split('\t');

                    string name = fields[0];
                    string[] sequences_inducing_cleavage = fields[1].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    string[] sequences_preventing_cleavage = fields[2].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    var cleavage_terminus = (OldLogicTerminus)Enum.Parse(typeof(OldLogicTerminus), fields[3], true);
                    var cleavage_specificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true);
                    string psi_ms_accession_number = fields[5];
                    string psi_ms_name = fields[6];
                    string site_regexp = fields[7];
                    var protease = new Protease(name, sequences_inducing_cleavage, sequences_preventing_cleavage, cleavage_terminus, cleavage_specificity, psi_ms_accession_number, psi_ms_name, site_regexp);
                    Add(protease);
                }
            }
        }

        #endregion Private Constructors

        #region Public Properties

        public static ProteaseDictionary Instance
        {
            get { return instance; }
        }

        #endregion Public Properties

        #region Internal Methods

        internal void Add(Protease protease)
        {
            Add(protease.Name, protease);
        }

        #endregion Internal Methods

    }
}