﻿using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;

namespace EngineLayer
{
    public class GlobalEngineLevelSettings
    {
        #region Public Fields

        public static readonly string elementsLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"elements.dat");

        #endregion Public Fields

        #region Private Fields

        private static readonly string unimodLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"unimod.xml");
        private static readonly string uniprotLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"ptmlist.txt");
        private static readonly string psiModLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"PSI-MOD.obo.xml");

        #endregion Private Fields

        #region Public Constructors

        static GlobalEngineLevelSettings()
        {
            try
            {
                UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
                UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation).ToList();
                PsiModDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(psiModLocation);
                var formalChargesDictionary = UsefulProteomicsDatabases.Loaders.GetFormalChargesDictionary(PsiModDeserialized);
                UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation, formalChargesDictionary).ToList();
            }
            catch (WebException)
            {
            }

            MetaMorpheusVersion = typeof(GlobalEngineLevelSettings).Assembly.GetName().Version.ToString();
            if (MetaMorpheusVersion.Equals("1.0.0.0"))
                MetaMorpheusVersion = "Not a release version";

            ProteaseDictionary = LoadProteaseDictionary();
            AllModsKnown = new List<Modification>();
        }

        #endregion Public Constructors

        #region Public Properties

        public static string MetaMorpheusVersion { get; }
        public static IEnumerable<Modification> UnimodDeserialized { get; }
        public static IEnumerable<Modification> UniprotDeseralized { get; }
        public static UsefulProteomicsDatabases.Generated.obo PsiModDeserialized { get; }

        public static Dictionary<string, Protease> ProteaseDictionary { get; }
        public static List<Modification> AllModsKnown { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static void AddMods(IEnumerable<Modification> enumerable)
        {
            foreach (var ye in enumerable)
            {
                if (string.IsNullOrEmpty(ye.modificationType) || string.IsNullOrEmpty(ye.id))
                    throw new MetaMorpheusException(ye.ToString() + Environment.NewLine + " has null or empty modification type");
                if (AllModsKnown.Any(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType) && !b.Equals(ye)))
                    throw new MetaMorpheusException("Modification id and type are equal, but some fields are not! Please modify/remove one of the modifications: " + Environment.NewLine + Environment.NewLine + ye.ToString() + Environment.NewLine + Environment.NewLine + " has same and id and modification type as " + Environment.NewLine + Environment.NewLine + AllModsKnown.First(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType)) + Environment.NewLine + Environment.NewLine);
                else if (AllModsKnown.Any(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType)))
                    continue;
                else
                    AllModsKnown.Add(ye);
            }
        }

        #endregion Public Methods

        #region Private Methods

        private static Dictionary<string, Protease> LoadProteaseDictionary()
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();
            using (StreamReader proteases = new StreamReader(Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", "proteases.tsv")))
            {
                proteases.ReadLine();

                while (proteases.Peek() != -1)
                {
                    string line = proteases.ReadLine();
                    string[] fields = line.Split('\t');

                    string name = fields[0];
                    string[] sequences_inducing_cleavage = fields[1].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    string[] sequences_preventing_cleavage = fields[2].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    var cleavage_terminus = (TerminusType)Enum.Parse(typeof(TerminusType), fields[3], true);
                    var cleavage_specificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true);
                    string psi_ms_accession_number = fields[5];
                    string psi_ms_name = fields[6];
                    string site_regexp = fields[7];
                    var protease = new Protease(name, sequences_inducing_cleavage, sequences_preventing_cleavage, cleavage_terminus, cleavage_specificity, psi_ms_accession_number, psi_ms_name, site_regexp);
                    dict.Add(protease.Name, protease);
                }
            }
            return dict;
        }

        #endregion Private Methods
    }
}