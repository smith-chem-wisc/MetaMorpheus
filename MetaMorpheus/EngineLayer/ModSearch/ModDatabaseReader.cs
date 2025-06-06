using Easy.Common.Extensions;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Security.Principal;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.ModSearch
{
    /// <summary>
    /// ModDatabaseReader is responsible for reading and processing modification data from a specified file.
    /// The loaded modifications include both glycans and regular modifications, which are then organized into a list of ModBox objects.
    /// </summary>
    public class ModDatabaseReader
    {
        /// Static variable to hold the modification dictionary
        public static List<Modification> ModDictionary { get; private set; } // A dictionary of all modifications in the MetaMorpheus based on "aListOfmods.txt"

        // Modification list
        public List<Modification> GlobalMods { get; private set; } // A list of all modifications, including glycans and regular modifications
        public List<Glycan> GlobalGlycans { get; private set; } // A list of glycans loaded from the database file
        public List<Modification> GlobalRegularMods { get; private set; } // A list of regular modifications excluding glycans

        public List<ModBox> ModBoxes { get; private set; } // A list of ModBox objects, each representing a combination of modifications
        public bool ToGenerateIons { get; } // Whether to generate ions for the modifications
        public SearchType SearchType { get; private set; } // The type of search being performed (e.g., O-glycan search, N-glycan search, regular mod search, mixed mod search)
        public int MaxModNum { get; } // The maximum number of modifications allowed on a peptide
        public int MaxGlycanNum { get; } // The maximum number of glycans allowed on a peptide

        public ModDatabaseReader(string filePath, int maxGlycanNum = 3, int maxModNum = 3, bool toGenerateIons = false, SearchType searchType = SearchType.RegularModSearch)
        {
            if (maxGlycanNum > maxModNum)
            {
                throw new ArgumentException("The maximum number of glycans cannot be greater than the maximum number of total modifications.");
            }

            MaxModNum = maxModNum;
            MaxGlycanNum = maxGlycanNum;
            ToGenerateIons = toGenerateIons;
            SearchType = searchType;
            Initized();
            ReadModDatabase(filePath, ToGenerateIons, searchType);
            MergeMod();
            BuildModBoxes();
        }

        /// <summary>
        /// Initializes the ModDictionary with modifications loaded from the specified file.
        /// </summary>
        public void Initized() 
        {
            ModDictionary =  new List<Modification>();
            ModDictionary = UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(Path.Combine(GlobalVariables.DataDir, @"Mods", @"Mods.txt"), out var errors).ToList();
            ModDictionary.Add(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(Path.Combine(GlobalVariables.DataDir, @"Mods", @"aListOfmods.txt"), out errors).ToList());
        }

        /// <summary>
        /// Read the modification database and store the modifications in GlobalGlycans and GlobalRegularMods.
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="ToGenerateIons"></param>
        /// <param name="searchType"></param>
        public void ReadModDatabase(string filePath, bool toGenerateIons, SearchType searchType) 
        {
            using (StreamReader lines = new StreamReader(filePath)) 
            {

                while (lines.Peek() != -1) 
                {
                    string line = lines.ReadLine();
                    if (line.Contains("HexNAc") && searchType == SearchType.O_GlycanSearch) // For now, we just search for O-glycan, then we will enable to search other type in the future
                    {
                        GlobalGlycans = GlycanDatabase.LoadKindGlycan(filePath, toGenerateIons, true).ToList();
                    }
                    else 
                    {
                        LoadRegularMod(line,lines);
                    }
                }
            }
        }

        /// <summary>
        /// Load regular modifications from the specified line and subsequent lines in the StreamReader.
        /// </summary>
        /// <param name="line"></param>
        /// <param name="lines"></param>
        public void LoadRegularMod(string line, StreamReader lines)
        {
            GlobalRegularMods = new List<Modification>();
            string name = line.Split('\t')[0];

            GlobalRegularMods.Add(ModDictionary.FirstOrDefault(p => p.IdWithMotif == name));
            while (lines.Peek() != -1)
            {
                line = lines.ReadLine();
                name = line.Split('\t')[0];
                GlobalRegularMods.Add(ModDictionary.FirstOrDefault(p => p.IdWithMotif == name));
            }

            GlobalRegularMods.RemoveAll(p => p == null);
        }

        /// <summary>
        /// Merge the glycans and regular modifications into a globalMods.
        /// </summary>
        public void MergeMod() 
        {
            GlobalMods = new List<Modification>();

            if (GlobalGlycans != null) 
            {
                foreach (var glycan in GlobalGlycans)
                {
                    Modification modificatinon = Glycan.OGlycanToModification(glycan);
                    GlobalMods.Add(modificatinon);
                }
            }

            if (GlobalRegularMods != null) 
            {
                foreach (var mod in GlobalRegularMods.Where(p => p != null))
                {
                    GlobalMods.Add(mod);
                }
            }
        }

        /// <summary>
        /// Builds ModBoxes from the global modifications, creating combinations of modifications based on the maximum number allowed.
        /// </summary>
        public void BuildModBoxes()
        {
            ModBoxes = new List<ModBox>();
            for (int i = 1; i < MaxModNum + 1; i++) 
            {
                foreach (var combination in ModBox.GetCombinationWithRept(Enumerable.Range(0, GlobalMods.Count()), i))
                {
                    var modForModBox = combination.Select(id => GlobalMods[id]).ToList();
                    if (modForModBox.Where(p => p.ModificationType == "O-Glycosylation").Count() > MaxGlycanNum) 
                    {
                        continue;
                    }
                    var modBox = new ModBox(modForModBox, combination.ToArray(), MaxModNum, false);
                    ModBoxes.Add(modBox);
                }
            }
            ModBoxes = ModBoxes.OrderBy(p => p.Mass).ToList();
        }
    }
}
