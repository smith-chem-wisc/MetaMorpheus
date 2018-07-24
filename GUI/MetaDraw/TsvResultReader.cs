using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using System.Text.RegularExpressions;
using Proteomics;
using MzLibUtil;
using System.Windows;
using Proteomics.ProteolyticDigestion;

namespace MetaMorpheusGUI
{
    public class TsvResultReader
    {
        private const string FULL_SEQUENCE_HEADER = "Full Sequence";
        private const string SCAN_NUMBER_HEADER = "Scan Number";
        private static char[] split = new char[] { '\t' };

        public static List<MetaDrawPsm> ReadTsv(string filePath)
        {
            List<MetaDrawPsm> psms = new List<MetaDrawPsm>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filePath);
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not read the file: " + e.Message);
                return psms;
            }

            int lineCount = 0;
            int fullSeqIndex = -1;
            int scanNumberIndex = -1;
            string line;
            while (reader.Peek() > 0)
            {
                lineCount++;

                line = reader.ReadLine();
                var spl = line.Split(split);

                if (lineCount == 1)
                {
                    fullSeqIndex = Array.IndexOf(spl, FULL_SEQUENCE_HEADER);
                    scanNumberIndex = Array.IndexOf(spl, SCAN_NUMBER_HEADER);
                }

                try
                {
                    int oneBasedScanNumber = int.Parse(spl[scanNumberIndex]);
                    string peptideSequence = spl[fullSeqIndex];
                    PeptideWithSetModifications peptide = ReadPeptideFromString(peptideSequence);

                    psms.Add(new MetaDrawPsm(oneBasedScanNumber, peptide));
                }
                catch (Exception)
                {
                    // TODO: write some kind of warning here?
                    // will skip ambiguous PSMs... and unreadable mods
                }
            }

            reader.Close();

            return psms;
        }

        private static PeptideWithSetModifications ReadPeptideFromString(string sequence)
        {
            // kind of weird... need to turn the peptide sequence into a protein and then digest it to get the peptide object...
            // TODO: make a FromString method in PeptideWithSetModifications to parse from string instead of doing this

            // read amino acid sequence and mods
            string aminoAcidSequence = "";
            int position = 0;

            Dictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>();
            int currentModPosition = 0;
            string currentMod = "";
            bool readingMod = false;
            
            for (int r = 0; r < sequence.Length; r++)
            {
                switch (sequence[r])
                {
                    case '[':
                        // new mod
                        readingMod = true;
                        currentModPosition = r;
                        break;
                    case ']':
                        // add the mod
                        var sp = currentMod.Split(new char[] { ':' });
                        string modType = sp[0];
                        string modId = sp[1];

                        var theMod = GlobalVariables.AllModsKnown.Where(v => v.id == modId && v.modificationType == modType).First() as ModificationWithMass;
                        if (oneBasedModifications.TryGetValue(currentModPosition, out var alreadyObservedMods))
                        {
                            alreadyObservedMods.Add(theMod);
                        }
                        else
                        {
                            oneBasedModifications.Add(currentModPosition, new List<Modification> { theMod });
                        }

                        readingMod = false;
                        currentMod = "";
                        break;
                    default:
                        // amino acid sequence or mod info
                        if (!readingMod)
                        {
                            aminoAcidSequence += sequence[r];
                            position++;
                        }
                        else
                        {
                            currentMod += sequence[r];
                        }
                        break;
                }
            }

            // create the peptide object w/ mods
            PeptideWithSetModifications peptide = new Protein(aminoAcidSequence, "", oneBasedModifications: oneBasedModifications)
                .Digest(new DigestionParams(protease: "top-down"), 
                new List<ModificationWithMass>(), new List<ModificationWithMass>()).Where(v => v.Sequence == sequence).First();

            return peptide;
        }
    }
}
