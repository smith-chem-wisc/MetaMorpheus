using FlashLFQ;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Easy.Common.Extensions;
using UsefulProteomicsDatabases;
using IO.Mgf;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using System.Text.RegularExpressions;
using Proteomics;

namespace EngineLayer.PsmTsv
{
    internal enum PsmFileType
    { MetaMorpheus, Morpheus, MaxQuant, MaxQuantEvidence, PeptideShaker, Generic, Percolator, Unknown }

    public class PsmGenericReader
    {
        internal static int _fileNameCol;
        internal static int _baseSequCol;
        internal static int _fullSequCol;
        internal static int _monoMassCol;
        internal static int _msmsRetnCol;
        internal static int _msmsScanCol;
        internal static int _chargeStCol;
        internal static int _protNameCol;
        internal static int _decoyCol;
        internal static int _qValueCol;
        internal static int _qValueNotchCol;

        // optional columns
        internal static int _geneNameCol;
        internal static int _organismCol;

        // Match between runs columns
        internal static int _matchScoreCol;
        internal static int _matchMzDeltaCol;
        internal static int _matchRtDeltaCol;
        internal static int _retentionLengthCol;
        internal static int _intensityCol;
        internal static int _mzCol;



        internal static Dictionary<string, double> _modSequenceToMonoMass;
        internal static Dictionary<string, FlashLFQ.ProteinGroup> allProteinGroups;
        internal static List<ScanHeaderInfo> _scanHeaderInfo = new List<ScanHeaderInfo>();

        //Delimiters refere to contents of one field, not the delimiter between fields
        internal static readonly Dictionary<PsmFileType, string[]> delimiters = new Dictionary<PsmFileType, string[]>
        {
            { PsmFileType.MetaMorpheus, new string[] { "|", " or " } },
            { PsmFileType.Morpheus, new string[] { ";" } },
            { PsmFileType.MaxQuant, new string[] { ";" } },
            { PsmFileType.MaxQuantEvidence, new string[] { ";" } },
            { PsmFileType.Percolator, new string[] { "," } },
            { PsmFileType.Generic, new string[] { ";" } },
            { PsmFileType.PeptideShaker, new string[] { ", " } },
        };

        /// <summary>
        /// Reads results files from various search software and creates a list of identifications
        /// Metamorpheus, Morpheus, MaxQuant, PeptideShaker, and Percolator are all supported
        /// </summary>
        /// <param name="filepath">Path to the results file</param>
        /// <param name="silent"></param>
        /// <param name="rawfiles">List of data files (as SpectraFileInfo) that are present in the results</param>
        /// <param name="mbrPeaks">Optional param, a list of chromatographic peaks that will be populated with
        /// any Match Betweens Runs peaks in a MaxQuant evidence file</param>
        /// <returns>A list of Identifications</returns>
        /// <exception cref="Exception"></exception>
        public static List<Identification> ReadPsms(string filepath, bool silent, List<SpectraFileInfo> rawfiles, List<ChromatographicPeak> mbrPeaks = null)
        {
            if (_modSequenceToMonoMass == null)
            {
                _modSequenceToMonoMass = new Dictionary<string, double>();
            }

            if (allProteinGroups == null)
            {
                allProteinGroups = new Dictionary<string, FlashLFQ.ProteinGroup>();
            }

            var rawFileDictionary = rawfiles.ToDictionary(p => p.FilenameWithoutExtension, v => v);
            List<Identification> flashLfqIdentifications = new();
            PsmFileType fileType = PsmFileType.Unknown;

            if (!silent)
            {
                Console.WriteLine("Opening PSM file " + filepath);
            }

            StreamReader reader;
            try
            {
                reader = new StreamReader(filepath);
            }
            catch (Exception e)
            {
                if (!silent)
                {
                    Console.WriteLine("Error reading file " + filepath + "\n" + e.Message);
                }

                return new List<Identification>();
            }

            List<string> inputPsms = File.ReadAllLines(filepath).ToList();
            string[] header = inputPsms[0].Split('\t');

            try
            {
                fileType = GetFileTypeFromHeader(inputPsms[0]);
                inputPsms.RemoveAt(0);
            }
            catch
            {
                throw new Exception("Could not interpret PSM header labels from file: " + filepath);
            }

            var psmsGroupedByFile = inputPsms.GroupBy(p => MzLibUtil.PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(p.Split('\t')[_fileNameCol])).ToList();

            foreach (var fileSpecificPsms in psmsGroupedByFile)
            {
                int myFileIndex = rawFileDictionary.Keys.ToList().IndexOf(fileSpecificPsms.Key);
                string fullFilePathWithExtension = rawFileDictionary[fileSpecificPsms.Key].FullFilePathWithExtension;
                List<Identification> myFileIndentifications = new();
                List<ScanHeaderInfo> scanHeaderInfo = new();
                switch (fileType)
                {
                    case (PsmFileType.Percolator):
                        scanHeaderInfo = scanHeaderInfo = ScanInfoRecovery.FileScanHeaderInfo(fullFilePathWithExtension);
                        foreach (var psm in fileSpecificPsms)
                        {
                            try
                            {
                                Identification id = GetPercolatorIdentification(psm, scanHeaderInfo, silent, rawFileDictionary);
                                if (id != null)
                                {
                                    myFileIndentifications.Add(id);
                                }
                            }
                            catch (Exception e)
                            {
                                if (!silent)
                                {
                                    Console.WriteLine("Problem reading line in the identification file" + "; " + e.Message);
                                }
                            }
                        }

                        break;

                    case (PsmFileType.MaxQuantEvidence):
                        foreach (var psm in fileSpecificPsms)
                        {
                            try
                            {
                                Identification id = GetIdentification(psm, silent, rawFileDictionary, fileType);
                                PeptideWithSetModifications pswm = GetPWSM(psm, silent, rawFileDictionary, fileType);
                                if (id != null)
                                {
                                    if (psm.Split('\t')[_matchScoreCol].IsNullOrEmptyOrWhiteSpace())
                                    {
                                        myFileIndentifications.Add(id);
                                    }
                                    else if (mbrPeaks != null)
                                    {
                                        mbrPeaks.Add(GetMbrPeak(id, psm, rawFileDictionary));
                                    }

                                }
                            }
                            catch (Exception e)
                            {
                                if (!silent)
                                {
                                    Console.WriteLine("Problem reading line in the identification file" + "; " + e.Message);
                                }
                            }
                        }

                        break;

                    default:
                        foreach (var psm in fileSpecificPsms)
                        {
                            try
                            {
                                Identification id = GetIdentification(psm, silent, rawFileDictionary, fileType);
                                if (id != null)
                                {
                                    myFileIndentifications.Add(id);
                                }
                            }
                            catch (Exception e)
                            {
                                if (!silent)
                                {
                                    Console.WriteLine("Problem reading line in the identification file" + "; " + e.Message);
                                }
                            }
                        }

                        break;
                }

                _scanHeaderInfo.AddRange(scanHeaderInfo);
                flashLfqIdentifications.AddRange(myFileIndentifications);
            }

            if (!silent)
            {
                Console.WriteLine("Done reading PSMs; found " + flashLfqIdentifications.Count);
            }

            reader.Close();
            return flashLfqIdentifications;
        }

        internal static ChromatographicPeak GetMbrPeak(Identification id, string psm,
            Dictionary<string, SpectraFileInfo> fileDictionary)
        {
            string[] psmSplit = psm.Split('\t');
            string fileName = psmSplit[_fileNameCol];
            ChromatographicPeak mbrPeak = new ChromatographicPeak(id, true, fileDictionary[fileName]);
            IndexedMassSpectralPeak msPeak = new IndexedMassSpectralPeak(
                Double.Parse(psmSplit[_mzCol]),
                Double.Parse(psmSplit[_intensityCol]),
                zeroBasedMs1ScanIndex: -1,
                Double.Parse(psmSplit[_msmsRetnCol]));
            FlashLFQ.IsotopicEnvelope envelope =
                new FlashLFQ.IsotopicEnvelope(msPeak, Int32.Parse(psmSplit[_chargeStCol]), msPeak.Intensity);
            mbrPeak.IsotopicEnvelopes.Add(envelope);
            return mbrPeak;
        }

        internal static PeptideWithSetModifications GetPWSM(string psmString, bool silent,
            Dictionary<string, SpectraFileInfo> fileDictionary, PsmFileType fileType)
        {


            return null;
        }

        public static string MaxQuantToPWSM(string mqFullSeq, List<Modification> fixedMods = null)
        {
            string mqTrimmedSeq = mqFullSeq.Trim('_');
            string fullSeq = null;

            if (FindVariableModsInMaxQuantSequence(mqTrimmedSeq, out var modDict))
            {
                // Matches all AA sequences except for the last
                string internalAAPattern = @"([A-Z]*)\(";
                Regex internalAARegex = new Regex(internalAAPattern);
                IEnumerable<string> internalAAs = internalAARegex.Matches(mqTrimmedSeq)
                    .Select(m => m.Groups[1].Value)
                    .Where(s => s.IsNotNullOrEmptyOrWhiteSpace());

                // The second  regular expression ( "\([A-Z]*" ) captures the last section
                // of the peptide. 
                string finalAAPattern = @"([A-Z]*)$";
                Regex finalAARegex = new Regex(finalAAPattern);
                var match = finalAARegex.Match(mqTrimmedSeq);
                string finalAAs = finalAARegex.Match(mqTrimmedSeq).Groups[1].Value;

                    var orderedMods = modDict.
                    OrderBy(kvp => kvp.Key).ToList();

                if (orderedMods[0].Key == 1)
                {
                    fullSeq = string.Concat(
                        orderedMods.
                            Select(kvp => kvp.Value).
                            Select(mod => '[' + mod.ModificationType + ":" + mod.IdWithMotif + ']').
                            Zip(internalAAs.Append(finalAAs), (first, second) => first + second)
                    );
                }
                else
                {
                    fullSeq = string.Concat(
                        internalAAs.
                            Zip(orderedMods, (first, second) => 
                                first + '[' + second.Value.ModificationType + ":" + second.Value.IdWithMotif + ']')
                            // Zip combines elements into pairs, which results in a dangling tail,
                            // which requires manually adding the last peptide sequence.
                            .Append(finalAAs)
                    );
                }
            }
            else
            {
                fullSeq = mqTrimmedSeq;
            }

            return AddFixedMods(fullSeq, fixedMods, out var modPositions);
        }

        /// <summary>
        /// Adds fixed modifications to a string containing a full peptide sequence. 
        /// </summary>
        /// <param name="fullSequence">The sequence to be modified</param>
        /// <param name="fixedMods">A list of fixed modifications.
        /// If null, defaults to Carbamidomethylation on C</param>
        /// <returns>String containing the full sequence in MetaMorpheus format</returns>
        public static string AddFixedMods(string fullSequence, List<Modification> fixedMods, 
            out Dictionary<int, Modification> modPositions)
        {
            modPositions = new();
            string fullSeqWithFixedMods = fullSequence;
            fixedMods ??= new List<Modification> { GlobalVariables.AllModsKnownDictionary["Carbamidomethyl on C"] };
            foreach (Modification mod in fixedMods)
            {
                string location = mod.Target.ToString();
                if (location.Length == 1) // Single AA target
                {
                    int modPosition = 1;

                    // Pattern employs a negative lookbehind - "(?<!\[|\:)" - as not to match existing modifications
                    // e.g. (PEPT[Common Biological:Phosphorylation on T]IDE), the C in Common will not be matched
                    string fixedModPattern = @"(?<!\[|\:)(" + location + ")";
                    Regex fixedModRegex = new(fixedModPattern);
                    string[] sequenceFragments = fixedModRegex.Split(fullSequence);

                    if (sequenceFragments.Length > 1)
                    {
                        StringBuilder sb = new();

                        foreach (string seqFragment in sequenceFragments)
                        {
                            if (seqFragment.Equals(location))
                            {
                                modPositions.Add(modPosition, mod);
                                sb.Append(seqFragment + '[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
                            }
                            else
                            {
                                sb.Append(seqFragment);
                            }

                            modPosition += seqFragment.Length;
                        }

                        fullSeqWithFixedMods = sb.ToString();
                    }
                }
                // Fixed mods with variable positions are currently unhandled
            }
            return fullSeqWithFixedMods;
        }


        /// <summary>
        /// Converts a MaxQuant sequence to a MetaMorpheus sequence
        /// </summary>
        /// <param name="fullSequence"> Full sequence of the peptide in question, taken from a MaxQuant file</param>
        /// <param name="modDict">A dictionary with integer positions as keys, Modifications as values. 1 represents the n-terminus</param>
        /// <returns> True if modifications were present, false otherwise </returns>
        public static bool FindVariableModsInMaxQuantSequence(string fullSeq, out Dictionary<int, Modification> modDict)
        {
            // use a regex to get all modifications
            string pattern = @"(.\([A-Z][a-z]*[\s]\([^\(]*\)\).)"; // Matches single modification + leading and trailing AAs
            Regex regex = new(pattern);
            MatchCollection matches = regex.Matches(fullSeq);

            // Here, int is a one based index (one represents the n-terminus) of the location of each mod
            // within the peptide
            if (matches.Any())
            {
                int currentPosition = 0;
                modDict = new();
                foreach (Match match in matches)
                {
                    GroupCollection group = match.Groups;
                    string modString = group[1].Value;
                    int startIndex = group[0].Index;
                    int captureLength = group[0].Length;

                    // int to add is startIndex - current position + 1 (one based indexing)
                    int positionToAddToDict = 1 + startIndex - currentPosition;
                    modDict.Add(positionToAddToDict, ParseMaxQuantMod(modString));

                    currentPosition += startIndex + captureLength;
                }
                return true;
            }

            modDict = null;
            return false;

        }

        /// <summary>
        /// Returns a Modification object that most closely matches the modification
        /// found in a MaxQuant results file
        /// </summary>
        /// <param name="modString">A MaxQuant modification, in the form of (Mod Name (Location))</param>
        /// <returns>A Modification object that most closely matches the input string</returns>
        public static Modification ParseMaxQuantMod(string modString)
        {
            if (MaxQuantToMetaModDictionary == null) MaxQuantToMetaModDictionary = new();

            if (MaxQuantToMetaModDictionary.TryGetValue(modString, out var metaMod))
            {
                return metaMod;
            }

            string modificationPattern = @"\((.*)\(";
            Regex modRegex = new(modificationPattern);
            string mod = modRegex.Match(modString).Groups[1].Value;

            string location = TranslateMaxQuantPosition(modString);

            // Search MM modifications to find the closest match
            metaMod = GlobalVariables.AllModsKnown.
                OrderBy(m => ComputeLevenshteinDistance(mod, m.OriginalId)).
                Where(m => m.Target.ToString().Equals(location)).
                FirstOrDefault();


            if (metaMod != null) MaxQuantToMetaModDictionary.Add(modString, metaMod);
            return metaMod;

        }

        public static Dictionary<string, Modification> MaxQuantToMetaModDictionary { get; set; }

        // Ripped straight from StackOverflow:
        // https://stackoverflow.com/questions/13793560/find-closest-match-to-input-string-in-a-list-of-strings
        /// <summary>
        /// Compute the distance between two strings.
        /// </summary>
        /// <returns> An int representing the string distance, where a value of 0 represents a perfect match</returns>
        public static int ComputeLevenshteinDistance(string s, string t)
        {
            int n = s.Length;
            int m = t.Length;
            int[,] d = new int[n + 1, m + 1];

            // Step 1
            if (n == 0)
            {
                return m;
            }

            if (m == 0)
            {
                return n;
            }

            // Step 2
            for (int i = 0; i <= n; d[i, 0] = i++)
            {
            }

            for (int j = 0; j <= m; d[0, j] = j++)
            {
            }

            // Step 3
            for (int i = 1; i <= n; i++)
            {
                //Step 4
                for (int j = 1; j <= m; j++)
                {
                    // Step 5
                    int cost = (t[j - 1] == s[i - 1]) ? 0 : 1;

                    // Step 6
                    d[i, j] = Math.Min(
                        Math.Min(d[i - 1, j] + 1, d[i, j - 1] + 1),
                        d[i - 1, j - 1] + cost);
                }
            }
            // Step 7
            return d[n, m];
        }

        /// <summary>
        /// Translates a position given by MaxQuant to the Metamorpheus format
        /// </summary>
        /// <param name="location"> MaxQuant position as a trimmed string</param>
        /// <returns> Position in the MetaMorpheus format</returns>
        public static string TranslateMaxQuantPosition(string modString)
        {
            string locationPattern = @"\(.*\(([^\)]*)\)";
            Regex locationRegex = new(locationPattern);
            string location = locationRegex.Match(modString).Groups[1].Value;

            // Assume that any 1 character location is a single amino acid
            if (location.Length == 1) return location; 
            // Any hyphenated location is a terminus
            else switch (location) 
            {
                case "Protein N-term":
                case "N-term":
                    return modString[modString.Length - 1].ToString();

                case "Protein C-term":
                case "C-term": 
                    return modString[0].ToString();

                default:
                    return "X";
            }


        }

        internal static Identification GetIdentification(string line, bool silent, Dictionary<string, SpectraFileInfo> rawFileDictionary, PsmFileType fileType)
        {
            var param = line.Split('\t');

            // only quantify PSMs below 1% FDR with MetaMorpheus/Morpheus results
            if (fileType == PsmFileType.MetaMorpheus && double.Parse(param[_qValueCol], CultureInfo.InvariantCulture) > 0.01)
            {
                return null;
            }
            else if (fileType == PsmFileType.Morpheus && double.Parse(param[_qValueCol], CultureInfo.InvariantCulture) > 1.00)
            {
                return null;
            }

            // only quantify PSMs below 1% notch FDR with MetaMorpheus/Morpheus results
            if (fileType == PsmFileType.MetaMorpheus && double.Parse(param[_qValueNotchCol], CultureInfo.InvariantCulture) > 0.01)
            {
                return null;
            }

            // skip decoys with MetaMorpheus/Morpheus results
            //TODO: what about decoys from other input types?
            if ((fileType == PsmFileType.MetaMorpheus || fileType == PsmFileType.Morpheus) &&
                param[_decoyCol].Contains("D"))
            {
                return null;
            }

            // spectrum file name
            string fileName = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(param[_fileNameCol]);

            // base sequence
            string baseSequence = param[_baseSequCol];

            // modified sequence
            string modSequence = param[_fullSequCol];

            // skip ambiguous sequence in MetaMorpheus output
            if (fileType == PsmFileType.MetaMorpheus && (modSequence.Contains(" or ") || modSequence.Contains("|") || modSequence.ToLowerInvariant().Contains("too long")))
            {
                return null;
            }

            // monoisotopic mass
            if (double.TryParse(param[_monoMassCol], NumberStyles.Number, CultureInfo.InvariantCulture, out double monoisotopicMass))
            {
                if (_modSequenceToMonoMass.TryGetValue(modSequence, out double storedMonoisotopicMass))
                {
                    if (storedMonoisotopicMass != monoisotopicMass)
                    {
                        if (!silent)
                        {
                            Console.WriteLine("Caution! PSM with could not be read. A peptide with the same modified sequence but a different monoisotopic mass has already been added." + "\n"
                                + line);
                        }

                        return null;
                    }
                }
                else
                {
                    _modSequenceToMonoMass.Add(modSequence, monoisotopicMass);
                }
            }
            else
            {
                if (!silent)
                {
                    Console.WriteLine("Caution! PSM with could not be read. Monoisotopic mass not interpretable." + "\n"
                                + line);
                }

                return null;
            }

            // retention time
            double ms2RetentionTime = -1;
            if (double.TryParse(param[_msmsRetnCol], NumberStyles.Number, CultureInfo.InvariantCulture, out double retentionTime))
            {
                ms2RetentionTime = retentionTime;
                if (fileType == PsmFileType.PeptideShaker)
                {
                    // peptide shaker RT is in seconds - convert to minutes
                    ms2RetentionTime = retentionTime / 60.0;
                }

                if (ms2RetentionTime < 0)
                {
                    if (!silent)
                    {
                        Console.WriteLine("Caution! PSM retention time was negative." + "\n"
                            + line);
                    }

                    return null;
                }
            }
            else
            {
                if (!silent)
                {
                    Console.WriteLine("PSM retention time was not interpretable." + "\n"
                        + line);
                }

                return null;
            }

            // charge state
            int chargeState;
            if (fileType == PsmFileType.PeptideShaker)
            {
                string chargeStringNumbersOnly = new String(param[_chargeStCol].Where(Char.IsDigit).ToArray());

                if (string.IsNullOrWhiteSpace(chargeStringNumbersOnly))
                {
                    if (!silent)
                    {
                        Console.WriteLine("PSM charge state was not interpretable." + "\n"
                            + line);
                    }

                    return null;
                }
                else
                {
                    if (!int.TryParse(chargeStringNumbersOnly, out chargeState))
                    {
                        if (!silent)
                        {
                            Console.WriteLine("PSM charge state was not interpretable." + "\n"
                            + line);
                        }

                        return null;
                    }
                }
            }
            else
            {
                if (!double.TryParse(param[_chargeStCol], NumberStyles.Number, CultureInfo.InvariantCulture, out double chargeStateDouble))
                {
                    if (!silent)
                    {
                        Console.WriteLine("PSM charge state was not interpretable." + "\n"
                            + line);
                    }

                    return null;
                }

                chargeState = (int)chargeStateDouble;
            }

            // protein groups
            // use all proteins listed
            string[] proteins = null;
            string[] genes = null;
            string[] organisms = null;
            List<FlashLFQ.ProteinGroup> proteinGroups = new List<FlashLFQ.ProteinGroup>();
            proteins = param[_protNameCol].Split(delimiters[fileType], StringSplitOptions.None);

            if (_geneNameCol >= 0)
            {
                genes = param[_geneNameCol].Split(delimiters[fileType], StringSplitOptions.None);
            }

            if (_organismCol >= 0)
            {
                organisms = param[_organismCol].Split(delimiters[fileType], StringSplitOptions.None);
            }

            for (int pr = 0; pr < proteins.Length; pr++)
            {
                string proteinName = proteins[pr];
                string gene = "";
                string organism = "";

                if (genes != null)
                {
                    if (genes.Length == 1)
                    {
                        gene = genes[0];
                    }
                    else if (genes.Length == proteins.Length)
                    {
                        gene = genes[pr];
                    }
                    else if (proteins.Length == 1)
                    {
                        gene = param[_geneNameCol];
                    }
                }

                if (organisms != null)
                {
                    if (organisms.Length == 1)
                    {
                        organism = organisms[0];
                    }
                    else if (organisms.Length == proteins.Length)
                    {
                        organism = organisms[pr];
                    }
                    else if (proteins.Length == 1)
                    {
                        organism = param[_organismCol];
                    }
                }

                if (allProteinGroups.TryGetValue(proteinName, out FlashLFQ.ProteinGroup pg))
                {
                    proteinGroups.Add(pg);
                }
                else
                {
                    FlashLFQ.ProteinGroup newPg = new FlashLFQ.ProteinGroup(proteinName, gene, organism);
                    allProteinGroups.Add(proteinName, newPg);
                    proteinGroups.Add(newPg);
                }
            }

            if (!rawFileDictionary.TryGetValue(fileName, out SpectraFileInfo spectraFileInfoToUse))
            {
                // skip PSMs for files with no spectrum data input
                return null;
            }

            // construct id
            return new Identification(spectraFileInfoToUse, baseSequence, modSequence, monoisotopicMass, ms2RetentionTime, chargeState, proteinGroups);
        }

        internal static Identification GetPercolatorIdentification(string line, List<ScanHeaderInfo> scanHeaderInfo, bool silent, Dictionary<string, SpectraFileInfo> rawFileDictionary)
        {
            var param = line.Split('\t');

            // spectrum file name
            string fileName = param[_fileNameCol];

            // base sequence
            string baseSequence = null;

            // modified sequence
            string modSequence = param[_fullSequCol];

            // skip ambiguous sequence in MetaMorpheus output
            if (modSequence.Contains("|") || modSequence.ToLowerInvariant().Contains("too long"))
            {
                return null;
            }

            // monoisotopic mass
            if (double.TryParse(param[_monoMassCol], NumberStyles.Number, CultureInfo.InvariantCulture, out double monoisotopicMass))
            {
                if (_modSequenceToMonoMass.TryGetValue(modSequence, out double storedMonoisotopicMass))
                {
                    if (storedMonoisotopicMass != monoisotopicMass)
                    {
                        if (!silent)
                        {
                            Console.WriteLine("Caution! PSM with could not be read. A peptide with the same modified sequence but a different monoisotopic mass has already been added." + "\n"
                                + line);
                        }
                        return null;
                    }
                }
                else
                {
                    _modSequenceToMonoMass.Add(modSequence, monoisotopicMass);
                }
            }
            else
            {
                if (!silent)
                {
                    Console.WriteLine("Caution! PSM with could not be read. Monoisotopic mass not interpretable." + "\n"
                                + line);
                }
                return null;
            }

            // retention time
            double ms2RetentionTime = -1;
            //percolator input files do not have retention times. So, we have to get them from the data file using the scan number.

            if (int.TryParse(param[_msmsScanCol], NumberStyles.Number, CultureInfo.InvariantCulture, out int scanNumber))
            {
                ms2RetentionTime = scanHeaderInfo.Where(i => PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(i.FileNameWithoutExtension) == PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(fileName) && i.ScanNumber == scanNumber).FirstOrDefault().RetentionTime;
            }

            // charge state
            int chargeState;

            if (!double.TryParse(param[_chargeStCol], NumberStyles.Number, CultureInfo.InvariantCulture, out double chargeStateDouble))
            {
                if (!silent)
                {
                    Console.WriteLine("Caution! PSM with could not be read. Charge state not interpretable." + "\n"
                            + line);
                }

                return null;
            }

            chargeState = (int)chargeStateDouble;

            // protein groups
            // use all proteins listed
            string[] proteins = null;
            string[] genes = null;
            string[] organisms = null;
            List<FlashLFQ.ProteinGroup> proteinGroups = new List<FlashLFQ.ProteinGroup>();
            proteins = param[_protNameCol].Split(delimiters[PsmFileType.Percolator], StringSplitOptions.None);

            if (_geneNameCol >= 0)
            {
                genes = param[_geneNameCol].Split(delimiters[PsmFileType.Percolator], StringSplitOptions.None);
            }

            if (_organismCol >= 0)
            {
                organisms = param[_organismCol].Split(delimiters[PsmFileType.Percolator], StringSplitOptions.None);
            }

            for (int pr = 0; pr < proteins.Length; pr++)
            {
                string proteinName = proteins[pr];
                string gene = "";
                string organism = "";

                if (genes != null)
                {
                    if (genes.Length == 1)
                    {
                        gene = genes[0];
                    }
                    else if (genes.Length == proteins.Length)
                    {
                        gene = genes[pr];
                    }
                    else if (proteins.Length == 1)
                    {
                        gene = param[_geneNameCol];
                    }
                }

                if (organisms != null)
                {
                    if (organisms.Length == 1)
                    {
                        organism = organisms[0];
                    }
                    else if (organisms.Length == proteins.Length)
                    {
                        organism = organisms[pr];
                    }
                    else if (proteins.Length == 1)
                    {
                        organism = param[_organismCol];
                    }
                }

                if (allProteinGroups.TryGetValue(proteinName, out FlashLFQ.ProteinGroup pg))
                {
                    proteinGroups.Add(pg);
                }
                else
                {
                    FlashLFQ.ProteinGroup newPg = new FlashLFQ.ProteinGroup(proteinName, gene, organism);
                    allProteinGroups.Add(proteinName, newPg);
                    proteinGroups.Add(newPg);
                }
            }

            if (!rawFileDictionary.TryGetValue(PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(fileName), out SpectraFileInfo spectraFileInfoToUse))
            {
                // skip PSMs for files with no spectrum data input
                return null;
            }

            return new Identification(spectraFileInfoToUse, baseSequence, modSequence, monoisotopicMass, ms2RetentionTime, chargeState, proteinGroups);
        }

        internal static PsmFileType GetFileTypeFromHeader(string header)
        {
            PsmFileType type = PsmFileType.Unknown;

            var split = header.Split('\t').Select(p => p.ToLowerInvariant()).ToArray();

            // MetaMorpheus MS/MS input
            if (split.Contains("File Name".ToLowerInvariant())
                        && split.Contains("Base Sequence".ToLowerInvariant())
                        && split.Contains("Full Sequence".ToLowerInvariant())
                        && split.Contains("Peptide Monoisotopic Mass".ToLowerInvariant())
                        && split.Contains("Scan Retention Time".ToLowerInvariant())
                        && split.Contains("Precursor Charge".ToLowerInvariant())
                        && split.Contains("Protein Accession".ToLowerInvariant())
                        && split.Contains("Decoy/Contaminant/Target".ToLowerInvariant())
                        && split.Contains("QValue".ToLowerInvariant())
                        && split.Contains("QValue Notch".ToLowerInvariant()))
            {
                _fileNameCol = Array.IndexOf(split, "File Name".ToLowerInvariant());
                _baseSequCol = Array.IndexOf(split, "Base Sequence".ToLowerInvariant());
                _fullSequCol = Array.IndexOf(split, "Full Sequence".ToLowerInvariant());
                _monoMassCol = Array.IndexOf(split, "Peptide Monoisotopic Mass".ToLowerInvariant());
                _msmsRetnCol = Array.IndexOf(split, "Scan Retention Time".ToLowerInvariant());
                _chargeStCol = Array.IndexOf(split, "Precursor Charge".ToLowerInvariant());
                _protNameCol = Array.IndexOf(split, "Protein Accession".ToLowerInvariant());
                _decoyCol = Array.IndexOf(split, "Decoy/Contaminant/Target".ToLowerInvariant());
                _qValueCol = Array.IndexOf(split, "QValue".ToLowerInvariant());
                _qValueNotchCol = Array.IndexOf(split, "QValue Notch".ToLowerInvariant());
                _geneNameCol = Array.IndexOf(split, "Gene Name".ToLowerInvariant());
                _organismCol = Array.IndexOf(split, "Organism Name".ToLowerInvariant());

                return PsmFileType.MetaMorpheus;
            }

            // Morpheus MS/MS input
            else if (split.Contains("Filename".ToLowerInvariant())
                && split.Contains("Base Peptide Sequence".ToLowerInvariant())
                && split.Contains("Peptide Sequence".ToLowerInvariant())
                && split.Contains("Theoretical Mass (Da)".ToLowerInvariant())
                && split.Contains("Retention Time (minutes)".ToLowerInvariant())
                && split.Contains("Precursor Charge".ToLowerInvariant())
                && split.Contains("Protein Description".ToLowerInvariant())
                && split.Contains("Decoy?".ToLowerInvariant())
                && split.Contains("Q-Value (%)".ToLowerInvariant()))
            {
                _fileNameCol = Array.IndexOf(split, "Filename".ToLowerInvariant());
                _baseSequCol = Array.IndexOf(split, "Base Peptide Sequence".ToLowerInvariant());
                _fullSequCol = Array.IndexOf(split, "Peptide Sequence".ToLowerInvariant());
                _monoMassCol = Array.IndexOf(split, "Theoretical Mass (Da)".ToLowerInvariant());
                _msmsRetnCol = Array.IndexOf(split, "Retention Time (minutes)".ToLowerInvariant());
                _chargeStCol = Array.IndexOf(split, "Precursor Charge".ToLowerInvariant());
                _protNameCol = Array.IndexOf(split, "Protein Description".ToLowerInvariant());
                _decoyCol = Array.IndexOf(split, "Decoy?".ToLowerInvariant());
                _qValueCol = Array.IndexOf(split, "Q-Value (%)".ToLowerInvariant());

                _geneNameCol = Array.IndexOf(split, "Gene Name".ToLowerInvariant()); // probably doesn't exist
                _organismCol = Array.IndexOf(split, "Organism Name".ToLowerInvariant());

                return PsmFileType.Morpheus;
            }

            // MaxQuant Evidence input (contains Match between runs information)
            else if (split.Contains("Raw file".ToLowerInvariant())
                     && split.Contains("Sequence".ToLowerInvariant())
                     && split.Contains("Modified sequence".ToLowerInvariant())
                     && split.Contains("Mass".ToLowerInvariant())
                     && split.Contains("Retention time".ToLowerInvariant())
                     && split.Contains("Charge".ToLowerInvariant())
                     && split.Contains("Proteins".ToLowerInvariant())
                     && split.Contains("Match score".ToLowerInvariant()))
            {
                _fileNameCol = Array.IndexOf(split, "Raw file".ToLowerInvariant());
                _baseSequCol = Array.IndexOf(split, "Sequence".ToLowerInvariant());
                _fullSequCol = Array.IndexOf(split, "Modified sequence".ToLowerInvariant());
                _monoMassCol = Array.IndexOf(split, "Mass".ToLowerInvariant());
                _msmsRetnCol = Array.IndexOf(split, "Retention time".ToLowerInvariant());
                _chargeStCol = Array.IndexOf(split, "Charge".ToLowerInvariant());
                _protNameCol = Array.IndexOf(split, "Proteins".ToLowerInvariant());

                // optional
                _geneNameCol = Array.IndexOf(split, "Gene Names".ToLowerInvariant());
                _organismCol = Array.IndexOf(split, "Organism Name".ToLowerInvariant());

                // MBR
                _matchScoreCol = Array.IndexOf(split, "Match score".ToLowerInvariant());
                _matchMzDeltaCol = Array.IndexOf(split, "Match m/z difference".ToLowerInvariant());
                _matchRtDeltaCol = Array.IndexOf(split, "Match time difference".ToLowerInvariant());
                _retentionLengthCol = Array.IndexOf(split, "Retention length".ToLowerInvariant());
                _msmsScanCol = Array.IndexOf(split, "MS/MS scan number".ToLowerInvariant());
                _intensityCol = Array.IndexOf(split, "Intensity".ToLowerInvariant());
                _mzCol = Array.IndexOf(split, "m/z".ToLowerInvariant());

                return PsmFileType.MaxQuantEvidence;
            }

            // MaxQuant MS/MS input
            else if (split.Contains("Raw file".ToLowerInvariant())
                && split.Contains("Sequence".ToLowerInvariant())
                && split.Contains("Modified sequence".ToLowerInvariant())
                && split.Contains("Mass".ToLowerInvariant())
                && split.Contains("Retention time".ToLowerInvariant())
                && split.Contains("Charge".ToLowerInvariant())
                && split.Contains("Proteins".ToLowerInvariant()))
            {
                _fileNameCol = Array.IndexOf(split, "Raw file".ToLowerInvariant());
                _baseSequCol = Array.IndexOf(split, "Sequence".ToLowerInvariant());
                _fullSequCol = Array.IndexOf(split, "Modified sequence".ToLowerInvariant());
                _monoMassCol = Array.IndexOf(split, "Mass".ToLowerInvariant());
                _msmsRetnCol = Array.IndexOf(split, "Retention time".ToLowerInvariant());
                _chargeStCol = Array.IndexOf(split, "Charge".ToLowerInvariant());
                _protNameCol = Array.IndexOf(split, "Proteins".ToLowerInvariant());
                _geneNameCol = Array.IndexOf(split, "Gene Names".ToLowerInvariant());

                _organismCol = Array.IndexOf(split, "Organism Name".ToLowerInvariant());

                return PsmFileType.MaxQuant;
            }

            // Peptide Shaker Input
            else if (split.Contains("Spectrum File".ToLowerInvariant())
                && split.Contains("Sequence".ToLowerInvariant())
                && split.Contains("Modified Sequence".ToLowerInvariant())
                && split.Contains("Theoretical Mass".ToLowerInvariant())
                && split.Contains("RT".ToLowerInvariant())
                && split.Contains("Identification Charge".ToLowerInvariant())
                && split.Contains("Protein(s)".ToLowerInvariant()))
            {
                _fileNameCol = Array.IndexOf(split, "Spectrum File".ToLowerInvariant());
                _baseSequCol = Array.IndexOf(split, "Sequence".ToLowerInvariant());
                _fullSequCol = Array.IndexOf(split, "Modified Sequence".ToLowerInvariant());
                _monoMassCol = Array.IndexOf(split, "Theoretical Mass".ToLowerInvariant());
                _msmsRetnCol = Array.IndexOf(split, "RT".ToLowerInvariant());
                _chargeStCol = Array.IndexOf(split, "Identification Charge".ToLowerInvariant());
                _protNameCol = Array.IndexOf(split, "Protein(s)".ToLowerInvariant());

                _geneNameCol = Array.IndexOf(split, "Gene Name".ToLowerInvariant()); // probably doesn't exist
                _organismCol = Array.IndexOf(split, "Organism Name".ToLowerInvariant());

                return PsmFileType.PeptideShaker;
            }

            // Percolator Input
            // Assume that no decoy are provided in this input
            else if (split.Contains("file_idx".ToLowerInvariant())
                && split.Contains("scan".ToLowerInvariant())
                && split.Contains("charge".ToLowerInvariant())
                && split.Contains("spectrum neutral mass".ToLowerInvariant()) //experimental neutral mass
                && split.Contains("peptide mass".ToLowerInvariant()) //theoretical neutral (uncharged) peptide mass
                && split.Contains("sequence".ToLowerInvariant())
                && split.Contains("protein id".ToLowerInvariant()))
            {
                _fileNameCol = Array.IndexOf(split, "file_idx".ToLowerInvariant());
                _fullSequCol = Array.IndexOf(split, "sequence".ToLowerInvariant());
                _monoMassCol = Array.IndexOf(split, "peptide mass".ToLowerInvariant()); //TODO: see if this needs to be theoretical or experimental mass AND if it is neutral or monoisotopic(H+)
                _msmsScanCol = Array.IndexOf(split, "scan".ToLowerInvariant());
                _chargeStCol = Array.IndexOf(split, "charge".ToLowerInvariant());
                _protNameCol = Array.IndexOf(split, "protein id".ToLowerInvariant());
                _qValueCol = Array.IndexOf(split, "percolator q-value".ToLowerInvariant());

                return PsmFileType.Percolator;
            }

            // Generic MS/MS input
            if (split.Contains("File Name".ToLowerInvariant())
                        && split.Contains("Base Sequence".ToLowerInvariant())
                        && split.Contains("Full Sequence".ToLowerInvariant())
                        && split.Contains("Peptide Monoisotopic Mass".ToLowerInvariant())
                        && split.Contains("Scan Retention Time".ToLowerInvariant())
                        && split.Contains("Precursor Charge".ToLowerInvariant())
                        && split.Contains("Protein Accession".ToLowerInvariant()))
            {
                _fileNameCol = Array.IndexOf(split, "File Name".ToLowerInvariant());
                _baseSequCol = Array.IndexOf(split, "Base Sequence".ToLowerInvariant());
                _fullSequCol = Array.IndexOf(split, "Full Sequence".ToLowerInvariant());
                _monoMassCol = Array.IndexOf(split, "Peptide Monoisotopic Mass".ToLowerInvariant());
                _msmsRetnCol = Array.IndexOf(split, "Scan Retention Time".ToLowerInvariant());
                _chargeStCol = Array.IndexOf(split, "Precursor Charge".ToLowerInvariant());
                _protNameCol = Array.IndexOf(split, "Protein Accession".ToLowerInvariant());

                _geneNameCol = Array.IndexOf(split, "Gene Name".ToLowerInvariant()); // probably doesn't exist
                _organismCol = Array.IndexOf(split, "Organism Name".ToLowerInvariant());

                return PsmFileType.Generic;
            }

            return type;
        }

    }

    public class ScanHeaderInfo
    {
        public ScanHeaderInfo(string fullFilePathWithExtension, string filename, int scanNumber, double retentionTime)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            FileNameWithoutExtension = filename;
            ScanNumber = scanNumber;
            RetentionTime = retentionTime;
        }
        public string FullFilePathWithExtension { get; private set; }
        public string FileNameWithoutExtension { get; private set; }
        public int ScanNumber { get; private set; }
        public double RetentionTime { get; private set; }
    }

    public class ScanInfoRecovery
    {
        private enum DataFileType
        { Thermo, mzML, mgf, unknown }

        public static List<ScanHeaderInfo> FileScanHeaderInfo(string fullFilePathWithExtension)
        {
            string filename = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(fullFilePathWithExtension);
            List<ScanHeaderInfo> scanHeaderInfoList = new();
            switch (GetDataFileType(fullFilePathWithExtension))
            {
                case DataFileType.Thermo:
                    ThermoRawFileReader staticRaw = ThermoRawFileReader.LoadAllStaticData(fullFilePathWithExtension);
                    foreach (MsDataScan item in staticRaw)
                    {
                        scanHeaderInfoList.Add(new ScanHeaderInfo(fullFilePathWithExtension, filename, item.OneBasedScanNumber, item.RetentionTime));
                    }
                    break;
                case DataFileType.mzML:
                    List<MsDataScan> mzmlDataScans = Mzml.LoadAllStaticData(fullFilePathWithExtension).GetAllScansList();
                    foreach (MsDataScan item in mzmlDataScans)
                    {
                        scanHeaderInfoList.Add(new ScanHeaderInfo(fullFilePathWithExtension, filename, item.OneBasedScanNumber, item.RetentionTime));
                    }
                    break;
                case DataFileType.mgf:
                    List<MsDataScan> mgfDataScans = Mgf.LoadAllStaticData(fullFilePathWithExtension).GetAllScansList();
                    foreach (MsDataScan item in mgfDataScans)
                    {
                        scanHeaderInfoList.Add(new ScanHeaderInfo(fullFilePathWithExtension, filename, item.OneBasedScanNumber, item.RetentionTime));
                    }
                    break;
                case DataFileType.unknown:
                default:
                    break;
            }
            return scanHeaderInfoList;
        }
        private static DataFileType GetDataFileType(string fullFilePathWithExtension)
        {
            string g = Path.GetExtension(fullFilePathWithExtension).ToLowerInvariant();
            return Path.GetExtension(fullFilePathWithExtension).ToLowerInvariant() switch
            {
                ".raw" => DataFileType.Thermo,
                ".mzml" => DataFileType.mzML,
                ".mgf" => DataFileType.mgf,
                _ => DataFileType.unknown,
            };
        }
    }
}
