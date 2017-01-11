using InternalLogic;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Xml;

namespace InternalLogicWithFileIO
{
    public abstract class MyTaskEngine : MyEngine
    {

        protected List<string> xMLdblist;
        protected List<string> rawDataAndResultslist;

        public static event EventHandler<List<string>> newDbsHandler;

        public static event EventHandler<List<string>> newSpectrasHandler;

        public MyTaskEnum taskType { get; internal set; }
        public bool IsMySelected { get; set; }
        public string output_folder { get; private set; }
        public int maxMissedCleavages { get; set; }
        public Protease protease { get; set; }
        public int maxModificationIsoforms { get; set; }
        public InitiatorMethionineBehavior initiatorMethionineBehavior { get; set; }
        public Tolerance productMassTolerance { get; set; }
        public bool bIons { get; set; }
        public bool yIons { get; set; }

        public List<ModListForSearch> listOfModListsForSearch { get; set; }
        public List<SearchModeFoSearch> searchModes { get; set; }

        protected static void GenerateModsFromStrings(List<string> listOfXMLdbs, List<MorpheusModification> modsKnown, out Dictionary<string, List<MorpheusModification>> modsToLocalize, out HashSet<string> modsInXMLtoTrim)
        {
            modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
            var modsInXML = ProteomeDatabaseReader.ReadXMLmodifications(listOfXMLdbs);
            modsInXMLtoTrim = new HashSet<string>(modsInXML);
            foreach (var knownMod in modsKnown)
                if (modsInXML.Contains(knownMod.NameInXML))
                {
                    if (modsToLocalize.ContainsKey(knownMod.NameInXML))
                        modsToLocalize[knownMod.NameInXML].Add(knownMod);
                    else
                        modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification>() { knownMod });
                    modsInXMLtoTrim.Remove(knownMod.NameInXML);
                }
        }

        public void setOutputFolder(string thisOutputPath)
        {
            this.output_folder = thisOutputPath;
        }

        public IEnumerable<Protein> getProteins(bool onTheFlyDecoys, IDictionary<string, List<MorpheusModification>> allModifications, string FileName)
        {
            using (XmlReader xml = XmlReader.Create(FileName))
            {
                string[] nodes = new string[6];

                string dataset = null;
                string accession = null;
                string name = null;
                string full_name = null;
                bool fragment = false;
                string organism = null;
                string gene_name = null;
                string protein_existence = null;
                string sequence_version = null;
                string sequence = null;
                string feature_type = null;
                string feature_description = null;
                int oneBasedfeature_position = -1;
                int oneBasedbeginPosition = -1;
                int oneBasedendPosition = -1;
                List<int> oneBasedBeginPositions = new List<int>();
                List<int> oneBasedEndPositions = new List<int>();
                List<string> peptideTypes = new List<string>();
                Dictionary<int, List<MorpheusModification>> oneBasedModifications = new Dictionary<int, List<MorpheusModification>>();
                int offset = 0;
                while (xml.Read())
                {
                    switch (xml.NodeType)
                    {
                        case XmlNodeType.Element:
                            nodes[xml.Depth] = xml.Name;
                            switch (xml.Name)
                            {
                                case "entry":
                                    dataset = xml.GetAttribute("dataset");
                                    break;

                                case "accession":
                                    if (accession == null)
                                    {
                                        accession = xml.ReadElementString();
                                    }
                                    break;

                                case "name":
                                    if (xml.Depth == 2)
                                    {
                                        name = xml.ReadElementString();
                                    }
                                    else if (nodes[2] == "gene")
                                    {
                                        if (gene_name == null)
                                        {
                                            gene_name = xml.ReadElementString();
                                        }
                                    }
                                    else if (nodes[2] == "organism")
                                    {
                                        if (organism == null)
                                        {
                                            organism = xml.ReadElementString();
                                        }
                                    }
                                    break;

                                case "fullName":
                                    if (full_name == null)
                                    {
                                        full_name = xml.ReadElementString();
                                    }
                                    break;

                                case "proteinExistence":
                                    protein_existence = xml.GetAttribute("type");
                                    break;

                                case "feature":
                                    feature_type = xml.GetAttribute("type");
                                    feature_description = xml.GetAttribute("description");
                                    break;

                                case "position":
                                    oneBasedfeature_position = int.Parse(xml.GetAttribute("position"));
                                    break;

                                case "begin":
                                    try
                                    {
                                        oneBasedbeginPosition = int.Parse(xml.GetAttribute("position"));
                                    }
                                    catch
                                    {
                                    }
                                    break;

                                case "end":
                                    try
                                    {
                                        oneBasedendPosition = int.Parse(xml.GetAttribute("position"));
                                    }
                                    catch
                                    {
                                    }
                                    break;

                                case "sequence":
                                    fragment = xml.GetAttribute("fragment") != null;
                                    sequence_version = xml.GetAttribute("version");
                                    sequence = xml.ReadElementString().Replace("\n", null);
                                    break;
                            }
                            break;

                        case XmlNodeType.EndElement:
                            switch (xml.Name)
                            {
                                case "feature":
                                    if (feature_type == "modified residue" && allModifications != null && !feature_description.Contains("variant") && allModifications.ContainsKey(feature_description))
                                    {
                                        List<MorpheusModification> residue_modifications;
                                        if (!oneBasedModifications.TryGetValue(oneBasedfeature_position, out residue_modifications))
                                        {
                                            residue_modifications = new List<MorpheusModification>();
                                            oneBasedModifications.Add(oneBasedfeature_position, residue_modifications);
                                        }
                                        int semicolon_index = feature_description.IndexOf(';');
                                        if (semicolon_index >= 0)
                                        {
                                            feature_description = feature_description.Substring(0, semicolon_index);
                                        }
                                        residue_modifications.AddRange(allModifications[feature_description]);
                                    }
                                    else if ((feature_type == "peptide" || feature_type == "propeptide" || feature_type == "chain") && oneBasedbeginPosition >= 0 && oneBasedendPosition >= 0)
                                    {
                                        oneBasedBeginPositions.Add(oneBasedbeginPosition);
                                        oneBasedEndPositions.Add(oneBasedendPosition);
                                        peptideTypes.Add(feature_type);
                                    }
                                    oneBasedbeginPosition = -1;
                                    oneBasedendPosition = -1;

                                    break;

                                case "entry":
                                    string dataset_abbreviation;
                                    if (dataset != null && dataset.Equals("Swiss-Prot", StringComparison.InvariantCultureIgnoreCase))
                                        dataset_abbreviation = "sp";
                                    else
                                        dataset_abbreviation = "uk";

                                    if (accession != null && sequence != null && (!SpecificProteinSelection.enabled || SpecificProteinSelection.ConsiderProtein(accession)))
                                    {
                                        Protein protein = new Protein(sequence, accession, dataset_abbreviation, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false);

                                        yield return protein;

                                        offset += protein.Length;
                                        if (onTheFlyDecoys)
                                        {
                                            char[] sequence_array = sequence.ToCharArray();
                                            Dictionary<int, List<MorpheusModification>> decoy_modifications = null;
                                            if (sequence.StartsWith("M"))
                                            {
                                                // Do not include the initiator methionine in reversal!!!
                                                Array.Reverse(sequence_array, 1, sequence.Length - 1);
                                                if (oneBasedModifications != null)
                                                {
                                                    decoy_modifications = new Dictionary<int, List<MorpheusModification>>(oneBasedModifications.Count);
                                                    foreach (KeyValuePair<int, List<MorpheusModification>> kvp in oneBasedModifications)
                                                    {
                                                        if (kvp.Key == 1)
                                                        {
                                                            decoy_modifications.Add(1, kvp.Value);
                                                        }
                                                        else if (kvp.Key > 1)
                                                        {
                                                            decoy_modifications.Add(sequence.Length - kvp.Key + 2, kvp.Value);
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                Array.Reverse(sequence_array);
                                                if (oneBasedModifications != null)
                                                {
                                                    decoy_modifications = new Dictionary<int, List<MorpheusModification>>(oneBasedModifications.Count);
                                                    foreach (KeyValuePair<int, List<MorpheusModification>> kvp in oneBasedModifications)
                                                    {
                                                        decoy_modifications.Add(sequence.Length - kvp.Key + 1, kvp.Value);
                                                    }
                                                }
                                            }
                                            string reversed_sequence = new string(sequence_array);
                                            int[] decoybeginPositions = new int[oneBasedBeginPositions.Count];
                                            int[] decoyendPositions = new int[oneBasedEndPositions.Count];
                                            string[] decoyBigPeptideTypes = new string[oneBasedEndPositions.Count];
                                            for (int i = 0; i < decoybeginPositions.Length; i++)
                                            {
                                                decoybeginPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedEndPositions[i] + 1;
                                                decoyendPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedBeginPositions[i] + 1;
                                                decoyBigPeptideTypes[oneBasedBeginPositions.Count - i - 1] = peptideTypes[i];
                                            }
                                            Protein decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, dataset_abbreviation, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true);
                                            yield return decoy_protein;
                                            offset += protein.Length;
                                        }
                                    }
                                    dataset = null;
                                    accession = null;
                                    name = null;
                                    full_name = null;
                                    fragment = false;
                                    organism = null;
                                    gene_name = null;
                                    protein_existence = null;
                                    sequence_version = null;
                                    sequence = null;
                                    feature_type = null;
                                    feature_description = null;
                                    oneBasedfeature_position = -1;
                                    oneBasedModifications = new Dictionary<int, List<MorpheusModification>>();

                                    oneBasedBeginPositions = new List<int>();
                                    oneBasedEndPositions = new List<int>();
                                    peptideTypes = new List<string>();
                                    break;
                            }
                            break;
                    }
                }
            }
        }
        public void GetPeptideAndFragmentIndices(out List<CompactPeptide> peptideIndex, out Dictionary<float, List<int>> fragmentIndexDict, List<ModListForSearch> collectionOfModLists, bool doFDRanalysis, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, List<Protein> hm, Protease protease, string output_folder)
        {
            #region Index file names

            string folderName = output_folder;
            StringBuilder indexFileSB = new StringBuilder();
            foreach (var heh in xMLdblist)
                indexFileSB.Append(Path.GetFileNameWithoutExtension(heh));
            if (doFDRanalysis)
                indexFileSB.Append("-WithDecoys");
            if (collectionOfModLists.Where(b => b.Fixed).Count() > 0)
            {
                indexFileSB.Append("-fixed");
                foreach (var heh in collectionOfModLists.Where(b => b.Fixed))
                    indexFileSB.Append("-" + Path.GetFileNameWithoutExtension(heh.FileName));
            }
            if (collectionOfModLists.Where(b => b.Fixed).Count() > 0)
            {
                indexFileSB.Append("-variable");
                foreach (var heh in collectionOfModLists.Where(b => b.Variable))
                    indexFileSB.Append("-" + Path.GetFileNameWithoutExtension(heh.FileName));
            }
            if (collectionOfModLists.Where(b => b.Localize).Count() > 0)
            {
                indexFileSB.Append("-localize");
                foreach (var heh in collectionOfModLists.Where(b => b.Localize))
                    indexFileSB.Append("-" + Path.GetFileNameWithoutExtension(heh.FileName));
            }

            string peptideIndexFile = Path.Combine(folderName, indexFileSB.ToString() + "-peptideIndex.ind");
            string fragmentIndexFile = Path.Combine(folderName, indexFileSB.ToString() + "-fragmentIndex.ind");

            #endregion Index file names

            if (!File.Exists(peptideIndexFile) || !File.Exists(fragmentIndexFile))
            {
                output("Generating indices...");

                IndexEngine indexEngine = new IndexEngine(hm, variableModifications, fixedModifications, localizeableModifications, protease);
                IndexResults indexResults = (IndexResults)indexEngine.Run();
                peptideIndex = indexResults.peptideIndex;
                fragmentIndexDict = indexResults.fragmentIndexDict;

                output("Writing peptide index...");
                writePeptideIndex(peptideIndex, peptideIndexFile);
                output("Writing fragment index...");
                writeFragmentIndexNetSerializer(fragmentIndexDict, fragmentIndexFile);
                output("Done Writing fragment index");
            }
            else
            {
                output("Reading peptide index...");
                peptideIndex = readPeptideIndex(peptideIndexFile);
                output("Reading fragment index...");
                fragmentIndexDict = readFragmentIndexNetSerializer(fragmentIndexFile);
            }
        }

        private IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            foreach (var ok in type.Assembly.GetTypes().Where(t => t.IsSubclassOf(type)))
                yield return ok;
            yield return type;
        }

        internal Dictionary<float, List<int>> readFragmentIndexNetSerializer(string fragmentIndexFile)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();

            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            Dictionary<float, List<int>> newPerson;
            using (var file = File.OpenRead(fragmentIndexFile))
                newPerson = (Dictionary<float, List<int>>)ser.Deserialize(file);

            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = string.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                ts.Hours, ts.Minutes, ts.Seconds,
                ts.Milliseconds / 10);
            output("Time to read fragment index with netSerializer: " + elapsedTime);

            return newPerson;
        }

        internal void writeFragmentIndexNetSerializer(Dictionary<float, List<int>> fragmentIndex, string fragmentIndexFile)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();

            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);

            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = string.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                ts.Hours, ts.Minutes, ts.Seconds,
                ts.Milliseconds / 10);
            output("Time to write fragment index with netserializer: " + elapsedTime);
        }

        internal void writePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFile))
            {
                ser.Serialize(file, peptideIndex);
            }
        }

        internal List<CompactPeptide> readPeptideIndex(string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);
            List<CompactPeptide> newPerson;
            using (var file = File.OpenRead(peptideIndexFile))
            {
                newPerson = (List<CompactPeptide>)ser.Deserialize(file);
            }

            return newPerson;
        }
    }
}