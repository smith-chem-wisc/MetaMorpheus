using InternalLogicEngineLayer;
using MathNet.Numerics.Distributions;
using OldInternalLogic;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Xml;

namespace InternalLogicTaskLayer
{
    public enum MyTaskEnum
    {
        Search,
        GPTMD,
        Calibrate
    }

    public abstract class MyTaskEngine : MyEngine
    {

        #region Public Fields

        public List<string> rawDataFilenameList;
        public List<string> xmlDbFilenameList;

        #endregion Public Fields

        #region Protected Constructors

        protected MyTaskEngine() : base(1)
        {
        }

        #endregion Protected Constructors

        #region Public Events

        public static event EventHandler<SingleTaskEventArgs> finishedSingleTaskHandler;

        public static event EventHandler<SingleFileEventArgs> finishedWritingFileHandler;

        public static event EventHandler<SingleTaskEventArgs> startingSingleTaskHander;

        #endregion Public Events

        #region Public Properties

        public MyTaskEnum taskType { get; internal set; }
        public bool bIons { get; set; }
        public InitiatorMethionineBehavior initiatorMethionineBehavior { get; set; }
        public bool IsMySelected { get; set; }
        public int maxMissedCleavages { get; set; }
        public int maxModificationIsoforms { get; set; }
        public string output_folder { get; set; }
        public Protease protease { get; set; }
        public bool yIons { get; set; }
        public int maxNumPeaksPerScan { get; set; }

        #endregion Public Properties

        #region Public Methods

        public new MyResults Run()
        {
            startingSingleTask();
            var paramsFileName = Path.Combine(output_folder, "params.txt");
            using (StreamWriter file = new StreamWriter(paramsFileName))
            {
                if (MyEngine.MetaMorpheusVersion.Equals("1.0.0.0"))
                    file.WriteLine("MetaMorpheus: Not a release version");
                else
                    file.WriteLine("MetaMorpheus: version " + MyEngine.MetaMorpheusVersion);
                file.Write(ToString());
            }
            SucessfullyFinishedWritingFile(paramsFileName);
            var heh = base.Run();
            var resultsFileName = Path.Combine(output_folder, "results.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                if (MyEngine.MetaMorpheusVersion.Equals("1.0.0.0"))
                    file.WriteLine("MetaMorpheus: Not a release version");
                else
                    file.WriteLine("MetaMorpheus: version " + MyEngine.MetaMorpheusVersion);
                file.Write(heh.ToString());
            }
            SucessfullyFinishedWritingFile(resultsFileName);
            finishedSingleTask();
            return heh;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(GetSpecificTaskInfo());
            sb.AppendLine(taskType.ToString());
            sb.AppendLine("Spectra files:");
            sb.AppendLine(string.Join(Environment.NewLine, rawDataFilenameList.Select(b => '\t' + b)));
            sb.AppendLine("XML files:");
            sb.AppendLine(string.Join(Environment.NewLine, xmlDbFilenameList.Select(b => '\t' + b)));
            sb.AppendLine("initiatorMethionineBehavior: " + initiatorMethionineBehavior);
            sb.AppendLine("maxMissedCleavages: " + maxMissedCleavages);
            sb.AppendLine("maxModificationIsoforms: " + maxModificationIsoforms);
            sb.AppendLine("output_folder: " + output_folder);
            sb.AppendLine("protease: " + protease);
            sb.AppendLine("bIons: " + bIons);
            sb.Append("yIons: " + yIons);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected static void MatchXMLmodsToKnownMods(List<string> listOfXMLdbs, List<MorpheusModification> modsKnown, out Dictionary<string, List<MorpheusModification>> modsToLocalize, out HashSet<string> modsInXMLtoTrim)
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
                        modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification> { knownMod });
                    modsInXMLtoTrim.Remove(knownMod.NameInXML);
                }
        }

        protected IEnumerable<Protein> getProteins(bool onTheFlyDecoys, IDictionary<string, List<MorpheusModification>> allModifications, string FileName)
        {
            using (XmlReader xml = XmlReader.Create(FileName))
            {
                string[] nodes = new string[6];

                string dataset = null;
                string accession = null;
                string name = null;
                string full_name = null;
                string organism = null;
                string gene_name = null;
                string sequence = null;
                string feature_type = null;
                string feature_description = null;
                int oneBasedfeature_position = -1;
                int oneBasedbeginPosition = -1;
                int oneBasedendPosition = -1;
                var oneBasedBeginPositions = new List<int>();
                var oneBasedEndPositions = new List<int>();
                var peptideTypes = new List<string>();
                var oneBasedModifications = new Dictionary<int, List<MorpheusModification>>();
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
                                    catch (ArgumentNullException)
                                    {
                                    }
                                    break;

                                case "end":
                                    try
                                    {
                                        oneBasedendPosition = int.Parse(xml.GetAttribute("position"));
                                    }
                                    catch (ArgumentNullException)
                                    {
                                    }
                                    break;

                                case "sequence":
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

                                    if (accession != null && sequence != null)
                                    {
                                        var protein = new Protein(sequence, accession, dataset_abbreviation, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false);

                                        yield return protein;

                                        offset += protein.Length;
                                        if (onTheFlyDecoys)
                                        {
                                            char[] sequence_array = sequence.ToCharArray();
                                            Dictionary<int, List<MorpheusModification>> decoy_modifications = null;
                                            if (sequence.StartsWith("M", StringComparison.InvariantCulture))
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
                                            var reversed_sequence = new string(sequence_array);
                                            int[] decoybeginPositions = new int[oneBasedBeginPositions.Count];
                                            int[] decoyendPositions = new int[oneBasedEndPositions.Count];
                                            string[] decoyBigPeptideTypes = new string[oneBasedEndPositions.Count];
                                            for (int i = 0; i < decoybeginPositions.Length; i++)
                                            {
                                                decoybeginPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedEndPositions[i] + 1;
                                                decoyendPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedBeginPositions[i] + 1;
                                                decoyBigPeptideTypes[oneBasedBeginPositions.Count - i - 1] = peptideTypes[i];
                                            }
                                            var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, dataset_abbreviation, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true);
                                            yield return decoy_protein;
                                            offset += protein.Length;
                                        }
                                    }
                                    dataset = null;
                                    accession = null;
                                    name = null;
                                    full_name = null;
                                    organism = null;
                                    gene_name = null;
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

        protected void WritePSMsToTSV(List<NewPsmWithFDR> items, string output_folder, string fileName)
        {
            var writtenFile = Path.Combine(output_folder, fileName + ".psmtsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(NewPsmWithFDR.GetTabSeparatedHeader());
                for (int i = 0; i < items.Count; i++)
                    output.WriteLine(items[i]);
            }
            SucessfullyFinishedWritingFile(writtenFile);
        }

        protected void WriteProteinGroupsToTSV(List<ProteinGroup> items, string output_folder, string fileName)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(output_folder, fileName + ".psmtsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(ProteinGroup.GetTabSeparatedHeader());
                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile);
            }
        }

        protected void WriteTree(BinTreeStructure myTreeStructure, string output_folder, string fileName)
        {
            var writtenFile = Path.Combine(output_folder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tNtermLocFrac\tCtermLocFrac\tUniprot");
                foreach (Bin bin in myTreeStructure.finalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.NlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile);
        }

        protected abstract string GetSpecificTaskInfo();

        protected void SucessfullyFinishedWritingFile(string path)
        {
            finishedWritingFileHandler?.Invoke(this, new SingleFileEventArgs(path));
        }

        #endregion Protected Methods

        #region Private Methods

        private void finishedSingleTask()
        {
            finishedSingleTaskHandler?.Invoke(this, new SingleTaskEventArgs(this));
        }

        private void startingSingleTask()
        {
            startingSingleTaskHander?.Invoke(this, new SingleTaskEventArgs(this));
        }

        #endregion Private Methods

    }
}