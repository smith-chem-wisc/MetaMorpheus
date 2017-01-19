using InternalLogicEngineLayer;
using MathNet.Numerics.Distributions;
using OldInternalLogic;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Xml;

namespace InternalLogicTaskLayer
{
    public enum MyTask
    {
        Search,
        Gptmd,
        Calibrate
    }

    public abstract class MyTaskEngine : MyEngine
    {

        #region Public Fields

        public List<string> rawDataFilenameList;
        public List<XmlForTask> xmlDbFilenameList;

        #endregion Public Fields

        #region Protected Constructors

        protected MyTaskEngine() : base(1)
        {
        }

        #endregion Protected Constructors

        #region Public Events

        public static event EventHandler<SingleTaskEventArgs> FinishedSingleTaskHandler;

        public static event EventHandler<SingleFileEventArgs> FinishedWritingFileHandler;

        public static event EventHandler<SingleTaskEventArgs> StartingSingleTaskHander;

        #endregion Public Events

        #region Public Properties

        public MyTask TaskType { get; internal set; }
        public bool BIons { get; set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        public bool IsMySelected { get; set; }
        public int MaxMissedCleavages { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public string OutputFolder { get; set; }
        public Protease Protease { get; set; }
        public bool YIons { get; set; }
        public int MaxNumPeaksPerScan { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected abstract string SpecificTaskInfo { get; }

        #endregion Protected Properties

        #region Public Methods

        public new MyResults Run()
        {
            startingSingleTask();
            var paramsFileName = Path.Combine(OutputFolder, "params.txt");
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
            var resultsFileName = Path.Combine(OutputFolder, "results.txt");
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
            sb.AppendLine(SpecificTaskInfo);
            sb.AppendLine(TaskType.ToString());
            sb.AppendLine("Spectra files:");
            sb.AppendLine(string.Join(Environment.NewLine, rawDataFilenameList.Select(b => '\t' + b)));
            sb.AppendLine("XML files:");
            sb.AppendLine(string.Join(Environment.NewLine, xmlDbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FileName)));
            sb.AppendLine("initiatorMethionineBehavior: " + InitiatorMethionineBehavior);
            sb.AppendLine("maxMissedCleavages: " + MaxMissedCleavages);
            sb.AppendLine("maxModificationIsoforms: " + MaxModificationIsoforms);
            sb.AppendLine("output_folder: " + OutputFolder);
            sb.AppendLine("protease: " + Protease);
            sb.AppendLine("bIons: " + BIons);
            sb.Append("yIons: " + YIons);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Internal Methods

        protected internal static void MatchXMLmodsToKnownMods(List<XmlForTask> listOfXMLdbs, List<MorpheusModification> modsKnown, out Dictionary<string, List<MorpheusModification>> modsToLocalize, out HashSet<string> modsInXMLtoTrim)
        {
            modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
            var modsInXML = ProteomeDatabaseReader.ReadXmlModifications(listOfXMLdbs.Select(b => b.FileName));
            modsInXMLtoTrim = new HashSet<string>(modsInXML);
            foreach (var knownMod in modsKnown)
                if (modsInXML.Contains(knownMod.NameInXml))
                {
                    if (modsToLocalize.ContainsKey(knownMod.NameInXml))
                        modsToLocalize[knownMod.NameInXml].Add(knownMod);
                    else
                        modsToLocalize.Add(knownMod.NameInXml, new List<MorpheusModification> { knownMod });
                    modsInXMLtoTrim.Remove(knownMod.NameInXml);
                }
        }

        protected internal void WritePsmsToTsv(List<NewPsmWithFdr> items, string outputFolder, string fileName)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".psmtsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(NewPsmWithFdr.TabSeparatedHeader);
                for (int i = 0; i < items.Count; i++)
                    output.WriteLine(items[i]);
            }
            SucessfullyFinishedWritingFile(writtenFile);
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected IEnumerable<Protein> GetProteins(bool onTheFlyDecoys, IDictionary<string, List<MorpheusModification>> allModifications, XmlForTask xmlForTask)
        {
            using (var stream = new FileStream(xmlForTask.FileName, FileMode.Open))
            {
                Stream uniprotXmlFileStream = stream;
                if (xmlForTask.FileName.EndsWith(".gz"))
                    uniprotXmlFileStream = new GZipStream(stream, CompressionMode.Decompress);
                using (XmlReader xml = XmlReader.Create(uniprotXmlFileStream))
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
                                        if (accession != null && sequence != null)
                                        {
                                            var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false, xmlForTask.IsContaminant);

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
                                                var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true, xmlForTask.IsContaminant);
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
        }

        protected void WriteProteinGroupsToTsv(List<ProteinGroup> items, string outputFolder, string fileName)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".prottsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(ProteinGroup.TabSeparatedHeader);
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
                foreach (Bin bin in myTreeStructure.FinalBins.OrderByDescending(b => b.Count))
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
                        + "\t" + bin.Mine
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

        protected void SucessfullyFinishedWritingFile(string path)
        {
            FinishedWritingFileHandler?.Invoke(this, new SingleFileEventArgs(path));
        }

        #endregion Protected Methods

        #region Private Methods

        private void finishedSingleTask()
        {
            FinishedSingleTaskHandler?.Invoke(this, new SingleTaskEventArgs(this));
        }

        private void startingSingleTask()
        {
            StartingSingleTaskHander?.Invoke(this, new SingleTaskEventArgs(this));
        }

        #endregion Private Methods

    }
}