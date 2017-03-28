using EngineLayer;
using EngineLayer.Analysis;
using MathNet.Numerics.Distributions;
using Chemistry;
using Spectra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Xml;
using MzIdentML;
using System.Xml.Serialization;

namespace TaskLayer
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
        public List<DbForTask> dbFilenameList;

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

        public static event EventHandler<StringEventArgs> StartingDataFileHandler;

        public static event EventHandler<StringEventArgs> FinishedDataFileHandler;

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
            sb.AppendLine(string.Join(Environment.NewLine, dbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FileName)));
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

        protected internal static void MatchXMLmodsToKnownMods(List<DbForTask> listOfDbs, List<MetaMorpheusModification> modsKnown, out Dictionary<string, List<MetaMorpheusModification>> modsToLocalize, out HashSet<string> modsInXMLtoTrim)
        {
            modsToLocalize = new Dictionary<string, List<MetaMorpheusModification>>();
            var modsInXML = ReadXmlModifications(listOfDbs.Select(b => b.FileName).Where(r => (!r.Contains(".fasta"))));
            modsInXMLtoTrim = new HashSet<string>(modsInXML);
            foreach (var knownMod in modsKnown)
                if (modsInXML.Contains(knownMod.NameInXml))
                {
                    if (modsToLocalize.ContainsKey(knownMod.NameInXml))
                        modsToLocalize[knownMod.NameInXml].Add(knownMod);
                    else
                        modsToLocalize.Add(knownMod.NameInXml, new List<MetaMorpheusModification> { knownMod });
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

        protected internal void WritePsmsToMzIdentmL(List<NewPsmWithFdr> items, string outputFolder, string fileName)
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML.Generated.MzIdentMLType));
            var _mzid = new mzIdentML.Generated.MzIdentMLType()
            {
                version = "1.2.0"
                
            };

            //cvlist: URLs of controlled vocabularies used within the file. add others?
            _mzid.cvList = new mzIdentML.Generated.cvType[2] { new mzIdentML.Generated.cvType()
            {
                id = "PSI-MS",
                fullName = "Proteomics Standards Initiative Mass Spectrometry Vocabularies",
                uri = "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",
                version= "4.0.9"
            }
            , new mzIdentML.Generated.cvType()
            {
                id = "MetaMorpheus",
                fullName = "MetaMorpheus controlled vocabulary", //used for modification names and MM scores
                uri ="https://github.com/smith-chem-wisc/MetaMorpheus",
            }
            }; 

            //analysis software list:software packages used
            _mzid.AnalysisSoftwareList = new mzIdentML.Generated.AnalysisSoftwareType[1] {  new mzIdentML.Generated.AnalysisSoftwareType()
            {
                id = "MetaMorpheus", 
                name = "MetaMorpheus",
                version = MyEngine.MetaMorpheusVersion,
                uri = "https://github.com/smith-chem-wisc/MetaMorpheus"
            }};

            //analysis sample collection: biological smpales analysed, annotated with CV terms
            _mzid.AnalysisCollection = new mzIdentML.Generated.AnalysisCollectionType();

            //sequence collection: database entries of protein/peptide sequences identified and modifications
            _mzid.SequenceCollection = new mzIdentML.Generated.SequenceCollectionType();
            _mzid.SequenceCollection.Peptide = new mzIdentML.Generated.PeptideType[items.SelectMany(p => p.thisPSM.peptidesWithSetModifications).Distinct().Count()];
            _mzid.SequenceCollection.DBSequence = new mzIdentML.Generated.DBSequenceType[items.SelectMany(p => p.thisPSM.peptidesWithSetModifications.Select(t => t.Protein)).Distinct().Count()];
            _mzid.SequenceCollection.PeptideEvidence = new mzIdentML.Generated.PeptideEvidenceType[items.Count()];
            int protein_index = 0; // protein overall
            int peptide_index = 0; // peptide overall
            int peptide_envidence = 0; // peptide evidence overall
            int protein_peptide_index = 0; // peptides for each protein
            int mod_count = 0; // mod count for each peptide
            Dictionary<Peptide, string> peptideEvidenceID = new Dictionary<Peptide, string>(); //key = peptide, value = peptide evidence id

            foreach (Protein protein in items.SelectMany(p => p.thisPSM.peptidesWithSetModifications).Select(p => p.Protein).Distinct())
            {
                protein_peptide_index = 0;
                string database_sequence_ref = "DBSeq_" + protein.Accession;
                _mzid.SequenceCollection.DBSequence[protein_index] = new mzIdentML.Generated.DBSequenceType()
                {
                    id = database_sequence_ref,
                    length = protein.Length,
                    searchDatabase_ref = "Uniprot",
                    accession = protein.Accession,
                    Seq = protein.BaseSequence,
                    cvParam = new mzIdentML.Generated.CVParamType[1] { new mzIdentML.Generated.CVParamType()
                    {
                        accession = "MS:1001088",
                        name = "protein description",
                        cvRef = "PSI-MS",
                        value = protein.FullDescription
                    }
                    }
                };
                foreach (PeptideWithSetModifications peptide in items.Select(p => p.thisPSM.peptidesWithSetModifications.First()).Where(p => p.Protein == protein).Distinct().ToList())
                {
                    string peptide_id = "peptide_" + protein_index + "_" + protein_peptide_index;
                    _mzid.SequenceCollection.PeptideEvidence[peptide_envidence] = new mzIdentML.Generated.PeptideEvidenceType()
                    {
                        isDecoy = peptide.Protein.IsDecoy,
                        peptide_ref = peptide_id,
                        id = "PE_" + protein_index + "_" + protein_peptide_index + "_" + peptide.Protein.Accession,
                        startSpecified = true,
                        endSpecified = true,
                        start = peptide.OneBasedStartResidueInProtein,
                        end = peptide.OneBasedEndResidueInProtein,
                        pre = peptide.PreviousAminoAcid.ToString(),
                        post = (peptide.OneBasedEndResidueInProtein < protein.BaseSequence.Length) ? protein[peptide.OneBasedEndResidueInProtein].ToString() : "-",
                        dBSequence_ref = database_sequence_ref
                    };
                    peptideEvidenceID.Add(peptide, "PE_" + protein_index + "_" + protein_peptide_index + "_" + peptide.Protein.Accession);
                    peptide_envidence++;

                    _mzid.SequenceCollection.Peptide[peptide_index] = new mzIdentML.Generated.PeptideType()
                    {
                        id = peptide_id,
                        PeptideSequence = peptide.BaseSequence,
                        Modification = new mzIdentML.Generated.ModificationType[peptide.allModsOneIsNterminus.Count]
                    };
                    mod_count = 0;
                    for (int position = 0; position < peptide.BaseSequence.Length + 1; position++)
                    {
                        if (peptide.allModsOneIsNterminus.ContainsKey(position))
                        {
                            MetaMorpheusModification mod = peptide.allModsOneIsNterminus[position];
                            _mzid.SequenceCollection.Peptide[peptide_index].Modification[mod_count] = new mzIdentML.Generated.ModificationType()
                            {
                                locationSpecified = true,
                                location = position,
                                monoisotopicMassDeltaSpecified = true,
                                monoisotopicMassDelta = mod.ObservedMassShift,
                                cvParam = new mzIdentML.Generated.CVParamType[1] { new mzIdentML.Generated.CVParamType()
                            {
                                cvRef = "MetaMorpheus", 
                                name = "MetaMorpheusModification",
                                accession = mod.NameInXml
                            } }
                            };
                            mod_count++;
                        }
                    }
                    peptide_index++;
                    protein_peptide_index++;
                }
                protein_index++;
            }

            //Analysis collection: the analyses performed to get the results, which map the input and output data sets. 
            _mzid.AnalysisCollection = new mzIdentML.Generated.AnalysisCollectionType()
            {
                SpectrumIdentification = new mzIdentML.Generated.SpectrumIdentificationType[1] { new mzIdentML.Generated.SpectrumIdentificationType()
                {
                    id = "SI",
                    spectrumIdentificationProtocol_ref = "SIP",
                    spectrumIdentificationList_ref = "SIL",
                    InputSpectra = new mzIdentML.Generated.InputSpectraType[items.Select(p => p.thisPSM.newPsm.fileName).Distinct().Count()]
                } }
            };
            int file_index = 0;
            foreach (string filename in items.Select(p => p.thisPSM.newPsm.fileName).Distinct().OrderBy(p => p).ToList())
            {
                _mzid.AnalysisCollection.SpectrumIdentification[0].InputSpectra[file_index] = new mzIdentML.Generated.InputSpectraType()
                {
                    spectraData_ref = "SD_" + file_index
                };
            }
            //Analysis protocol collection: parameters
            _mzid.AnalysisProtocolCollection = new mzIdentML.Generated.AnalysisProtocolCollectionType();
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol = new mzIdentML.Generated.SpectrumIdentificationProtocolType[1] { new mzIdentML.Generated.SpectrumIdentificationProtocolType() };
            //add later? - already in params.txt file, may want here also

            //DataColection
            _mzid.DataCollection = new mzIdentML.Generated.DataCollectionType();
            _mzid.DataCollection.AnalysisData = new mzIdentML.Generated.AnalysisDataType();
            _mzid.DataCollection.Inputs = new mzIdentML.Generated.InputsType();
            _mzid.DataCollection.Inputs.SpectraData = new mzIdentML.Generated.SpectraDataType[1]
            { new mzIdentML.Generated.SpectraDataType()
            {
             FileFormat = new mzIdentML.Generated.FileFormatType()
             {
                cvParam = new mzIdentML.Generated.CVParamType()
                {
                    name = "mzML format"
                }
             }
            }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList = new mzIdentML.Generated.SpectrumIdentificationListType[items.Select(p => p.thisPSM.newPsm.fileName).Distinct().Count()];

            file_index = 0; //file index
            foreach (string filename in items.Select(p => p.thisPSM.newPsm.fileName).Distinct())
            {
                int spectrum_id_index = 0;
                _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[file_index] = new mzIdentML.Generated.SpectrumIdentificationListType()
                { SpectrumIdentificationResult = new mzIdentML.Generated.SpectrumIdentificationResultType[items.Where(p => p.thisPSM.newPsm.fileName == filename).Select(p => p.thisPSM.newPsm.scanNumber).Distinct().Count()] };

                foreach (int scan in items.Where(p => p.thisPSM.newPsm.fileName == filename).Select(p => p.thisPSM.newPsm.scanNumber).Distinct())
                {
                    //spectrum identification result: all identifications made from searching one spectrum (could be more than one if DIA or coisolation)
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[file_index].SpectrumIdentificationResult[spectrum_id_index] = new mzIdentML.Generated.SpectrumIdentificationResultType()
                    {
                        spectrumID = ("scan=" + scan).ToString(),
                        SpectrumIdentificationItem = new mzIdentML.Generated.SpectrumIdentificationItemType[items.Where(p => p.thisPSM.newPsm.fileName == filename && p.thisPSM.newPsm.scanNumber == scan).Count()]
                    };

                    int psm_index = 0;
                    //spectrum identification item: one PSM
                    foreach (NewPsmWithFdr psm in items.Where(p => p.thisPSM.newPsm.fileName == filename && p.thisPSM.newPsm.scanNumber == scan))
                    {
                        _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[file_index].SpectrumIdentificationResult[spectrum_id_index].SpectrumIdentificationItem[psm_index] = new mzIdentML.Generated.SpectrumIdentificationItemType()
                        {
                            experimentalMassToCharge = psm.thisPSM.newPsm.scanPrecursorMZ,
                            calculatedMassToCharge = psm.thisPSM.PeptideMonoisotopicMass.ToMassToChargeRatio(psm.thisPSM.newPsm.scanPrecursorCharge),
                            calculatedMassToChargeSpecified = true,
                            chargeState = psm.thisPSM.newPsm.scanPrecursorCharge,
                            cvParam = new mzIdentML.Generated.CVParamType[2] { new mzIdentML.Generated.CVParamType()
                         {
                         name = "Score",
                        value = psm.thisPSM.Score.ToString(),
                        cvRef = "MetaMorpheus"
                        },
                        new mzIdentML.Generated.CVParamType()
                                {
                                    accession = "MS:1001868",
                                    name = "distinct peptide-level q-value",
                                    cvRef = "PSI-MS",
                                    value = psm.qValue.ToString()
                                }
                         },
                            PeptideEvidenceRef = new mzIdentML.Generated.PeptideEvidenceRefType[1] { new mzIdentML.Generated.PeptideEvidenceRefType()
                        {
                            peptideEvidence_ref = peptideEvidenceID[psm.thisPSM.peptidesWithSetModifications.First()]
                        }
                        },
                            Fragmentation = new mzIdentML.Generated.IonTypeType[psm.thisPSM.newPsm.matchedIonsList.Count]
                        };
                        int ion_list_count = 0;
                        foreach (var ion_list in psm.thisPSM.newPsm.matchedIonsList)
                        {
                            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[file_index].SpectrumIdentificationResult[spectrum_id_index].SpectrumIdentificationItem[psm_index].Fragmentation[ion_list_count] = new mzIdentML.Generated.IonTypeType()
                            {
                                cvParam = new mzIdentML.Generated.CVParamType()
                                {
                                    name = "frag: " + ion_list.Key
                                },
                                FragmentArray = new mzIdentML.Generated.FragmentArrayType[1] { new mzIdentML.Generated.FragmentArrayType() { values = new float[ion_list.Value.Length] ,
                               measure_ref = "m_mass",
                         } }
                            };

                            for (int i = 0; i < ion_list.Value.Length; i++)
                            {
                                _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[file_index].SpectrumIdentificationResult[spectrum_id_index].SpectrumIdentificationItem[psm_index].Fragmentation[ion_list_count].FragmentArray[0].values[i] = Math.Abs((float)ion_list.Value[i]);
                            }
                            ion_list_count++;
                        }
                        psm_index++;
                    }
                    spectrum_id_index++;
                }
                file_index++;
            }

            TextWriter writer = new StreamWriter(Path.Combine(outputFolder, fileName + ".mzid"));
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected IEnumerable<Protein> GetProteins(bool onTheFlyDecoys, IDictionary<string, List<MetaMorpheusModification>> allModifications, DbForTask dbForTask)
        {
            using (var stream = new FileStream(dbForTask.FileName, FileMode.Open))
            {
                Stream uniprotXmlFileStream = stream;
                if (dbForTask.FileName.EndsWith(".gz"))
                    uniprotXmlFileStream = new GZipStream(stream, CompressionMode.Decompress);

                string[] nodes = new string[6];

                string accession = null;
                string name = null;
                string full_name = null;
                string sequence = null;
                string feature_type = null;
                string feature_description = null;
                int oneBasedfeature_position = -1;
                int oneBasedbeginPosition = -1;
                int oneBasedendPosition = -1;
                var oneBasedBeginPositions = new List<int>();
                var oneBasedEndPositions = new List<int>();
                var peptideTypes = new List<string>();
                var oneBasedModifications = new Dictionary<int, List<MetaMorpheusModification>>();
                int offset = 0;

                // xml db
                if (!dbForTask.FileName.EndsWith(".fasta"))
                {
                    using (XmlReader xml = XmlReader.Create(uniprotXmlFileStream))
                    {
                        while (xml.Read())
                        {
                            switch (xml.NodeType)
                            {
                                case XmlNodeType.Element:
                                    nodes[xml.Depth] = xml.Name;
                                    switch (xml.Name)
                                    {
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
                                                List<MetaMorpheusModification> residue_modifications;
                                                if (!oneBasedModifications.TryGetValue(oneBasedfeature_position, out residue_modifications))
                                                {
                                                    residue_modifications = new List<MetaMorpheusModification>();
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
                                                var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false, dbForTask.IsContaminant);

                                                yield return protein;

                                                offset += protein.Length;
                                                if (onTheFlyDecoys)
                                                {
                                                    char[] sequence_array = sequence.ToCharArray();
                                                    Dictionary<int, List<MetaMorpheusModification>> decoy_modifications = null;
                                                    if (sequence.StartsWith("M", StringComparison.InvariantCulture))
                                                    {
                                                        // Do not include the initiator methionine in reversal!!!
                                                        Array.Reverse(sequence_array, 1, sequence.Length - 1);
                                                        if (oneBasedModifications != null)
                                                        {
                                                            decoy_modifications = new Dictionary<int, List<MetaMorpheusModification>>(oneBasedModifications.Count);
                                                            foreach (KeyValuePair<int, List<MetaMorpheusModification>> kvp in oneBasedModifications)
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
                                                            decoy_modifications = new Dictionary<int, List<MetaMorpheusModification>>(oneBasedModifications.Count);
                                                            foreach (KeyValuePair<int, List<MetaMorpheusModification>> kvp in oneBasedModifications)
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
                                                    var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true, dbForTask.IsContaminant);
                                                    yield return decoy_protein;
                                                    offset += protein.Length;
                                                }
                                            }
                                            accession = null;
                                            name = null;
                                            full_name = null;
                                            sequence = null;
                                            feature_type = null;
                                            feature_description = null;
                                            oneBasedfeature_position = -1;
                                            oneBasedModifications = new Dictionary<int, List<MetaMorpheusModification>>();

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

                // fasta db
                else
                {
                    StreamReader fasta = new StreamReader(stream);

                    while (true)
                    {
                        string line = fasta.ReadLine();

                        if (line.StartsWith(">"))
                        {
                            // fasta protein only has accession, fullname, sequence (no mods)
                            string[] delimiters = { ">", "|", " OS=" };
                            string[] output = line.Split(delimiters, StringSplitOptions.None);
                            if (output.Length > 4)
                            {
                                accession = output[2];
                                name = accession;
                                full_name = output[3];
                            }
                            else
                            {
                                // can't read protein description
                                full_name = line.Substring(1);
                                accession = "";
                            }

                            // new protein
                            sequence = "";
                        }
                        else
                        {
                            sequence += line.Trim();
                        }

                        if (fasta.Peek() == '>' || fasta.Peek() == -1)
                        {
                            if (accession != null && sequence != null)
                            {
                                var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false, dbForTask.IsContaminant);
                                yield return protein;

                                if (onTheFlyDecoys)
                                {
                                    char[] sequence_array = sequence.ToCharArray();
                                    Dictionary<int, List<MetaMorpheusModification>> decoy_modifications = null;
                                    if (sequence.StartsWith("M", StringComparison.InvariantCulture))
                                    {
                                        // Do not include the initiator methionine in reversal!!!
                                        Array.Reverse(sequence_array, 1, sequence.Length - 1);
                                        if (oneBasedModifications != null)
                                        {
                                            decoy_modifications = new Dictionary<int, List<MetaMorpheusModification>>(oneBasedModifications.Count);
                                            foreach (KeyValuePair<int, List<MetaMorpheusModification>> kvp in oneBasedModifications)
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
                                            decoy_modifications = new Dictionary<int, List<MetaMorpheusModification>>(oneBasedModifications.Count);
                                            foreach (KeyValuePair<int, List<MetaMorpheusModification>> kvp in oneBasedModifications)
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
                                    var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true, dbForTask.IsContaminant);
                                    yield return decoy_protein;
                                    offset += protein.Length;
                                }
                            }
                        }

                        // no input left
                        if (fasta.Peek() == -1)
                        {
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
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

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
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tNtermLocFrac\tCtermLocFrac\tFracWithSingle\tOverlappingFrac\tUniprot");
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
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile);
        }

        protected void SucessfullyFinishedWritingFile(string path)
        {
            FinishedWritingFileHandler?.Invoke(this, new SingleFileEventArgs(path));
        }

        protected void StartingDataFile(string v)
        {
            StartingDataFileHandler?.Invoke(this, new StringEventArgs(v));
        }

        protected void FinishedDataFile(string v)
        {
            FinishedDataFileHandler?.Invoke(this, new StringEventArgs(v));
        }

        #endregion Protected Methods

        #region Private Methods

        private static HashSet<string> ReadXmlModifications(IEnumerable<string> uniProtXmlProteomeDatabaseFilepaths)
        {
            var modifications_in_database = new HashSet<string>();
            foreach (var uniProtXmlProteomeDatabaseFilepath in uniProtXmlProteomeDatabaseFilepaths)
                using (var stream = new FileStream(uniProtXmlProteomeDatabaseFilepath, FileMode.Open))
                {
                    Stream uniprotXmlFileStream = stream;
                    if (uniProtXmlProteomeDatabaseFilepath.EndsWith(".gz", StringComparison.OrdinalIgnoreCase))
                        uniprotXmlFileStream = new GZipStream(stream, CompressionMode.Decompress);
                    using (XmlReader xml = XmlReader.Create(uniprotXmlFileStream))
                        while (xml.ReadToFollowing("feature"))
                            if (xml.GetAttribute("type") == "modified residue")
                            {
                                string description = xml.GetAttribute("description");
                                if (!description.Contains("variant"))
                                {
                                    int semicolon_index = description.IndexOf(';');
                                    if (semicolon_index >= 0)
                                        description = description.Substring(0, semicolon_index);
                                    modifications_in_database.Add(description);
                                }
                            }
                }
            return modifications_in_database;
        }

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