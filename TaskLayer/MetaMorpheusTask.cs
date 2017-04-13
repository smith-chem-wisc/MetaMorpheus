using EngineLayer;
using EngineLayer.Analysis;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml.Serialization;
using UsefulProteomicsDatabases;
using Chemistry;

namespace TaskLayer
{
    public enum MyTask
    {
        Search,
        Gptmd,
        Calibrate
    }

    public abstract class MetaMorpheusTask
    {

        #region Public Fields

        public static readonly TomlConfig tomlConfig = TomlConfig.Create(cfg => cfg
                        .ConfigureType<Tolerance>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => new Tolerance(tmlString.Value))))
                        .ConfigureType<SearchMode>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => MetaMorpheusTask.ParseSearchMode(tmlString.Value))))
                        .ConfigureType<Protease>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => GlobalTaskLevelSettings.ProteaseDictionary[tmlString.Value])))
                        .ConfigureType<List<Tuple<string, string>>>(type => type
                             .WithConversionFor<TomlTableArray>(convert => convert
                             .FromToml(tml => tml.Items.Select(b => new Tuple<string, string>(b.Values.First().Get<string>(), b.Values.Last().Get<string>())).ToList()))));

        #endregion Public Fields

        #region Protected Fields

        protected MyTaskResults myTaskResults;

        #endregion Protected Fields

        #region Public Constructors

        public MetaMorpheusTask(MyTask taskType)
        {
            this.TaskType = taskType;
        }

        #endregion Public Constructors

        #region Public Events

        public static event EventHandler<SingleTaskEventArgs> FinishedSingleTaskHandler;

        public static event EventHandler<SingleFileEventArgs> FinishedWritingFileHandler;

        public static event EventHandler<SingleTaskEventArgs> StartingSingleTaskHander;

        public static event EventHandler<StringEventArgs> StartingDataFileHandler;

        public static event EventHandler<StringEventArgs> FinishedDataFileHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<StringEventArgs> NewCollectionHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        #endregion Public Events

        #region Public Properties

        public MyTask TaskType { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static SearchMode ParseSearchMode(string text)
        {
            SearchMode ye = null;

            var split = text.Split(' ');

            switch (split[1])
            {
                case "dot":
                    ToleranceUnit tu = ToleranceUnit.PPM;
                    if (split[3].ToUpperInvariant().Equals("PPM"))
                        tu = ToleranceUnit.PPM;
                    else if (split[3].ToUpperInvariant().Equals("DA"))
                        tu = ToleranceUnit.Absolute;
                    else
                        break;

                    var massShifts = Array.ConvertAll(split[4].Split(','), Double.Parse);
                    var newString = split[2].Replace("±", "");
                    var toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                    ye = new DotSearchMode(split[0], massShifts, new Tolerance(tu, toleranceValue));
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                    ye = new IntervalSearchMode(split[0], doubleRanges);
                    break;

                case "OpenSearch":
                    ye = new OpenSearchMode();
                    break;

                case "daltonsAroundZero":
                    ye = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                case "ppmAroundZero":
                    ye = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                default:
                    throw new Exception("Could not parse search mode string");
            }
            return ye;
        }

        public MyTaskResults RunTask(string output_folder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawDataFilenameList, string taskId)
        {
            StartingSingleTask(taskId);
            var paramsFileName = Path.Combine(output_folder, "prose.txt");
            using (StreamWriter file = new StreamWriter(paramsFileName))
            {
                file.WriteLine("MetaMorpheus version "
                    + (GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ? "NOT A RELEASE" : GlobalEngineLevelSettings.MetaMorpheusVersion)
                    + " is used to run a "
                    + this.TaskType
                    + " task on "
                    + currentRawDataFilenameList.Count
                    + " spectra files.");

                file.WriteLine(ToString());

                file.WriteLine();
                file.WriteLine("taskId: " + taskId);
                file.WriteLine("Spectra files:");
                file.WriteLine(string.Join(Environment.NewLine, currentRawDataFilenameList.Select(b => '\t' + b)));
                file.WriteLine("XML files:");
                file.Write(string.Join(Environment.NewLine, currentXmlDbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FileName)));
            }
            SucessfullyFinishedWritingFile(paramsFileName, new List<string> { taskId });

            // TOML
            var tomlFileName = Path.Combine(output_folder, GetType().Name + "config.toml");

            Toml.WriteFile(this, tomlFileName, tomlConfig);
            SucessfullyFinishedWritingFile(tomlFileName, new List<string> { taskId });

            MetaMorpheusEngine.FinishedSingleEngineHandler += SingleEngineHandlerInTask;

#if !DEBUG
            try
            {
#endif
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            RunSpecific(output_folder, currentXmlDbFilenameList, currentRawDataFilenameList, taskId);
            stopWatch.Stop();
            myTaskResults.Time = stopWatch.Elapsed;
            var resultsFileName = Path.Combine(output_folder, "results.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine(GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion);
                file.Write(myTaskResults.ToString());
            }
            SucessfullyFinishedWritingFile(resultsFileName, new List<string> { taskId });
            FinishedSingleTask(taskId);
#if !DEBUG
            }
            catch (Exception e)
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
                var resultsFileName = Path.Combine(output_folder, "results.txt");
                using (StreamWriter file = new StreamWriter(resultsFileName))
                {
                    file.WriteLine(GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion);
                    file.Write("e: " + e);
                    file.Write("e.Message: " + e.Message);
                    file.Write("e.InnerException: " + e.InnerException);
                    file.Write("e.Source: " + e.Source);
                    file.Write("e.StackTrace: " + e.StackTrace);
                    file.Write("e.TargetSite: " + e.TargetSite);
                }
                throw;
            }
#endif

            MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
            return myTaskResults;
        }

        #endregion Public Methods

        #region Protected Internal Methods

        protected internal void WritePsms(List<NewPsmWithFdr> items, List<DbForTask> dbFilenameList, string outputfolder, string filename, List<string> nestedIds)
        {
            WritePsmsToTsv(items, outputfolder, filename, nestedIds);
            WritePsmsToMzIdentmL(items, dbFilenameList, outputfolder, filename, nestedIds);
        }
        #endregion Protected Internal Methods

        #region Protected Methods

        protected void WritePsmsToTsv(List<NewPsmWithFdr> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".psmtsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(NewPsmWithFdr.TabSeparatedHeader);
                for (int i = 0; i < items.Count; i++)
                    output.WriteLine(items[i]);
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        protected internal void WritePsmsToMzIdentmL(List<NewPsmWithFdr> items, List<DbForTask> dbFilenameList, string outputFolder, string fileName, List<string> nestedIds)
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
                version = GlobalEngineLevelSettings.MetaMorpheusVersion, 
                uri = "https://github.com/smith-chem-wisc/MetaMorpheus"
            }};
            //sequence collection: database entries of protein/peptide sequences identified and modifications
            _mzid.SequenceCollection = new mzIdentML.Generated.SequenceCollectionType();
            _mzid.SequenceCollection.Peptide = new mzIdentML.Generated.PeptideType[items.Select(p => p.thisPSM.PeptidesWithSetModifications.First()).Distinct().Count()];
            _mzid.SequenceCollection.DBSequence = new mzIdentML.Generated.DBSequenceType[items.SelectMany(p => p.thisPSM.PeptidesWithSetModifications.Select(t => t.Protein)).Distinct().Count()];
            _mzid.SequenceCollection.PeptideEvidence = new mzIdentML.Generated.PeptideEvidenceType[items.SelectMany(p => p.thisPSM.PeptidesWithSetModifications).Distinct().Count()];
            _mzid.DataCollection = new mzIdentML.Generated.DataCollectionType();
            _mzid.DataCollection.AnalysisData = new mzIdentML.Generated.AnalysisDataType()
            {
                SpectrumIdentificationList = new mzIdentML.Generated.SpectrumIdentificationListType[1]
            { new mzIdentML.Generated.SpectrumIdentificationListType()
            { id = "SIL_1",
             SpectrumIdentificationResult = new mzIdentML.Generated.SpectrumIdentificationResultType[items.Select(p => p.thisPSM.newPsm.scanNumber).Distinct().Count()],},
            },
                //include protein groups?? 
                //ProteinDetectionList = new mzIdentML.Generated.ProteinDetectionListType[1]
                //{
                //    new mzIdentML.Generated.ProteinDetectionListType()
                //    {
                //        id = "PDL_1",
                //        ProteinAmbiguityGroup =
                //    }
                //}
            };

            //DataColection
            _mzid.DataCollection.Inputs = new mzIdentML.Generated.InputsType()
            {
                SpectraData = new mzIdentML.Generated.SpectraDataType[items.Select(p => p.thisPSM.newPsm.fileName).Distinct().Count()],
                SearchDatabase = new mzIdentML.Generated.SearchDatabaseType[dbFilenameList.Count] //ACCESS # OF DATABASES USED?
            };

            int protein_index = 0; // protein overall
            int peptide_index = 0; // peptide overall
            int peptide_envidence = 0; // peptide evidence overall
            int mod_count = 0; // mod count for each peptide
            int file_index = 0; //file index

            Dictionary<Protein, string> databaseRef = new Dictionary<Protein, string>(); //key protein, value database ID
            Dictionary<PeptideWithSetModifications, int> peptideRef = new Dictionary<PeptideWithSetModifications, int>(); //key peptide, value = peptide reference 
            Dictionary<PeptideWithSetModifications, string> peptideEvidenceRef = new Dictionary<PeptideWithSetModifications, string>(); //key peptide evidence, value = peptide evidence reference
            Dictionary<string, int> fileRef = new Dictionary<string, int>(); //key filename, value = file index
            Dictionary<string, List<int>> scanRef = new Dictionary<string, List<int>>(); //key = filename, value = list of scan ints 
            Dictionary<string, List<double>> spectrumIdentificationItemRef = new Dictionary<string, List<double>>(); //key = filename_scan, value = list of spectrum identification items for that scan

            //Analysis collection: the analyses performed to get the results, which map the input and output data sets. 
            _mzid.AnalysisCollection = new mzIdentML.Generated.AnalysisCollectionType()
            {
                SpectrumIdentification = new mzIdentML.Generated.SpectrumIdentificationType[1] { new mzIdentML.Generated.SpectrumIdentificationType()
                {
                    id = "SI",
                    spectrumIdentificationProtocol_ref = "SIP",
                    spectrumIdentificationList_ref = "SIL_1",
                    activityDate = DateTime.Now,
                    SearchDatabaseRef = new mzIdentML.Generated.SearchDatabaseRefType[1] { new mzIdentML.Generated.SearchDatabaseRefType() {searchDatabase_ref = "Uniprot" } },
                    InputSpectra = new mzIdentML.Generated.InputSpectraType[items.Select(p => p.thisPSM.newPsm.fileName).Distinct().Count()],
                } },
                ProteinDetection = new mzIdentML.Generated.ProteinDetectionType()
                {
                    id = "PD_1",
                    proteinDetectionProtocol_ref = "PDP_1",
                    proteinDetectionList_ref = "PDL_1",
                    InputSpectrumIdentifications = new mzIdentML.Generated.InputSpectrumIdentificationsType[1] {new mzIdentML.Generated.InputSpectrumIdentificationsType()
                    { spectrumIdentificationList_ref = "SIL_1" }
                }
                }
            };
            foreach (string filename in items.Select(p => p.thisPSM.newPsm.fileName).Distinct().OrderBy(p => p).ToList())
            {

                int spectrum_id_index = 0;
                _mzid.AnalysisCollection.SpectrumIdentification[0].InputSpectra[file_index] = new mzIdentML.Generated.InputSpectraType()
                {
                    spectraData_ref = "SD_" + file_index
                };
                _mzid.DataCollection.Inputs.SpectraData[file_index] = new mzIdentML.Generated.SpectraDataType()
                {
                    id = "SD_" + file_index,
                    location = "file:///" + filename + ".mzml",
                    FileFormat = new mzIdentML.Generated.FileFormatType()
                    {
                        cvParam = new mzIdentML.Generated.CVParamType()
                        {
                            accession = "MS:1000584",
                            cvRef = "PSI-MS",
                            name = "mzML format"
                        }
                    }
                };
                foreach (int scan in items.Where(p => p.thisPSM.newPsm.fileName == filename).Select(p => p.thisPSM.newPsm.scanNumber).Distinct())
                {
                    //spectrum identification result: all identifications made from searching one spectrum (could be more than one if DIA or coisolation)
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[spectrum_id_index] = new mzIdentML.Generated.SpectrumIdentificationResultType()
                    {
                        spectrumID = ("scan=" + scan).ToString(),
                        SpectrumIdentificationItem = new mzIdentML.Generated.SpectrumIdentificationItemType[items.Where(p => p.thisPSM.newPsm.fileName == filename && p.thisPSM.newPsm.scanNumber == scan).Count()]
                    };
                    spectrum_id_index++;
                    spectrumIdentificationItemRef.Add(filename + "_" + scan, new List<double>());
                }
                fileRef.Add(filename, file_index);
                scanRef.Add(filename, new List<int>());
                file_index++;
            }

            foreach (NewPsmWithFdr psm in items)
            {
                //new peptide object per psm (all peptides in peptidesWithSetModifications have same sequence and modifications, only need one)
                if (!peptideRef.ContainsKey(psm.thisPSM.PeptidesWithSetModifications.First()))
                {
                    PeptideWithSetModifications peptide = psm.thisPSM.PeptidesWithSetModifications.First();
                    string peptide_id = "peptide_" + peptide_index;
                    peptideRef.Add(peptide, peptide_index);
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
                            ModificationWithMass mod = peptide.allModsOneIsNterminus[position];
                            _mzid.SequenceCollection.Peptide[peptide_index].Modification[mod_count] = new mzIdentML.Generated.ModificationType()
                            {
                                locationSpecified = true,
                                location = position,
                                monoisotopicMassDeltaSpecified = true,
                                monoisotopicMassDelta = mod.monoisotopicMass,
                                cvParam = new mzIdentML.Generated.CVParamType[1] { new mzIdentML.Generated.CVParamType()
                            {
                                cvRef = "MetaMorpheus",
                                name = "MetaMorpheusModification",
                                accession = mod.id
                            } }
                            };
                            mod_count++;
                        }
                    }
                    peptide_index++;
                }

                foreach (PeptideWithSetModifications peptide_evidence in psm.thisPSM.PeptidesWithSetModifications)
                {
                    //if DB references doesn't contain protein, add it
                    Protein protein = peptide_evidence.Protein;
                    if (!databaseRef.ContainsKey(protein))
                    {
                        string database_sequence_ref = "DBSeq_" + protein.Accession;
                        databaseRef.Add(protein, database_sequence_ref);
                        _mzid.SequenceCollection.DBSequence[protein_index] = new mzIdentML.Generated.DBSequenceType()
                        {
                            id = database_sequence_ref,
                            length = protein.Length,
                            searchDatabase_ref = "Uniprot", //?
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
                        protein_index++;
                    }

                    //peptide evdience for each peptide in peptidesWithSetModifications
                    if (!peptideEvidenceRef.ContainsKey(peptide_evidence))
                    {
                        string id = "PE_" + peptideRef[psm.thisPSM.PeptidesWithSetModifications.First()] + "_" + peptide_envidence;
                        _mzid.SequenceCollection.PeptideEvidence[peptide_envidence] = new mzIdentML.Generated.PeptideEvidenceType()
                        {
                            isDecoy = peptide_evidence.Protein.IsDecoy,
                            peptide_ref = "peptide_" + peptideRef[psm.thisPSM.PeptidesWithSetModifications.First()],
                            id = id,
                            startSpecified = true,
                            endSpecified = true,
                            start = peptide_evidence.OneBasedStartResidueInProtein,
                            end = peptide_evidence.OneBasedEndResidueInProtein,
                            pre = peptide_evidence.PreviousAminoAcid.ToString(),
                            post = (peptide_evidence.OneBasedEndResidueInProtein < protein.BaseSequence.Length) ? protein[peptide_evidence.OneBasedEndResidueInProtein].ToString() : "-",
                            dBSequence_ref = databaseRef[peptide_evidence.Protein]
                        };
                        peptideEvidenceRef.Add(peptide_evidence, id);
                        peptide_envidence++;
                    }
                }

                if (!scanRef[psm.thisPSM.newPsm.fileName].Contains(psm.thisPSM.newPsm.scanNumber))
                {
                    scanRef[psm.thisPSM.newPsm.fileName].Add(psm.thisPSM.newPsm.scanNumber);
                }
                int scan_ref = scanRef[psm.thisPSM.newPsm.fileName].IndexOf(psm.thisPSM.newPsm.scanNumber);
                spectrumIdentificationItemRef[psm.thisPSM.newPsm.fileName + "_" + psm.thisPSM.newPsm.scanNumber].Add(psm.thisPSM.newPsm.scanPrecursorMZ);
                int spectrum_id_result = spectrumIdentificationItemRef[psm.thisPSM.newPsm.fileName + "_" + psm.thisPSM.newPsm.scanNumber].IndexOf(psm.thisPSM.newPsm.scanPrecursorMZ);
                _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scanRef[psm.thisPSM.newPsm.fileName].IndexOf(psm.thisPSM.newPsm.scanNumber)].SpectrumIdentificationItem[spectrum_id_result] = new mzIdentML.Generated.SpectrumIdentificationItemType()
                {
                    experimentalMassToCharge = psm.thisPSM.newPsm.scanPrecursorMZ,
                    calculatedMassToCharge = psm.thisPSM.PeptideMonoisotopicMass.ToMz(psm.thisPSM.newPsm.scanPrecursorCharge),
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
                                    accession = "MS:1002354",
                                    name = "PSM-level q-value",
                                    cvRef = "PSI-MS",
                                    value = psm.QValue.ToString()
                                }
                         },
                    PeptideEvidenceRef = new mzIdentML.Generated.PeptideEvidenceRefType[1] { new mzIdentML.Generated.PeptideEvidenceRefType()
                        {
                            peptideEvidence_ref = peptideEvidenceRef[psm.thisPSM.PeptidesWithSetModifications.First()]
                        }
                        },
                    Fragmentation = new mzIdentML.Generated.IonTypeType[psm.thisPSM.newPsm.matchedIonsListPositiveIsMatch.Count]
                };
                int ion_list_count = 0;
                foreach (var ion_list in psm.thisPSM.newPsm.matchedIonsListPositiveIsMatch)
                {
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scanRef[psm.thisPSM.newPsm.fileName].IndexOf(psm.thisPSM.newPsm.scanNumber)].SpectrumIdentificationItem[spectrum_id_result].Fragmentation[ion_list_count] = new mzIdentML.Generated.IonTypeType()
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
                        _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_ref].SpectrumIdentificationItem[spectrum_id_result].Fragmentation[ion_list_count].FragmentArray[0].values[i] = Math.Abs((float)ion_list.Value[i]);
                    }
                    ion_list_count++;
                }
            }

            //Analysis protocol collection: parameters
            _mzid.AnalysisProtocolCollection = new mzIdentML.Generated.AnalysisProtocolCollectionType();
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol = new mzIdentML.Generated.SpectrumIdentificationProtocolType[1] { new mzIdentML.Generated.SpectrumIdentificationProtocolType()
            {
                id = "SIP", analysisSoftware_ref = "MetaMorpheus"
            } };
            //add later? - already in params.txt file, may want here also

            TextWriter writer = new StreamWriter(Path.Combine(outputFolder, fileName + ".mzid"));
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
            SucessfullyFinishedWritingFile(Path.Combine(outputFolder, fileName + ".mzid"), nestedIds);
        }



        protected static List<Protein> LoadProteinDb(string fileName, bool generateDecoys, List<ModificationWithMass> localizeableModifications, bool isContaminant, IEnumerable<string> dbRefTypesToKeep, out Dictionary<string, Modification> um)
        {
            if (Path.GetExtension(fileName).Equals(".fasta"))
            {
                um = null;
                return ProteinDbLoader.LoadProteinFasta(fileName, generateDecoys, isContaminant, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_gene_expression);
            }
            else
                return ProteinDbLoader.LoadProteinXML(fileName, generateDecoys, localizeableModifications, isContaminant, null, out um);
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyTaskResults RunSpecific(string output_folder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawDataFilenameList, string taskId);

        protected void WriteProteinGroupsToTsv(List<ProteinGroup> items, string outputFolder, string fileName, List<string> nestedIds)
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

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        protected void WriteTree(BinTreeStructure myTreeStructure, string output_folder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(output_folder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
                foreach (Bin bin in myTreeStructure.FinalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F4", CultureInfo.InvariantCulture)
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
                        + "\t" + bin.UnimodDiffs
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.protNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.pepNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.pepClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.protClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.MedianLength).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        protected void SucessfullyFinishedWritingFile(string path, List<string> nestedIDs)
        {
            FinishedWritingFileHandler?.Invoke(this, new SingleFileEventArgs(path, nestedIDs));
        }

        protected void StartingDataFile(string v, List<string> nestedIDs)
        {
            StartingDataFileHandler?.Invoke(this, new StringEventArgs(v, nestedIDs));
        }

        protected void FinishedDataFile(string v, List<string> nestedIDs)
        {
            FinishedDataFileHandler?.Invoke(this, new StringEventArgs(v, nestedIDs));
        }

        protected void Status(string v, List<string> nestedIds)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void Warn(string v, List<string> nestedIds)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void NewCollection(string v, List<string> nestedIds)
        {
            NewCollectionHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        #endregion Protected Methods

        #region Private Methods

        private void SingleEngineHandlerInTask(object sender, SingleEngineFinishedEventArgs e)
        {
            myTaskResults.AddResultText(e.ToString());
        }

        private void FinishedSingleTask(string taskId)
        {
            FinishedSingleTaskHandler?.Invoke(this, new SingleTaskEventArgs(taskId));
        }

        private void StartingSingleTask(string taskId)
        {
            StartingSingleTaskHander?.Invoke(this, new SingleTaskEventArgs(taskId));
        }

        #endregion Private Methods

    }
}