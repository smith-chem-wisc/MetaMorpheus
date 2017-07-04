using Chemistry;
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

        public static readonly TomlSettings tomlConfig = TomlSettings.Create(cfg => cfg
                        .ConfigureType<Tolerance>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => new Tolerance(tmlString.Value))))
                        .ConfigureType<MassDiffAcceptor>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => MetaMorpheusTask.ParseSearchMode(tmlString.Value))))
                        .ConfigureType<Protease>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => GlobalTaskLevelSettings.ProteaseDictionary[tmlString.Value])))
                        .ConfigureType<List<Tuple<string, string>>>(type => type
                             .WithConversionFor<TomlString>(convert => convert
                                 .ToToml(custom => string.Join("\t\t", custom.Select(b => b.Item1 + "\t" + b.Item2)))
                                 .FromToml(tmlString => GetModsFromString(tmlString.Value)))));

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

        public int? MaxDegreeOfParallelism { get; set; }
        public bool LocalizeAll { get; set; }
        public List<Tuple<string, string>> ListOfModsFixed { get; set; }
        public List<Tuple<string, string>> ListOfModsVariable { get; set; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; set; }
        public MyTask TaskType { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor ye = null;

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
                    ye = new DotMassDiffAcceptor(split[0], massShifts, new Tolerance(tu, toleranceValue));
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                    ye = new IntervalMassDiffAcceptor(split[0], doubleRanges);
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

            #region Write Prose

            {
                var proseFilePath = Path.Combine(output_folder, "prose.txt");
                using (StreamWriter file = new StreamWriter(proseFilePath))
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
                SucessfullyFinishedWritingFile(proseFilePath, new List<string> { taskId });
            }

            #endregion Write Prose

            #region write TOML

            {
                var tomlFileName = Path.Combine(output_folder, GetType().Name + "config.toml");
                Toml.WriteFile(this, tomlFileName, tomlConfig);
                SucessfullyFinishedWritingFile(tomlFileName, new List<string> { taskId });
            }

            #endregion write TOML

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

        protected internal void WritePsmsToTsv(IEnumerable<PsmParent> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".psmtsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(PsmParent.GetTabSeparatedHeader());
                foreach (var heh in items)
                    output.WriteLine(heh);
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        protected internal void WriteMzidentml(IEnumerable<PsmParent> items, List<ProteinGroup> groups, List<ModificationWithMass> variableMods, List<ModificationWithMass> fixedMods, List<Protease> proteases, double threshold, MassDiffAcceptor searchMode, Tolerance productTolerance, int missedCleavages, string outputFolder, string fileName, List<string> nestedIds)
        {
            List<PeptideWithSetModifications> peptides = items.SelectMany(i => i.Pli.PeptidesWithSetModifications).Distinct().ToList();
            List<Protein> proteins = peptides.Select(p => p.Protein).Distinct().ToList();
            List<string> filenames = items.Select(i => i.FullFilePath).Distinct().ToList();

            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML.Generated.MzIdentMLType));
            var _mzid = new mzIdentML.Generated.MzIdentMLType()
            {
                version = "1.1.1",
            };

            //cvlist: URLs of controlled vocabularies used within the file.
            _mzid.cvList = new mzIdentML.Generated.cvType[2] { new mzIdentML.Generated.cvType()
            {
                id = "PSI-MS",
                fullName = "Proteomics Standards Initiative Mass Spectrometry Vocabularies",
                uri = "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",
                version= "4.0.9"
            },
            new mzIdentML.Generated.cvType()
            {
                id = "UO",
                fullName = "UNIT-ONTOLOGY",
                uri = "http://www.unimod.org/obo/unimod.obo"
            }
            };

            _mzid.AnalysisSoftwareList = new mzIdentML.Generated.AnalysisSoftwareType[1] {  new mzIdentML.Generated.AnalysisSoftwareType()
            {
                id = "AS_MetaMorpheus",
                name = "MetaMorpheus",
                version = GlobalEngineLevelSettings.MetaMorpheusVersion,
                uri = "https://github.com/smith-chem-wisc/MetaMorpheus",
                SoftwareName = new mzIdentML.Generated.ParamType()
                {
                    Item = new mzIdentML.Generated.CVParamType()
                    {
                        //using Morpheus's accession until we get MetaMorpheus entered
                        accession = "MS:1002661",
                        name = "Morpheus",
                        cvRef = "PSI-MS"
                    }
                }
            }};
            _mzid.SequenceCollection = new mzIdentML.Generated.SequenceCollectionType()
            {
                Peptide = new mzIdentML.Generated.PeptideType[peptides.Count()],
                DBSequence = new mzIdentML.Generated.DBSequenceType[proteins.Count()],
                PeptideEvidence = new mzIdentML.Generated.PeptideEvidenceType[peptides.Count()]
            };
            int protein_index = 0;
            foreach (Protein protein in proteins)
            {
                _mzid.SequenceCollection.DBSequence[protein_index] = new mzIdentML.Generated.DBSequenceType()
                {
                    id = "DBS_" + protein.Accession,
                    lengthSpecified = true,
                    length = protein.Length,
                    searchDatabase_ref = "SDB_1", //TODO: SPECIFIC DATABASE THE PROTEIN CAME FROM? NULL FOR FOR PROTEINS
                    accession = protein.Accession,
                    Seq = protein.BaseSequence,
                    cvParam = new mzIdentML.Generated.CVParamType[1]
                    {
                        new mzIdentML.Generated.CVParamType()
                        {
                            accession = "MS:1001088",
                            name = "protein description",
                            cvRef = "PSI-MS",
                            value = protein.FullDescription
                        }
                    },
                    name = protein.Name
                };
                protein_index++;
            }

            _mzid.DataCollection = new mzIdentML.Generated.DataCollectionType()
            {
                AnalysisData = new mzIdentML.Generated.AnalysisDataType()
                {
                    SpectrumIdentificationList = new mzIdentML.Generated.SpectrumIdentificationListType[1]
                    {
                        new mzIdentML.Generated.SpectrumIdentificationListType()
                        {
                            id = "SIL",
                            SpectrumIdentificationResult = new mzIdentML.Generated.SpectrumIdentificationResultType[items.Count()]
                        }
                    }
                },
                Inputs = new mzIdentML.Generated.InputsType()
                {
                    //TODO: SEPARATE SEARCH DATABASES IF MULTIPLE DATABASES ENTERED
                    SearchDatabase = new mzIdentML.Generated.SearchDatabaseType[1]
                    {
                        new mzIdentML.Generated.SearchDatabaseType()
                        {
                         id = "SDB_1",
                            FileFormat = new mzIdentML.Generated.FileFormatType()
                            {
                                cvParam = new mzIdentML.Generated.CVParamType()
                                {
                                    accession = "MS:1002660",
                                    name = "UniProtKB XML sequence format",
                                    cvRef = "PSI-MS"
                                }
                            },
                            DatabaseName = new mzIdentML.Generated.ParamType()
                            {
                                Item = new mzIdentML.Generated.CVParamType()
                                {
                                    accession = "MS:1001073",
                                    name = "database type amino acid",
                                    cvRef = "PSI-MS"
                                }
                            },
                        }
                    },
                    SpectraData = new mzIdentML.Generated.SpectraDataType[filenames.Count]
                }
            };

            _mzid.AnalysisCollection = new mzIdentML.Generated.AnalysisCollectionType()
            {
                SpectrumIdentification = new mzIdentML.Generated.SpectrumIdentificationType[1]
                {
                    new mzIdentML.Generated.SpectrumIdentificationType()
                    {
                        id = "SI",
                        spectrumIdentificationList_ref = "SIL",
                        spectrumIdentificationProtocol_ref = "SIP",
                        InputSpectra = new mzIdentML.Generated.InputSpectraType[filenames.Count],
                        SearchDatabaseRef = new mzIdentML.Generated.SearchDatabaseRefType[proteins.Count]
                    }
                }
            };

            Dictionary<string, int> spectral_ids = new Dictionary<string, int>(); //key is datafile, value is datafile's id
            int spectra_data_id = 0;
            foreach (string data_filepath in filenames)
            {
                bool thermoRawFile = Path.GetExtension(data_filepath) == ".raw";
                string spectral_data_id = "SD_" + spectra_data_id;
                spectral_ids.Add(data_filepath, spectra_data_id);
                _mzid.AnalysisCollection.SpectrumIdentification[0].InputSpectra[spectra_data_id] = new mzIdentML.Generated.InputSpectraType()
                {
                    spectraData_ref = spectral_data_id
                };
                _mzid.DataCollection.Inputs.SpectraData[spectra_data_id] = new mzIdentML.Generated.SpectraDataType()
                {
                    id = spectral_data_id,
                    name = Path.GetFileNameWithoutExtension(data_filepath),
                    location = Path.GetPathRoot(data_filepath),
                    FileFormat = new mzIdentML.Generated.FileFormatType()
                    {
                        cvParam = new mzIdentML.Generated.CVParamType()
                        {
                            accession = thermoRawFile ? "MS:1000563" : "MS:1000584",
                            name = thermoRawFile ? "Thermo RAW format" : "mzML format",
                            cvRef = "PSI-MS"
                        }
                    },
                    SpectrumIDFormat = new mzIdentML.Generated.SpectrumIDFormatType()
                    {
                        cvParam = new mzIdentML.Generated.CVParamType()
                        {
                            accession = thermoRawFile ? "MS:1000768" : "MS:1001530",
                            name = thermoRawFile ? "Thermo nativeID format" : "mzML unique identifier",
                            cvRef = "PSI-MS"
                        }
                    }
                };
                spectra_data_id++;
            }

            int sir_id = 0;
            int pe_index = 0;
            int p_index = 0;
            Dictionary<PeptideWithSetModifications, Tuple<int, int, List<string>>> peptide_ids = new Dictionary<PeptideWithSetModifications, Tuple<int, int, List<string>>>(); //key is peptide, value is <peptide id for that peptide, peptide evidence id>
            Dictionary<Tuple<string, int>, Tuple<int, int>> psm_per_scan = new Dictionary<Tuple<string, int>, Tuple<int, int>>(); //key is <filename, scan numer> value is <scan result id, scan item id #'s (could be more than one ID per scan)>
            foreach (PsmParent psm in items)
            {
                PeptideWithSetModifications peptide = psm.Pli.PeptidesWithSetModifications.OrderBy(p => p.PeptideDescription).First();
                Tuple<int, int, List<string>> peptide_id;
                //if first peptide on list hasn't been added, add peptide and peptide evidence
                if (!peptide_ids.TryGetValue(peptide, out peptide_id))
                {
                    peptide_id = new Tuple<int, int, List<string>>(p_index, 0, new List<string>());
                    p_index++;
                    _mzid.SequenceCollection.Peptide[peptide_id.Item1] = new mzIdentML.Generated.PeptideType()
                    {
                        PeptideSequence = peptide.BaseSequence,
                        id = "P_" + peptide_id.Item1,
                        Modification = new mzIdentML.Generated.ModificationType[peptide.NumMods]
                    };
                    int mod_id = 0;
                    foreach (KeyValuePair<int, ModificationWithMass> mod in peptide.allModsOneIsNterminus)
                    {
                        _mzid.SequenceCollection.Peptide[peptide_id.Item1].Modification[mod_id] = new mzIdentML.Generated.ModificationType()
                        {
                            location = mod.Key - 1,
                            locationSpecified = true,
                            monoisotopicMassDelta = mod.Value.monoisotopicMass,
                            residues = new string[1] { mod.Value.motif.Motif.ToString() },
                            monoisotopicMassDeltaSpecified = true,
                            cvParam = new mzIdentML.Generated.CVParamType[1]
                            {
                            new mzIdentML.Generated.CVParamType()
                            {
                                cvRef = "PSI-MS",
                                name = "unknown modification",
                                accession = "MS:1001460",
                                value = mod.Value.id
                            }
                            }
                        };
                        mod_id++;
                    }

                    foreach (PeptideWithSetModifications peptide_evidence in psm.Pli.PeptidesWithSetModifications)
                    {
                        _mzid.SequenceCollection.PeptideEvidence[pe_index] = new mzIdentML.Generated.PeptideEvidenceType()
                        {
                            id = "PE_" + pe_index,
                            peptide_ref = "P_" + peptide_id.Item1,
                            dBSequence_ref = "DBS_" + peptide_evidence.Protein.Accession,
                            isDecoy = peptide_evidence.Protein.IsDecoy,
                            startSpecified = true,
                            start = peptide.OneBasedStartResidueInProtein,
                            endSpecified = true,
                            end = peptide.OneBasedEndResidueInProtein,
                            pre = peptide.PreviousAminoAcid.ToString(),
                            post = (peptide.OneBasedEndResidueInProtein < peptide.Protein.BaseSequence.Length) ? peptide.Protein[peptide.OneBasedEndResidueInProtein].ToString() : "-",
                        };
                        peptide_ids.Add(peptide_evidence, new Tuple<int, int, List<string>>(peptide_id.Item1, pe_index, new List<string>()));
                        pe_index++;
                    }
                }

                Tuple<int, int> scan_result_scan_item;
                if (!psm_per_scan.TryGetValue(new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber), out scan_result_scan_item)) //check to see if scan has already been added
                {
                    scan_result_scan_item = new Tuple<int, int>(sir_id, 0);
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1] = new mzIdentML.Generated.SpectrumIdentificationResultType()
                    {
                        id = "SIR_" + scan_result_scan_item.Item1,
                        spectraData_ref = "SD_" + spectral_ids[psm.FullFilePath].ToString(),
                        spectrumID = psm.ScanNumber.ToString(),
                        SpectrumIdentificationItem = new mzIdentML.Generated.SpectrumIdentificationItemType[500]
                    };
                    psm_per_scan.Add(new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber), scan_result_scan_item);
                    sir_id++;
                }
                else
                {
                    psm_per_scan[new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber)] = new Tuple<int, int>(scan_result_scan_item.Item1, scan_result_scan_item.Item2 + 1);
                    scan_result_scan_item = psm_per_scan[new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber)];
                }
                foreach (PeptideWithSetModifications p in psm.Pli.PeptidesWithSetModifications)
                {
                    peptide_ids[p].Item3.Add("SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2);
                }
                _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2] = new mzIdentML.Generated.SpectrumIdentificationItemType()
                {
                    chargeState = psm.ScanPrecursorCharge,
                    id = "SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2,
                    experimentalMassToCharge = psm.ScanPrecursorMonoisotopicPeak.Mz,
                    calculatedMassToCharge = psm.Pli.PeptideMonoisotopicMass.ToMz(psm.ScanPrecursorCharge),
                    calculatedMassToChargeSpecified = true,
                    passThreshold = psm.FdrInfo.QValue <= threshold,
                    rank = 1,
                    peptide_ref = "P_" + peptide_id.Item1,
                    PeptideEvidenceRef = new mzIdentML.Generated.PeptideEvidenceRefType[psm.Pli.PeptidesWithSetModifications.Count],
                    cvParam = new mzIdentML.Generated.CVParamType[2]
                    {
                        new mzIdentML.Generated.CVParamType()
                        {
                            name = "Morpheus:Morpheus score",
                            cvRef = "PSI-MS",
                            accession = "MS:1002662",
                            value = psm.Score.ToString()
                        },
                        new mzIdentML.Generated.CVParamType()
                        {
                            accession = "MS:1002354",
                            name = "PSM-level q-value",
                            cvRef = "PSI-MS",
                            value = psm.FdrInfo.QValue.ToString()
                        }
                    }
                };

                int pe = 0;
                foreach (PeptideWithSetModifications p in psm.Pli.PeptidesWithSetModifications)
                {
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2].PeptideEvidenceRef[pe]
                        = new mzIdentML.Generated.PeptideEvidenceRefType()
                        {
                            peptideEvidence_ref = "PE_" + peptide_ids[p].Item2
                        };
                    pe++;
                }
            }

            _mzid.AnalysisProtocolCollection = new mzIdentML.Generated.AnalysisProtocolCollectionType()
            {
                SpectrumIdentificationProtocol = new mzIdentML.Generated.SpectrumIdentificationProtocolType[1]
                {
                    new mzIdentML.Generated.SpectrumIdentificationProtocolType()
                    {
                        id = "SIP",
                        analysisSoftware_ref = "AS_MetaMorpheus",
                        SearchType = new mzIdentML.Generated.ParamType()
                        {
                            Item = new mzIdentML.Generated.CVParamType()
                            {
                                accession = "MS:1001083",
                                name = "ms-ms search",
                                cvRef = "PSI-MS"
                            }
                        },
                        AdditionalSearchParams = new mzIdentML.Generated.ParamListType()
                        {
                            //TODO: ADD SEARCH PARAMS?
                            Items = new mzIdentML.Generated.AbstractParamType[2]
                            {
                                new mzIdentML.Generated.CVParamType()
                                {
                                    accession = "MS:1001211",
                                    cvRef = "PSI-MS",
                                    name = "parent mass type mono"
                                },
                                new mzIdentML.Generated.CVParamType()
                                {
                                    accession = "MS:1001255",
                                    name = "fragment mass type mono",
                                    cvRef = "PSI-MS"
                                },
                            }
                        },
                        ModificationParams = new mzIdentML.Generated.SearchModificationType[fixedMods.Count() + variableMods.Count()],
                        Enzymes = new mzIdentML.Generated.EnzymesType()
                        {
                            Enzyme = new mzIdentML.Generated.EnzymeType[proteases.Count]
                        },
                        FragmentTolerance = new mzIdentML.Generated.CVParamType[2]
                        {
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession = "MS:1001412",
                                name = "search tolerance plus value",
                                value = productTolerance.Value.ToString(),
                                cvRef = "PSI-MS",

                                unitCvRef = "UO"
                            },
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession = "MS:1001413",
                                name = "search tolerance minus value",
                                value = productTolerance.Value.ToString(),
                                cvRef = "PSI-MS",
                                unitAccession = productTolerance.Unit == ToleranceUnit.PPM? "UO:0000169": "UO:0000221",
                                unitName = productTolerance.Unit == ToleranceUnit.PPM? "parts per million" : "dalton" ,
                                unitCvRef = "UO"
                            }
                        },
                        ParentTolerance = new mzIdentML.Generated.CVParamType[1]
                        {
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession = "MS1001411",
                                name = "search tolerance specification",
                                cvRef = "PSI-MS",
                                value = searchMode.FileNameAddition
                            }
                        },
                        Threshold = new mzIdentML.Generated.ParamListType()
                        {
                            Items = new mzIdentML.Generated.CVParamType[1]
                            {
                                new mzIdentML.Generated.CVParamType()
                                {
                                    accession = "MS:1001448",
                                    name = "pep:FDR threshold",
                                    cvRef = "PSI-MS",
                                    value = threshold.ToString()
                                }
                            }
                        }
                    }
                }
            };

            int protease_index = 0;
            foreach (Protease protease in proteases)
            {
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].Enzymes.Enzyme[protease_index] = new mzIdentML.Generated.EnzymeType()
                {
                    id = "E_" + protease_index,
                    name = protease.Name,
                    semiSpecific = protease.CleavageSpecificity == CleavageSpecificity.Semi,
                    missedCleavages = missedCleavages,
                    SiteRegexp = protease.SiteRegexp,
                    EnzymeName = new mzIdentML.Generated.ParamListType()
                    {
                        Items = new mzIdentML.Generated.AbstractParamType[1]
                        {
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession = protease.PsiMsAccessionNumber,
                                name = protease.PsiMsName,
                                cvRef = "PSI-MS"
                            }
                        }
                    }
                };
                protease_index++;
            }

            int mod_index = 0;
            foreach (ModificationWithMass mod in fixedMods)
            {
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = new mzIdentML.Generated.SearchModificationType()
                {
                    fixedMod = true,
                    massDelta = mod.monoisotopicMass,
                    residues = mod.motif.Motif,
                };
                mod_index++;
            }
            foreach (ModificationWithMass mod in variableMods)
            {
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = new mzIdentML.Generated.SearchModificationType()
                {
                    fixedMod = false,
                    massDelta = mod.monoisotopicMass,
                    residues = mod.motif.Motif,
                };
                mod_index++;
            }

            _mzid.AnalysisProtocolCollection.ProteinDetectionProtocol = new mzIdentML.Generated.ProteinDetectionProtocolType()
            {
                id = "PDP",
                analysisSoftware_ref = "AS_MetaMorpheus",
                Threshold = new mzIdentML.Generated.ParamListType()
                {
                    Items = new mzIdentML.Generated.CVParamType[1]
                    {
                        new mzIdentML.Generated.CVParamType()
                        {
                            accession = "MS:1001448",
                            name = "pep:FDR threshold",
                            cvRef = "PSI-MS",
                            value = threshold.ToString()
                        }
                    }
                }
            };

            //TODO: SEARCH DATABASE FOR EACH PROTEIN DATABASE INPUT
            _mzid.AnalysisCollection.SpectrumIdentification[0].SearchDatabaseRef[0] = new mzIdentML.Generated.SearchDatabaseRefType()
            {
                searchDatabase_ref = "SDB_1"
            };

            if (groups != null)
            {
                _mzid.DataCollection.AnalysisData.ProteinDetectionList = new mzIdentML.Generated.ProteinDetectionListType()
                {
                    id = "PDL",
                    ProteinAmbiguityGroup = new mzIdentML.Generated.ProteinAmbiguityGroupType[groups.Count]
                };

                int group_id = 0;
                foreach (ProteinGroup proteinGroup in groups)
                {
                    _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id] = new mzIdentML.Generated.ProteinAmbiguityGroupType()
                    {
                        id = "PAG_" + group_id,
                        ProteinDetectionHypothesis = new mzIdentML.Generated.ProteinDetectionHypothesisType[proteinGroup.Proteins.Count]
                    };
                    int protein_id = 0;
                    foreach (Protein protein in proteinGroup.Proteins)
                    {
                        _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[protein_id] = new mzIdentML.Generated.ProteinDetectionHypothesisType()
                        {
                            id = "PDH_" + protein_id,
                            dBSequence_ref = "DBS_" + protein.Accession,
                            passThreshold = proteinGroup.QValue <= threshold,
                            PeptideHypothesis = new mzIdentML.Generated.PeptideHypothesisType[proteinGroup.AllPeptides.Count],
                            cvParam = new mzIdentML.Generated.CVParamType[4]
                            {
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession = "MS:1002663",
                                name = "Morpheus:summed Morpheus score",
                                cvRef = "PSI-MS",
                                value = proteinGroup.ProteinGroupScore.ToString()
                            },
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession = "MS1002373",
                                name = "protein group-level q-value",
                                value = proteinGroup.QValue.ToString()
                            },
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession =  "MS:1001093",
                                name = "sequence coverage",
                                cvRef = "PSI-MS",
                                value = proteinGroup.SequenceCoveragePercent.ToString()
                            },
                            new mzIdentML.Generated.CVParamType()
                            {
                                accession = "MS:1001097",
                                name = "distinct peptide sequences",
                                cvRef = "PSI-MS",
                                value = proteinGroup.UniquePeptides.Count.ToString()
                            }
                            }
                        };
                        int peptide_id = 0;
                        foreach (PeptideWithSetModifications peptide in proteinGroup.AllPeptides)
                        {
                            _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[protein_id].PeptideHypothesis[peptide_id] = new mzIdentML.Generated.PeptideHypothesisType()
                            {
                                peptideEvidence_ref = "PE_" + peptide_ids[peptide].Item2,
                                SpectrumIdentificationItemRef = new mzIdentML.Generated.SpectrumIdentificationItemRefType[peptide_ids[peptide].Item3.Count]
                            };
                            int i = 0;
                            foreach (string sii in peptide_ids[peptide].Item3)
                            {
                                _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[protein_id].PeptideHypothesis[peptide_id].SpectrumIdentificationItemRef[i] = new mzIdentML.Generated.SpectrumIdentificationItemRefType()
                                {
                                    spectrumIdentificationItem_ref = sii
                                };
                                i++;
                            }
                            peptide_id++;
                        }
                        protein_id++;
                    }
                    group_id++;
                }
            }
            TextWriter writer = new StreamWriter(Path.Combine(outputFolder, fileName + ".mzid"));
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
            SucessfullyFinishedWritingFile(Path.Combine(outputFolder, fileName + ".mzid"), nestedIds);
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected static List<Protein> LoadProteinDb(string fileName, bool generateDecoys, List<ModificationWithMass> localizeableModifications, bool isContaminant, out Dictionary<string, Modification> um)
        {
            if (Path.GetExtension(fileName).Equals(".fasta"))
            {
                um = null;
                return ProteinDbLoader.LoadProteinFasta(fileName, generateDecoys, isContaminant, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_gene_expression);
            }
            else
                return ProteinDbLoader.LoadProteinXML(fileName, generateDecoys, localizeableModifications, isContaminant, new List<string>(), out um);
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId);

        protected void WriteProteinGroupsToTsv(List<ProteinGroup> items, string outputFolder, string strippedFileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, strippedFileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(items.First().TabSeparatedHeader);
                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        protected void WritePeptideQuantificationResultsToTsv(List<FlashLFQ.FlashLFQSummedFeatureGroup> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    //output.WriteLine(items.First().TabSeparatedHeader);
                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        protected void WritePeakQuantificationResultsToTsv(List<FlashLFQ.FlashLFQFeature> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    //output.WriteLine(items.First().TabSeparatedHeader);
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

        private static List<Tuple<string, string>> GetModsFromString(string value)
        {
            return value.Split(new string[] { "\t\t" }, StringSplitOptions.None).Select(b => new Tuple<string, string>(b.Split('\t').First(), b.Split('\t').Last())).ToList();
        }

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