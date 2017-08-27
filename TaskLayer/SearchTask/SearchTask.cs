using Chemistry;
using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Serialization;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {
        #region Private Fields

        private const double binTolInDaltons = 0.003;

        private FlashLFQEngine FlashLfqEngine;

        #endregion Private Fields

        #region Public Constructors

        public SearchTask() : base(MyTask.Search)
        {
            CommonParameters = new CommonParameters();

            SearchParameters = new SearchParameters();

            FlashLfqEngine = new FlashLFQEngine();
        }

        #endregion Public Constructors

        #region Public Properties

        public SearchParameters SearchParameters { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static void WriteMzidentml(IEnumerable<Psm> items, List<ProteinGroup> groups, List<ModificationWithMass> variableMods, List<ModificationWithMass> fixedMods, List<Protease> proteases, double threshold, MassDiffAcceptor searchMode, Tolerance productTolerance, int missedCleavages, string outputPath)
        {
            List<PeptideWithSetModifications> peptides = items.SelectMany(i => i.MostProbableProteinInfo.PeptidesWithSetModifications).Distinct().ToList();
            List<Protein> proteins = peptides.Select(p => p.Protein).Distinct().ToList();
            List<string> filenames = items.Select(i => i.FullFilePath).Distinct().ToList();
            Dictionary<string, string> database_reference = new Dictionary<string, string>();
            List<string> databases = proteins.Select(p => p.DatabaseFilePath).Distinct().ToList();
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType()
            {
                version = "1.1.0",
                id = "",
            };

            //cvlist: URLs of controlled vocabularies used within the file.
            _mzid.cvList = new mzIdentML110.Generated.cvType[3] { new mzIdentML110.Generated.cvType()
            {
                id = "PSI-MS",
                fullName = "Proteomics Standards Initiative Mass Spectrometry Vocabularies",
                uri = "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",
                version= "4.0.9"
            },
             new mzIdentML110.Generated.cvType()
            {
                id = "PSI-MOD",
                fullName = "Proteomics Standards Initiative Modification Vocabularies",
                uri = "http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/mod/data/PSI-MOD.obo",
                version= "1.2"
            },
            new mzIdentML110.Generated.cvType()
            {
                id = "UO",
                fullName = "UNIT-ONTOLOGY",
                uri = "http://www.unimod.org/obo/unimod.obo"
            }
            };

            _mzid.AnalysisSoftwareList = new mzIdentML110.Generated.AnalysisSoftwareType[1] {  new mzIdentML110.Generated.AnalysisSoftwareType()
            {
                id = "AS_MetaMorpheus",
                name = "MetaMorpheus",
                version = GlobalEngineLevelSettings.MetaMorpheusVersion,
                uri = "https://github.com/smith-chem-wisc/MetaMorpheus",
                SoftwareName = new mzIdentML110.Generated.ParamType()
                {
                    Item = new mzIdentML110.Generated.CVParamType()
                    {
                        //using Morpheus's accession until we get MetaMorpheus entered
                        accession = "MS:1002661",
                        name = "Morpheus",
                        cvRef = "PSI-MS"
                    }
                }
            }};
            _mzid.DataCollection = new mzIdentML110.Generated.DataCollectionType()
            {
                AnalysisData = new mzIdentML110.Generated.AnalysisDataType()
                {
                    SpectrumIdentificationList = new mzIdentML110.Generated.SpectrumIdentificationListType[1]
        {
                        new mzIdentML110.Generated.SpectrumIdentificationListType()
                        {
                            id = "SIL",
                            SpectrumIdentificationResult = new mzIdentML110.Generated.SpectrumIdentificationResultType[items.Count()]
                        }
        }
                },
                Inputs = new mzIdentML110.Generated.InputsType()
                {
                    SearchDatabase = new mzIdentML110.Generated.SearchDatabaseType[databases.Count()],
                    SpectraData = new mzIdentML110.Generated.SpectraDataType[filenames.Count]
                }
            };

            _mzid.SequenceCollection = new mzIdentML110.Generated.SequenceCollectionType()
            {
                Peptide = new mzIdentML110.Generated.PeptideType[peptides.Count],
                DBSequence = new mzIdentML110.Generated.DBSequenceType[proteins.Count],
                PeptideEvidence = new mzIdentML110.Generated.PeptideEvidenceType[peptides.Count]
            };

            _mzid.AnalysisCollection = new mzIdentML110.Generated.AnalysisCollectionType()
            {
                SpectrumIdentification = new mzIdentML110.Generated.SpectrumIdentificationType[1]
                {
                    new mzIdentML110.Generated.SpectrumIdentificationType()
                    {
                        id = "SI",
                        spectrumIdentificationList_ref = "SIL",
                        spectrumIdentificationProtocol_ref = "SIP",
                        InputSpectra = new mzIdentML110.Generated.InputSpectraType[filenames.Count],
                        SearchDatabaseRef = new mzIdentML110.Generated.SearchDatabaseRefType[databases.Count]
                    }
                }
            };

            int database_index = 0;
            foreach (string database in databases)
            {
                _mzid.DataCollection.Inputs.SearchDatabase[database_index] = new mzIdentML110.Generated.SearchDatabaseType()
                {
                    id = "SDB_" + database_index,
                    location = database,
                    DatabaseName = new mzIdentML110.Generated.ParamType()
                    {
                        Item = new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1001073",
                            name = "database type amino acid",
                            cvRef = "PSI-MS"
                        }
                    }
                };
                database_reference.Add(database, "SDB_" + database_index);
                _mzid.AnalysisCollection.SpectrumIdentification[0].SearchDatabaseRef[database_index] = new mzIdentML110.Generated.SearchDatabaseRefType()
                {
                    searchDatabase_ref = "SDB_" + database_index
                };
                database_index++;
            }

            int protein_index = 0;
            foreach (Protein protein in proteins)
            {
                _mzid.SequenceCollection.DBSequence[protein_index] = new mzIdentML110.Generated.DBSequenceType()
                {
                    id = "DBS_" + protein.Accession,
                    lengthSpecified = true,
                    length = protein.Length,
                    searchDatabase_ref = database_reference[protein.DatabaseFilePath],
                    accession = protein.Accession,
                    Seq = protein.BaseSequence,
                    cvParam = new mzIdentML110.Generated.CVParamType[1]
                    {
                        new mzIdentML110.Generated.CVParamType()
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

            Dictionary<string, int> spectral_ids = new Dictionary<string, int>(); //key is datafile, value is datafile's id
            int spectra_data_id = 0;
            foreach (string data_filepath in filenames)
            {
                bool thermoRawFile = Path.GetExtension(data_filepath) == ".raw";
                string spectral_data_id = "SD_" + spectra_data_id;
                spectral_ids.Add(data_filepath, spectra_data_id);
                _mzid.AnalysisCollection.SpectrumIdentification[0].InputSpectra[spectra_data_id] = new mzIdentML110.Generated.InputSpectraType()
                {
                    spectraData_ref = spectral_data_id
                };
                _mzid.DataCollection.Inputs.SpectraData[spectra_data_id] = new mzIdentML110.Generated.SpectraDataType()
                {
                    id = spectral_data_id,
                    name = Path.GetFileNameWithoutExtension(data_filepath),
                    location = data_filepath,
                    FileFormat = new mzIdentML110.Generated.FileFormatType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
                        {
                            accession = thermoRawFile ? "MS:1000563" : "MS:1000584",
                            name = thermoRawFile ? "Thermo RAW format" : "mzML format",
                            cvRef = "PSI-MS"
                        }
                    },
                    SpectrumIDFormat = new mzIdentML110.Generated.SpectrumIDFormatType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
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
            foreach (Psm psm in items)
            {
                PeptideWithSetModifications peptide = psm.MostProbableProteinInfo.PeptidesWithSetModifications.First();
                //if first peptide on list hasn't been added, add peptide and peptide evidence
                if (!peptide_ids.TryGetValue(peptide, out Tuple<int, int, List<string>> peptide_id))
                {
                    peptide_id = new Tuple<int, int, List<string>>(p_index, 0, new List<string>());
                    p_index++;
                    _mzid.SequenceCollection.Peptide[peptide_id.Item1] = new mzIdentML110.Generated.PeptideType()
                    {
                        PeptideSequence = peptide.BaseSequence,
                        id = "P_" + peptide_id.Item1,
                        Modification = new mzIdentML110.Generated.ModificationType[peptide.NumMods]
                    };
                    int mod_id = 0;
                    foreach (KeyValuePair<int, ModificationWithMass> mod in peptide.allModsOneIsNterminus)
                    {
                        UsefulProteomicsDatabases.Generated.oboTerm psimod = null;
                        string name;
                        if (mod.Value.linksToOtherDbs.ContainsKey("PSI-MOD")) psimod = GlobalEngineLevelSettings.PsiModDeserialized.Items.OfType<UsefulProteomicsDatabases.Generated.oboTerm>().Where(m => m.id == mod.Value.linksToOtherDbs["PSI-MOD"].First()).FirstOrDefault();
                        name = psimod != null ? psimod.name : mod.Value.id;

                        _mzid.SequenceCollection.Peptide[peptide_id.Item1].Modification[mod_id] = new mzIdentML110.Generated.ModificationType()
                        {
                            location = mod.Key - 1,
                            locationSpecified = true,
                            monoisotopicMassDelta = mod.Value.monoisotopicMass,
                            residues = new string[1] { mod.Value.motif.ToString() },
                            monoisotopicMassDeltaSpecified = true,
                            cvParam = new mzIdentML110.Generated.CVParamType[1]
                            {
                            new mzIdentML110.Generated.CVParamType()
                            {
                                cvRef =  mod.Value.linksToOtherDbs.ContainsKey("PSI-MOD")? "PSI-MOD" : "PSI-MS",
                                name = mod.Value.linksToOtherDbs.ContainsKey("PSI-MOD")? name : "unknown modification",
                                accession =mod.Value.linksToOtherDbs.ContainsKey("PSI-MOD")? mod.Value.linksToOtherDbs["PSI-MOD"].First() : "MS:1001460",
                                value = mod.Value.linksToOtherDbs.ContainsKey("PSI-MOD")?  "" : mod.Value.id //give id of mod if unknown modification
                            }
                            }
                        };
                        mod_id++;
                    }

                    foreach (PeptideWithSetModifications peptide_evidence in psm.MostProbableProteinInfo.PeptidesWithSetModifications)
                    {
                        _mzid.SequenceCollection.PeptideEvidence[pe_index] = new mzIdentML110.Generated.PeptideEvidenceType()
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

                if (!psm_per_scan.TryGetValue(new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber), out Tuple<int, int> scan_result_scan_item)) //check to see if scan has already been added
                {
                    scan_result_scan_item = new Tuple<int, int>(sir_id, 0);
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1] = new mzIdentML110.Generated.SpectrumIdentificationResultType()
                    {
                        id = "SIR_" + scan_result_scan_item.Item1,
                        spectraData_ref = "SD_" + spectral_ids[psm.FullFilePath].ToString(),
                        spectrumID = psm.ScanNumber.ToString(),
                        SpectrumIdentificationItem = new mzIdentML110.Generated.SpectrumIdentificationItemType[500]
                    };
                    psm_per_scan.Add(new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber), scan_result_scan_item);
                    sir_id++;
                }
                else
                {
                    psm_per_scan[new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber)] = new Tuple<int, int>(scan_result_scan_item.Item1, scan_result_scan_item.Item2 + 1);
                    scan_result_scan_item = psm_per_scan[new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber)];
                }
                foreach (PeptideWithSetModifications p in psm.MostProbableProteinInfo.PeptidesWithSetModifications)
                {
                    peptide_ids[p].Item3.Add("SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2);
                }
                _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2] = new mzIdentML110.Generated.SpectrumIdentificationItemType()
                {
                    chargeState = psm.ScanPrecursorCharge,
                    id = "SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2,
                    experimentalMassToCharge = Math.Round(psm.ScanPrecursorMonoisotopicPeak.Mz, 5),
                    passThreshold = psm.FdrInfo.QValue <= threshold,
                    peptide_ref = "P_" + peptide_id.Item1,
                    PeptideEvidenceRef = new mzIdentML110.Generated.PeptideEvidenceRefType[psm.MostProbableProteinInfo.PeptidesWithSetModifications.Count],
                    cvParam = new mzIdentML110.Generated.CVParamType[2]
                    {
                        new mzIdentML110.Generated.CVParamType()
                        {
                            name = "Morpheus:Morpheus score",
                            cvRef = "PSI-MS",
                            accession = "MS:1002662",
                            value = psm.Score.ToString()
                        },
                        new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1002354",
                            name = "PSM-level q-value",
                            cvRef = "PSI-MS",
                            value = psm.FdrInfo.QValue.ToString()
                        }
                    }
                };
                if (psm.PeptideMonisotopicMass.HasValue)
                {
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2].calculatedMassToCharge = Math.Round(psm.PeptideMonisotopicMass.Value.ToMz(psm.ScanPrecursorCharge), 5);
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2].calculatedMassToChargeSpecified = true;
                }

                int pe = 0;
                foreach (PeptideWithSetModifications p in psm.MostProbableProteinInfo.PeptidesWithSetModifications)
                {
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2].PeptideEvidenceRef[pe]
                        = new mzIdentML110.Generated.PeptideEvidenceRefType()
                        {
                            peptideEvidence_ref = "PE_" + peptide_ids[p].Item2
                        };
                    pe++;
                }
            }

            _mzid.AnalysisProtocolCollection = new mzIdentML110.Generated.AnalysisProtocolCollectionType()
            {
                SpectrumIdentificationProtocol = new mzIdentML110.Generated.SpectrumIdentificationProtocolType[1]
                {
                    new mzIdentML110.Generated.SpectrumIdentificationProtocolType()
                    {
                        id = "SIP",
                        analysisSoftware_ref = "AS_MetaMorpheus",
                        SearchType = new mzIdentML110.Generated.ParamType()
                        {
                            Item = new mzIdentML110.Generated.CVParamType()
                            {
                                accession = "MS:1001083",
                                name = "ms-ms search",
                                cvRef = "PSI-MS"
                            }
                        },
                        AdditionalSearchParams = new mzIdentML110.Generated.ParamListType()
                        {
                            //TODO: ADD SEARCH PARAMS?
                            Items = new mzIdentML110.Generated.AbstractParamType[2]
                            {
                                new mzIdentML110.Generated.CVParamType()
                                {
                                    accession = "MS:1001211",
                                    cvRef = "PSI-MS",
                                    name = "parent mass type mono"
                                },
                                new mzIdentML110.Generated.CVParamType()
                                {
                                    accession = "MS:1001255",
                                    name = "fragment mass type mono",
                                    cvRef = "PSI-MS"
                                },
                            }
                        },
                        ModificationParams = new mzIdentML110.Generated.SearchModificationType[fixedMods.Count + variableMods.Count],
                        Enzymes = new mzIdentML110.Generated.EnzymesType()
                        {
                            Enzyme = new mzIdentML110.Generated.EnzymeType[proteases.Count]
                        },
                        FragmentTolerance = new mzIdentML110.Generated.CVParamType[2]
                        {
                            new mzIdentML110.Generated.CVParamType()
                            {
                                accession = "MS:1001412",
                                name = "search tolerance plus value",
                                value = productTolerance.Value.ToString(),
                                cvRef = "PSI-MS",

                                unitCvRef = "UO"
                            },
                            new mzIdentML110.Generated.CVParamType()
                            {
                                accession = "MS:1001413",
                                name = "search tolerance minus value",
                                value = productTolerance.Value.ToString(),
                                cvRef = "PSI-MS",
                                unitAccession = productTolerance is PpmTolerance? "UO:0000169": "UO:0000221",
                                unitName = productTolerance is PpmTolerance? "parts per million" : "dalton" ,
                                unitCvRef = "UO"
                            }
                        },
                        ParentTolerance = new mzIdentML110.Generated.CVParamType[1]
                        {
                            new mzIdentML110.Generated.CVParamType()
                            {
                                accession = "MS1001411",
                                name = "search tolerance specification",
                                cvRef = "PSI-MS",
                                value = searchMode.FileNameAddition
                            }
                        },
                        Threshold = new mzIdentML110.Generated.ParamListType()
                        {
                            Items = new mzIdentML110.Generated.CVParamType[1]
                            {
                                new mzIdentML110.Generated.CVParamType()
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
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].Enzymes.Enzyme[protease_index] = new mzIdentML110.Generated.EnzymeType()
                {
                    id = "E_" + protease_index,
                    name = protease.Name,
                    semiSpecific = protease.CleavageSpecificity == CleavageSpecificity.Semi,
                    missedCleavages = missedCleavages,
                    SiteRegexp = protease.SiteRegexp,
                    EnzymeName = new mzIdentML110.Generated.ParamListType()
                    {
                        Items = new mzIdentML110.Generated.AbstractParamType[1]
                        {
                            new mzIdentML110.Generated.CVParamType()
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
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = new mzIdentML110.Generated.SearchModificationType()
                {
                    fixedMod = true,
                    massDelta = (float)mod.monoisotopicMass,
                    residues = mod.motif.ToString(),
                };
                mod_index++;
            }
            foreach (ModificationWithMass mod in variableMods)
            {
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = new mzIdentML110.Generated.SearchModificationType()
                {
                    fixedMod = false,
                    massDelta = (float)mod.monoisotopicMass,
                    residues = mod.motif.ToString(),
                };
                mod_index++;
            }

            _mzid.AnalysisProtocolCollection.ProteinDetectionProtocol = new mzIdentML110.Generated.ProteinDetectionProtocolType()
            {
                id = "PDP",
                analysisSoftware_ref = "AS_MetaMorpheus",
                Threshold = new mzIdentML110.Generated.ParamListType()
                {
                    Items = new mzIdentML110.Generated.CVParamType[1]
                    {
                        new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1001448",
                            name = "pep:FDR threshold",
                            cvRef = "PSI-MS",
                            value = threshold.ToString()
                        }
                    }
                }
            };

            if (groups != null)
            {
                _mzid.DataCollection.AnalysisData.ProteinDetectionList = new mzIdentML110.Generated.ProteinDetectionListType()
                {
                    id = "PDL",
                    ProteinAmbiguityGroup = new mzIdentML110.Generated.ProteinAmbiguityGroupType[groups.Count]
                };

                int group_id = 0;
                int protein_id = 0;
                foreach (ProteinGroup proteinGroup in groups)
                {
                    _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id] = new mzIdentML110.Generated.ProteinAmbiguityGroupType()
                    {
                        id = "PAG_" + group_id,
                        ProteinDetectionHypothesis = new mzIdentML110.Generated.ProteinDetectionHypothesisType[proteinGroup.Proteins.Count]
                    };
                    int pag_protein_index = 0;
                    foreach (Protein protein in proteinGroup.Proteins)
                    {
                        _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[pag_protein_index] = new mzIdentML110.Generated.ProteinDetectionHypothesisType()
                        {
                            id = "PDH_" + protein_id,
                            dBSequence_ref = "DBS_" + protein.Accession,
                            passThreshold = proteinGroup.QValue <= threshold,
                            PeptideHypothesis = new mzIdentML110.Generated.PeptideHypothesisType[proteinGroup.AllPeptides.Count],
                            cvParam = new mzIdentML110.Generated.CVParamType[4]
                            {
                            new mzIdentML110.Generated.CVParamType()
                            {
                                accession = "MS:1002663",
                                name = "Morpheus:summed Morpheus score",
                                cvRef = "PSI-MS",
                                value = proteinGroup.ProteinGroupScore.ToString()
                            },
                            new mzIdentML110.Generated.CVParamType()
                            {
                                accession = "MS1002373",
                                name = "protein group-level q-value",
                                value = proteinGroup.QValue.ToString()
                            },
                            new mzIdentML110.Generated.CVParamType()
                            {
                                accession =  "MS:1001093",
                                name = "sequence coverage",
                                cvRef = "PSI-MS",
                                value = proteinGroup.SequenceCoveragePercent.First().ToString()
                            },
                            new mzIdentML110.Generated.CVParamType()
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
                            if (peptide.Protein == protein)
                            {
                                _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[pag_protein_index].PeptideHypothesis[peptide_id] = new mzIdentML110.Generated.PeptideHypothesisType()
                                {
                                    peptideEvidence_ref = "PE_" + peptide_ids[peptide].Item2,
                                    SpectrumIdentificationItemRef = new mzIdentML110.Generated.SpectrumIdentificationItemRefType[peptide_ids[peptide].Item3.Count],
                                };
                                int i = 0;
                                foreach (string sii in peptide_ids[peptide].Item3)
                                {
                                    _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[pag_protein_index].PeptideHypothesis[peptide_id].SpectrumIdentificationItemRef[i] = new mzIdentML110.Generated.SpectrumIdentificationItemRefType()
                                    {
                                        spectrumIdentificationItem_ref = sii
                                    };
                                    i++;
                                }
                                peptide_id++;
                            }
                        }
                        pag_protein_index++;
                        protein_id++;
                    }
                    group_id++;
                }
            }
            XmlWriter writer = XmlWriter.Create(outputPath, settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        public static CommonParameters SetAllFileSpecificCommonParams(CommonParameters commonParams, FileSpecificSettings currentFileSpecificSettings)
        {
            if (currentFileSpecificSettings == null)
                return commonParams;

            CommonParameters returnParams = new CommonParameters
            {
                DoPrecursorDeconvolution = currentFileSpecificSettings.DoPrecursorDeconvolution ?? commonParams.DoPrecursorDeconvolution,

                UseProvidedPrecursorInfo = currentFileSpecificSettings.UseProvidedPrecursorInfo ?? commonParams.UseProvidedPrecursorInfo,

                DeconvolutionIntensityRatio = currentFileSpecificSettings.DeconvolutionIntensityRatio ?? commonParams.DeconvolutionIntensityRatio,

                DeconvolutionMaxAssumedChargeState = currentFileSpecificSettings.DeconvolutionMaxAssumedChargeState ?? commonParams.DeconvolutionMaxAssumedChargeState,

                DeconvolutionMassTolerance = currentFileSpecificSettings.DeconvolutionMassTolerance ?? commonParams.DeconvolutionMassTolerance,

                TotalPartitions = currentFileSpecificSettings.TotalPartitions ?? commonParams.TotalPartitions,

                ProductMassTolerance = currentFileSpecificSettings.ProductMassTolerance ?? commonParams.ProductMassTolerance,

                ConserveMemory = currentFileSpecificSettings.ConserveMemory ?? commonParams.ConserveMemory,

                ScoreCutoff = currentFileSpecificSettings.ScoreCutoff ?? commonParams.ScoreCutoff,

                ReportAllAmbiguity = currentFileSpecificSettings.ReportAllAmbiguity ?? commonParams.ReportAllAmbiguity,

                ExcelCompatible = currentFileSpecificSettings.ExcelCompatible ?? commonParams.ExcelCompatible 
            };

            returnParams.DigestionParams.InitiatorMethionineBehavior = currentFileSpecificSettings.InitiatorMethionineBehavior.Equals(InitiatorMethionineBehavior.Undefined) ? commonParams.DigestionParams.InitiatorMethionineBehavior : currentFileSpecificSettings.InitiatorMethionineBehavior;
            returnParams.DigestionParams.MaxMissedCleavages = currentFileSpecificSettings.MaxMissedCleavages ?? commonParams.DigestionParams.MaxMissedCleavages;
            returnParams.DigestionParams.MinPeptideLength = currentFileSpecificSettings.MinPeptideLength ?? commonParams.DigestionParams.MinPeptideLength;
            returnParams.DigestionParams.MaxPeptideLength = currentFileSpecificSettings.MaxPeptideLength ?? commonParams.DigestionParams.MaxPeptideLength;
            returnParams.DigestionParams.MaxModificationIsoforms = currentFileSpecificSettings.MaxModificationIsoforms ?? commonParams.DigestionParams.MaxModificationIsoforms;
            returnParams.DigestionParams.Protease = currentFileSpecificSettings.Protease ?? commonParams.DigestionParams.Protease;
            returnParams.DigestionParams.MaxModsForPeptide = currentFileSpecificSettings.Max_mods_for_peptide ?? commonParams.DigestionParams.MaxModsForPeptide;

            return returnParams;
        }

        #endregion Public Methods

        #region Protected Methods

        protected static void WriteTree(BinTreeStructure myTreeStructure, string writtenFile)
        {
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
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this);
            List<Psm>[] allPsms = new List<Psm>[SearchParameters.MassDiffAcceptors.Count];
            for (int searchModeIndex = 0; searchModeIndex < SearchParameters.MassDiffAcceptors.Count; searchModeIndex++)
                allPsms[searchModeIndex] = new List<Psm>();

            Status("Loading modifications...", taskId);

            #region Load modifications

            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (CommonParameters.LocalizeAll)
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            #endregion Load modifications

            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons && SearchParameters.AddCompIons)
                ionTypes.Add(ProductType.B);
            else if (CommonParameters.BIons)
                ionTypes.Add(ProductType.BnoB1ions);
            if (CommonParameters.YIons)
                ionTypes.Add(ProductType.Y);
            if (CommonParameters.ZdotIons)
                ionTypes.Add(ProductType.Zdot);
            if (CommonParameters.CIons)
                ionTypes.Add(ProductType.C);
            TerminusType terminusType = ProductTypeToTerminusType.IdentifyTerminusType(ionTypes);
            Status("Loading proteins...", new List<string> { taskId });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, SearchParameters.SearchDecoy, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();

            proseCreatedWhileRunning.Append("The following search settings were used: ");
            proseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            if (CommonParameters.DigestionParams.MaxPeptideLength == null)
            {
                proseCreatedWhileRunning.Append("maximum peptide length = unspecified; ");
            }
            else
            {
                proseCreatedWhileRunning.Append("maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            }
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            proseCreatedWhileRunning.Append("parent mass tolerance(s) = {" + String.Join("; ", SearchParameters.MassDiffAcceptors.Select(m => m.ToProseString())) + "}; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");
            proseCreatedWhileRunning.Append("report all ambiguity = " + CommonParameters.ReportAllAmbiguity + "; ");


            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxDegreeOfParallelism.Value;
            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);

            HashSet<DigestionParams> ListOfDigestionParams = GetListOfDistinctDigestionParams(CommonParameters, fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b)));
            ListOfDigestionParams.ToList().ForEach(d => d.Protease.IsSemiSpecific(SearchParameters.SearchType, terminusType));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);
            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                Psm[][] fileSpecificPsms = new Psm[SearchParameters.MassDiffAcceptors.Count()][];

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                for (int aede = 0; aede < SearchParameters.MassDiffAcceptors.Count; aede++)
                    fileSpecificPsms[aede] = new Psm[arrayOfMs2ScansSortedByMass.Length];

                if (SearchParameters.SearchType == SearchType.Modern || SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                    {
                        List<CompactPeptide> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                        float[] keys = null;
                        List<int>[] fragmentIndex = null;

                        #region Generate indices for modern search

                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, ionTypes, currentPartition, SearchParameters.SearchDecoy, ListOfDigestionParams, CommonParameters.TotalPartitions, new List<string> { taskId });
                        Dictionary<float, List<int>> fragmentIndexDict;
                        lock (indexLock)
                        {
                            string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);

                            if (pathToFolderWithIndices == null)
                            {
                                var output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                                Status("Writing params...", new List<string> { taskId });
                                var paramsFile = Path.Combine(output_folderForIndices, "indexEngine.params");
                                WriteIndexEngineParams(indexEngine, paramsFile);
                                SucessfullyFinishedWritingFile(paramsFile, new List<string> { taskId });

                                Status("Running Index Engine...", new List<string> { taskId });
                                var indexResults = (IndexingResults)indexEngine.Run();
                                peptideIndex = indexResults.PeptideIndex;
                                fragmentIndexDict = indexResults.FragmentIndexDict;

                                Status("Writing peptide index...", new List<string> { taskId });
                                var peptideIndexFile = Path.Combine(output_folderForIndices, "peptideIndex.ind");
                                WritePeptideIndex(peptideIndex, peptideIndexFile);
                                SucessfullyFinishedWritingFile(peptideIndexFile, new List<string> { taskId });

                                Status("Writing fragment index...", new List<string> { taskId });
                                var fragmentIndexFile = Path.Combine(output_folderForIndices, "fragmentIndex.ind");
                                WriteFragmentIndexNetSerializer(fragmentIndexDict, fragmentIndexFile);
                                SucessfullyFinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
                            }
                            else
                            {
                                Status("Reading peptide index...", new List<string> { taskId });
                                var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
                                var ser = new NetSerializer.Serializer(messageTypes);
                                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
                                    peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);

                                Status("Reading fragment index...", new List<string> { taskId });
                                messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
                                ser = new NetSerializer.Serializer(messageTypes);
                                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
                                    fragmentIndexDict = (Dictionary<float, List<int>>)ser.Deserialize(file);
                            }
                        }
                        keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
                        fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

                        #endregion Generate indices for modern search

                        Status("Searching files...", taskId);
                        if (SearchParameters.SearchType == SearchType.NonSpecific)
                            new NonSpecificEnzymeEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ionTypes, currentPartition, combinedParams, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, thisId).Run();
                        else//if(SearchType==SearchType.Modern)
                            new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ionTypes, currentPartition, combinedParams, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, thisId).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                    }
                }
                else //If classic search
                {
                    for (int aede = 0; aede < SearchParameters.MassDiffAcceptors.Count; aede++)
                        fileSpecificPsms[aede] = new Psm[arrayOfMs2ScansSortedByMass.Length];

                    Status("Starting search...", thisId);
                    new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ionTypes, SearchParameters.MassDiffAcceptors, SearchParameters.AddCompIons, combinedParams, thisId).Run();

                    myFileManager.DoneWithFile(origDataFile);

                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }

                lock (psmLock)
                {
                    for (int searchModeIndex = 0; searchModeIndex < SearchParameters.MassDiffAcceptors.Count; searchModeIndex++)
                        allPsms[searchModeIndex].AddRange(fileSpecificPsms[searchModeIndex]);
                }

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            });
            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            // Group and order psms

            Status("Matching peptides to proteins...", taskId);
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            if (SearchParameters.SearchType == SearchType.NonSpecific)
            {
                NonSpecificEnzymeSequencesToActualPeptides sequencesToActualProteinPeptidesEngine = new NonSpecificEnzymeSequencesToActualPeptides(allPsms, proteinList, fixedModifications, variableModifications, terminusType, ListOfDigestionParams, SearchParameters.MassDiffAcceptors, new List<string> { taskId });
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;
            }
            else
            {
                SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, terminusType, ListOfDigestionParams, new List<string> { taskId });
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;
            }

            ProteinParsimonyResults proteinAnalysisResults = null;
            if (SearchParameters.DoParsimony)
                proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, SearchParameters.ModPeptidesAreUnique, new List<string> { taskId }).Run());

            Status("Resolving most probable peptide...", new List<string> { taskId });
            for (int j = 0; j < SearchParameters.MassDiffAcceptors.Count; j++)
            {
                if (allPsms[j] != null)
                {
                    foreach (var huh in allPsms[j])
                    {
                        if (huh != null && huh.MostProbableProteinInfo == null)
                            huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
                    }
                }
            }

            Status("Ordering and grouping psms...", taskId);
            for (int j = 0; j < allPsms.Length; j++)
                allPsms[j] = allPsms[j].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            Status("Running FDR analysis...", taskId);
            var fdrAnalysisResults = new FdrAnalysisEngine(allPsms, SearchParameters.MassDiffAcceptors, new List<string> { taskId }).Run();

            List<ProteinGroup>[] proteinGroupsHere = new List<ProteinGroup>[SearchParameters.MassDiffAcceptors.Count];

            if (SearchParameters.DoParsimony)
            {
                for (int j = 0; j < SearchParameters.MassDiffAcceptors.Count; j++)
                {
                    var ressdf = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, allPsms[j], SearchParameters.MassDiffAcceptors, SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreUnique, new List<string> { taskId }).Run();
                    proteinGroupsHere[j] = ressdf.sortedAndScoredProteinGroups;
                }
            }

            if (SearchParameters.DoLocalizationAnalysis)
            {
                Status("Running localization analysis...", taskId);
                Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    Status("Running localization analysis...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile);
                    var localizationEngine = new LocalizationEngine(allPsms.SelectMany(b => b).Where(b => b != null && b.FullFilePath.Equals(origDataFile)).ToList(), ionTypes, myMsDataFile, combinedParams.ProductMassTolerance, new List<string> { taskId, "Individual Spectra Files", origDataFile }, this.SearchParameters.AddCompIons);
                    localizationEngine.Run();
                    myFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                });
                Status("Grouping by matched ions...", taskId);
                for (int j = 0; j < allPsms.Length; j++)
                    allPsms[j] = allPsms[j].GroupBy(b => new Tuple<string, int, MatchedIonMassesListPositiveIsMatch>(b.FullFilePath, b.ScanNumber, b.MatchedIonDictPositiveIsMatch)).Select(b => b.First()).ToList();
            }

            new ModificationAnalysisEngine(allPsms, SearchParameters.MassDiffAcceptors.Count, new List<string> { taskId }).Run();

            if (SearchParameters.DoQuantification)
            {
                // use FlashLFQ to quantify peaks across all files
                Status("Quantifying...", taskId);
                FlashLfqEngine.PassFilePaths(currentRawFileList.ToArray());

                if (!FlashLfqEngine.ReadPeriodicTable(GlobalEngineLevelSettings.elementsLocation))
                    throw new MetaMorpheusException("Quantification error - could not find periodic table file");

                if (!FlashLfqEngine.ParseArgs(new string[] {
                        "--ppm " + SearchParameters.QuantifyPpmTol,
                        "--sil true",
                        "--pau false",
                        "--mbr " + SearchParameters.MatchBetweenRuns }
                    ))
                    throw new MetaMorpheusException("Quantification error - Could not pass parameters to quantification engine");

                var psmsBelowOnePercentFdr = allPsms.SelectMany(v => v).Where(p => p.FdrInfo.QValue < 0.01 && !p.IsDecoy);
                foreach (var psm in psmsBelowOnePercentFdr.Where(b => b.FullSequence != null && b.PeptideMonisotopicMass != null))
                    FlashLfqEngine.AddIdentification(Path.GetFileNameWithoutExtension(psm.FullFilePath), psm.BaseSequence, psm.FullSequence, psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, string.Join("|", psm.MostProbableProteinInfo.PeptidesWithSetModifications.Select(v => v.Protein.Accession).Distinct().OrderBy(v => v)));

                FlashLfqEngine.ConstructBinsFromIdentifications();

                Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    Status("Quantifying...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile);
                    FlashLfqEngine.Quantify(myMsDataFile, origDataFile);
                    myFileManager.DoneWithFile(origDataFile);
                    GC.Collect();
                    ReportProgress(new ProgressEventArgs(100, "Done quantifying!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                });

                if (FlashLfqEngine.mbr)
                    FlashLfqEngine.RetentionTimeCalibrationAndErrorCheckMatchedFeatures();

                // assign quantities to PSMs
                Dictionary<string, List<Psm>> baseseqToPsm = new Dictionary<string, List<Psm>>();
                List<Psm> list;
                foreach (var psm in psmsBelowOnePercentFdr.Where(b => b.BaseSequence != null))
                {
                    if (baseseqToPsm.TryGetValue(psm.BaseSequence, out list))
                        list.Add(psm);
                    else
                        baseseqToPsm.Add(psm.BaseSequence, new List<Psm>() { psm });
                }

                var summedPeaks = FlashLfqEngine.SumFeatures(FlashLfqEngine.allFeaturesByFile.SelectMany(p => p).ToList(), "BaseSequence");
                foreach (var summedPeak in summedPeaks)
                {
                    if (baseseqToPsm.TryGetValue(summedPeak.BaseSequence, out list))
                    {
                        var psmsForThisBaseSeqAndFile = list.GroupBy(p => p.FullFilePath);
                        foreach (var file in psmsForThisBaseSeqAndFile)
                        {
                            int j = Array.IndexOf(FlashLfqEngine.filePaths, file.Key);

                            foreach (var psm in file)
                                psm.QuantIntensity = summedPeak.intensitiesByFile[j];
                        }
                    }
                }
            }

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            if (SearchParameters.DoHistogramAnalysis)
            {
                for (int j = 0; j < SearchParameters.MassDiffAcceptors.Count; j++)
                {
                    var limitedpsms_with_fdr = allPsms[j].Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                    if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                    {
                        Status("Running histogram analysis...", new List<string> { taskId });
                        var myTreeStructure = new BinTreeStructure();
                        myTreeStructure.GenerateBins(limitedpsms_with_fdr, binTolInDaltons);
                        var writtenFile = Path.Combine(OutputFolder, "aggregate_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ".mytsv");
                        WriteTree(myTreeStructure, writtenFile);
                        SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });
                    }
                }
            }

            // Now that we are done with fdr analysis and localization analysis, can write the results!
            Status("Writing results...", taskId);
            for (int j = 0; j < SearchParameters.MassDiffAcceptors.Count; j++)
            {
                {
                    var writtenFile = Path.Combine(OutputFolder, "aggregatePSMs_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ".psmtsv");
                    WritePsmsToTsv(allPsms[j], writtenFile);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });
                    myTaskResults.AddNiceText("All target PSMS within 1% FDR " + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ": " + allPsms[j].Count(a => a.FdrInfo.QValue <= .01 && a.IsDecoy == false));
                }

                var uniquePeptides = allPsms[j].GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                {
                    var writtenFile = Path.Combine(OutputFolder, "aggregateUniquePeptides_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ".psmtsv");
                    WritePsmsToTsv(uniquePeptides, writtenFile);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });
                    myTaskResults.AddNiceText("Unique peptides within 1% FDR " + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ": " + uniquePeptides.Count(a => a.FdrInfo.QValue <= .01 && a.IsDecoy == false));
                }

                var psmsGroupedByFile = allPsms[j].GroupBy(p => p.FullFilePath);

                // individual psm files (with global psm fdr, global parsimony)
                foreach (var group in psmsGroupedByFile)
                {
                    var psmsForThisFile = group.ToList();

                    var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                    {
                        var writtenFile = Path.Combine(OutputFolder, strippedFileName + "_PSMs_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ".psmtsv");
                        WritePsmsToTsv(psmsForThisFile, writtenFile);
                        SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                        myTaskResults.AddNiceText("PSMs within 1% FDR " + SearchParameters.MassDiffAcceptors[j].FileNameAddition + " " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue <= .01 && a.IsDecoy == false));
                    }

                    {
                        var uniquePeptidesForFile = psmsForThisFile.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                        var writtenFile = Path.Combine(OutputFolder, strippedFileName + "_UniquePeptides_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ".psmtsv");
                        WritePsmsToTsv(uniquePeptidesForFile, writtenFile);
                        SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                        myTaskResults.AddNiceText("Unique peptides within 1% FDR " + SearchParameters.MassDiffAcceptors[j].FileNameAddition + " " + strippedFileName + ": " + uniquePeptidesForFile.Count(a => a.FdrInfo.QValue <= .01 && a.IsDecoy == false));
                    }
                }

                if (SearchParameters.DoParsimony)
                {
                    // aggregate protein group file
                    foreach (var pg in proteinGroupsHere[j])
                    {
                        if (pg.ProteinGroupScore != 0)
                            pg.AggregateQuantifyHelper(currentRawFileList);
                    }

                    WriteProteinGroupsToTsv(proteinGroupsHere[j], OutputFolder, "aggregateProteinGroups_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId }, psmsGroupedByFile.Select(b => b.Key).ToList());

                    // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                    foreach (var fullFilePath in currentRawFileList)
                    {
                        List<Psm> psmsForThisFile = psmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();

                        var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                        var subsetProteinGroupsForThisFile = new List<ProteinGroup>();
                        foreach (var pg in proteinGroupsHere[j])
                        {
                            var subsetPg = pg.ConstructSubsetProteinGroup(fullFilePath);
                            subsetPg.Score();

                            if (subsetPg.ProteinGroupScore != 0)
                            {
                                subsetPg.CalculateSequenceCoverage();
                                subsetPg.Quantify();
                                subsetProteinGroupsForThisFile.Add(subsetPg);
                            }
                        }
                        new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile, SearchParameters.MassDiffAcceptors, SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreUnique, new List<string> { taskId, "Individual Spectra Files", fullFilePath }).Run();
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, OutputFolder, strippedFileName + "_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + "_ProteinGroups", new List<string> { taskId, "Individual Spectra Files", fullFilePath }, null);

                        //Status("Writing mzid...", new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                        //var mzidFilePath = Path.Combine(OutputFolder, strippedFileName + "_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + ".mzid");
                        //WriteMzidentml(psmsForThisFile, subsetProteinGroupsForThisFile, variableModifications, fixedModifications, new List<Protease> { CommonParameters.DigestionParams.Protease }, 0.01, SearchParameters.MassDiffAcceptors[j], CommonParameters.ProductMassTolerance, CommonParameters.DigestionParams.MaxMissedCleavages, mzidFilePath);
                        //SucessfullyFinishedWritingFile(mzidFilePath, new List<string> { taskId, "Individual Spectra Files", fullFilePath });

                        ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", fullFilePath }));
                    }
                }

                if (SearchParameters.DoQuantification)
                {
                    foreach (var fullFilePath in currentRawFileList)
                    {
                        var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);
                        var peaksForThisFile = FlashLfqEngine.allFeaturesByFile[Array.IndexOf(FlashLfqEngine.filePaths, fullFilePath)];

                        WritePeakQuantificationResultsToTsv(peaksForThisFile, OutputFolder, strippedFileName + "_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition + "_QuantifiedPeaks", new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                    }

                    var summedPeaksByPeptide = FlashLfqEngine.SumFeatures(FlashLfqEngine.allFeaturesByFile.SelectMany(p => p).ToList(), "BaseSequence");
                    WritePeptideQuantificationResultsToTsv(summedPeaksByPeptide.ToList(), OutputFolder, "aggregateQuantifiedPeptides_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });

                    summedPeaksByPeptide = FlashLfqEngine.SumFeatures(FlashLfqEngine.allFeaturesByFile.SelectMany(p => p).ToList(), "FullSequence");
                    WritePeptideQuantificationResultsToTsv(summedPeaksByPeptide.ToList(), OutputFolder, "aggregateQuantifiedPeptidesByFullSeq_" + SearchParameters.MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });
                }
            }

            if (SearchParameters.WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { taskId });

                List<Modification> modificationsToAlwaysKeep = new List<Modification>();
                if (SearchParameters.KeepAllUniprotMods)
                    modificationsToAlwaysKeep.AddRange(GlobalEngineLevelSettings.AllModsKnown.Where(b => b.modificationType.Equals("Uniprot")));

                var goodPsmsForEachProtein = allPsms.SelectMany(b => b).Where(b => b.FdrInfo.QValueNotch < 0.01 && !b.IsDecoy && b.FullSequence != null && b.ProteinAccesion != null).GroupBy(b => b.CompactPeptides.First().Value.Item2.First().Protein).ToDictionary(b => b.Key);

                foreach (var protein in proteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        HashSet<Tuple<int, ModificationWithMass>> modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>();
                        if (goodPsmsForEachProtein.ContainsKey(protein))
                            modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>(goodPsmsForEachProtein[protein].SelectMany(b => b.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Select(c => new Tuple<int, ModificationWithMass>(GetOneBasedIndexInProtein(c.Key, b.MostProbableProteinInfo.PeptidesWithSetModifications.First()), c.Value))));

                        IDictionary<int, List<Modification>> modsToWrite = new Dictionary<int, List<Modification>>();
                        foreach (var modd in protein.OneBasedPossibleLocalizedModifications)
                            foreach (var mod in modd.Value)
                            {
                                if (modificationsToAlwaysKeep.Contains(mod as Modification)
                                    || modsObservedOnThisProtein.Contains(new Tuple<int, ModificationWithMass>(modd.Key, mod as ModificationWithMass)))
                                {
                                    if (!modsToWrite.ContainsKey(modd.Key))
                                        modsToWrite.Add(modd.Key, new List<Modification> { mod });
                                    else
                                        modsToWrite[modd.Key].Add(mod);
                                }
                            }
                        protein.OneBasedPossibleLocalizedModifications.Clear();
                        foreach (var kvp in modsToWrite)
                            protein.OneBasedPossibleLocalizedModifications.Add(kvp);
                    }
                }

                //writes all proteins
                if (dbFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
                }
                if (dbFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
                }

                //writes only detected proteins
                if (dbFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), goodPsmsForEachProtein.Keys.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
                }
                if (dbFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), goodPsmsForEachProtein.Keys.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
                }
            }

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            yield return type;
        }

        private static bool SameSettings(string pathToOldParamsFile, IndexingEngine indexEngine)
        {
            using (StreamReader reader = new StreamReader(pathToOldParamsFile))
                if (reader.ReadToEnd().Equals(indexEngine.ToString()))
                    return true;
            return false;
        }

        private static void WritePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFile))
            {
                ser.Serialize(file, peptideIndex);
            }
        }

        private static void WriteFragmentIndexNetSerializer(Dictionary<float, List<int>> fragmentIndex, string fragmentIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);
        }

        private static string GetExistingFolderWithIndices(IndexingEngine indexEngine, List<DbForTask> dbFilenameList)
        {
            // In every database location...
            foreach (var ok in dbFilenameList)
            {
                var baseDir = Path.GetDirectoryName(ok.FilePath);
                var directory = new DirectoryInfo(baseDir);
                DirectoryInfo[] directories = directory.GetDirectories();

                // Look at every subdirectory...
                foreach (DirectoryInfo possibleFolder in directories)
                {
                    if (File.Exists(Path.Combine(possibleFolder.FullName, "indexEngine.params")) &&
                        File.Exists(Path.Combine(possibleFolder.FullName, "peptideIndex.ind")) &&
                        File.Exists(Path.Combine(possibleFolder.FullName, "fragmentIndex.ind")) &&
                        SameSettings(Path.Combine(possibleFolder.FullName, "indexEngine.params"), indexEngine))
                        return possibleFolder.FullName;
                }
            }
            return null;
        }

        private static void WriteIndexEngineParams(IndexingEngine indexEngine, string fileName)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
        }

        private static string GenerateOutputFolderForIndices(List<DbForTask> dbFilenameList)
        {
            var folder = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FilePath), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            if (!Directory.Exists(folder))
                Directory.CreateDirectory(folder);
            return folder;
        }

        private static int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
                return peptideWithSetModifications.OneBasedStartResidueInProtein;
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
                return peptideWithSetModifications.OneBasedEndResidueInProtein;
            return peptideWithSetModifications.OneBasedStartResidueInProtein + oneIsNterminus - 2;
        }

        private HashSet<DigestionParams> GetListOfDistinctDigestionParams(CommonParameters commonParameters, IEnumerable<CommonParameters> enumerable)
        {
            HashSet<DigestionParams> okay = new HashSet<DigestionParams>
            {
                commonParameters.DigestionParams
            };

            foreach (var hah in enumerable)
                okay.Add(hah.DigestionParams);

            return okay;
        }

        private void WriteProteinGroupsToTsv(List<ProteinGroup> items, string outputFolder, string strippedFileName, List<string> nestedIds, List<string> FileNames)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, strippedFileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(ProteinGroup.GetTabSeparatedHeader(FileNames));
                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        private void WritePeptideQuantificationResultsToTsv(List<FlashLFQ.FlashLFQSummedFeatureGroup> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(FlashLFQ.FlashLFQSummedFeatureGroup.TabSeparatedHeader);

                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        private void WritePeakQuantificationResultsToTsv(List<FlashLFQ.FlashLFQFeature> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(FlashLFQ.FlashLFQFeature.TabSeparatedHeader);

                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }
                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        #endregion Private Methods
    }
}