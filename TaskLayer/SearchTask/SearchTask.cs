using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using EngineLayer.Indexing;
using EngineLayer.Localization;
using EngineLayer.ModernSearch;
using EngineLayer.ModificationAnalysis;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Serialization;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {
        #region Public Constructors

        public SearchTask() : base(MyTask.Search)
        {
            CommonParameters = new CommonParameters();

            SearchParameters = new SearchParameters();
        }

        #endregion Public Constructors

        #region Public Properties

        public SearchParameters SearchParameters { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static void WriteMzidentml(IEnumerable<Psm> items, List<EngineLayer.ProteinGroup> groups, List<ModificationWithMass> variableMods, List<ModificationWithMass> fixedMods, List<Protease> proteases, double threshold, Tolerance productTolerance, Tolerance parentTolerance, int missedCleavages, string outputPath)
        {
            List<PeptideWithSetModifications> peptides = items.SelectMany(i => i.CompactPeptides.SelectMany(c => c.Value.Item2)).Distinct().ToList();
            List<Protein> proteins = peptides.Select(p => p.Protein).Distinct().ToList();
            List<string> filenames = items.Select(i => i.FullFilePath).Distinct().ToList();
            Dictionary<string, string> database_reference = new Dictionary<string, string>();
            List<string> databases = proteins.Select(p => p.DatabaseFilePath).Distinct().ToList();

            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType110()
            {
                version = "1.1.0",
                id = "",
            };

            _mzid.Provider = new mzIdentML110.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML110.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML110.Generated.RoleType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            _mzid.AuditCollection = new mzIdentML110.Generated.AbstractContactType[2];

            _mzid.AuditCollection[0] = new mzIdentML110.Generated.PersonType()
            {
                id = "UWMadisonSmithGroupPerson",
                cvParam = new mzIdentML110.Generated.CVParamType[2]
                {
                    new mzIdentML110.Generated.CVParamType()
                    {
                        accession="MS:1000589",
                        name ="contact email",
                        cvRef ="PSI-MS",
                        value ="mm_support@chem.wisc.edu"
                    },

                       new mzIdentML110.Generated.CVParamType()
                    {
                        accession="MS:1000590",
                        name ="affiliation name",
                        cvRef ="PSI-MS",
                        value ="UWMadisonSmithGroup"
                    }
                }
            };

            _mzid.AuditCollection[1] = new mzIdentML110.Generated.OrganizationType()
            {
                id = "UWMadisonSmithGroup",

                cvParam = new mzIdentML110.Generated.CVParamType[2]
                {
                    new mzIdentML110.Generated.CVParamType()
                    {
                        accession="MS:1000589",
                        name ="contact email",
                        cvRef ="PSI-MS",
                        value ="mm_support@chem.wisc.edu"
                    },

                     new mzIdentML110.Generated.CVParamType()
                    {
                        accession="MS:1000590",
                        name ="affiliation name",
                        cvRef ="PSI-MS",
                        value ="UWMadisonSmithGroup"
                    }
                }
            };

            //cvlist: URLs of controlled vocabularies used within the file.
            _mzid.cvList = new mzIdentML110.Generated.cvType[4] { new mzIdentML110.Generated.cvType()
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
                id = "UNIMOD",
                fullName = "UNIT-ONTOLOGY",
                uri = "http://www.unimod.org/obo/unimod.obo"
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
                version = GlobalVariables.MetaMorpheusVersion,
                uri = "https://github.com/smith-chem-wisc/MetaMorpheus",
                SoftwareName = new mzIdentML110.Generated.ParamType()
                {
                    Item = new mzIdentML110.Generated.CVParamType
                    {
                        accession = "MS:1002826",
                        name = "MetaMorpheus",
                        cvRef = "PSI-MS"
                    }
                },

                ContactRole = new mzIdentML110.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML110.Generated.RoleType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
                        {
                             accession = "MS:1001267",
                             name="software vendor",
                             cvRef="PSI-MS"
                        }
                    }
                }
            }};
            _mzid.DataCollection = new mzIdentML110.Generated.DataCollectionType
            {
                AnalysisData = new mzIdentML110.Generated.AnalysisDataType()
                {
                    SpectrumIdentificationList = new mzIdentML110.Generated.SpectrumIdentificationListType[1]
        {
                        new mzIdentML110.Generated.SpectrumIdentificationListType
                        {
                            id = "SIL",
                            SpectrumIdentificationResult = new mzIdentML110.Generated.SpectrumIdentificationResultType[items.Count()]
                        }
        }
                },
                Inputs = new mzIdentML110.Generated.InputsType
                {
                    SearchDatabase = new mzIdentML110.Generated.SearchDatabaseType[databases.Count()],
                    SpectraData = new mzIdentML110.Generated.SpectraDataType[filenames.Count]
                }
            };

            _mzid.SequenceCollection = new mzIdentML110.Generated.SequenceCollectionType
            {
                Peptide = new mzIdentML110.Generated.PeptideType[peptides.Count],
                DBSequence = new mzIdentML110.Generated.DBSequenceType[proteins.Count],
                PeptideEvidence = new mzIdentML110.Generated.PeptideEvidenceType[peptides.Count]
            };

            _mzid.AnalysisCollection = new mzIdentML110.Generated.AnalysisCollectionType
            {
                SpectrumIdentification = new mzIdentML110.Generated.SpectrumIdentificationType[1]
                {
                    new mzIdentML110.Generated.SpectrumIdentificationType
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
                    DatabaseName = new mzIdentML110.Generated.ParamType
                    {
                        Item = new mzIdentML110.Generated.CVParamType
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
                _mzid.SequenceCollection.DBSequence[protein_index] = new mzIdentML110.Generated.DBSequenceType
                {
                    id = "DBS_" + protein.Accession,
                    lengthSpecified = true,
                    length = protein.Length,
                    searchDatabase_ref = database_reference[protein.DatabaseFilePath],
                    accession = protein.Accession,
                    Seq = protein.BaseSequence,
                    cvParam = new mzIdentML110.Generated.CVParamType[1]
                    {
                        new mzIdentML110.Generated.CVParamType
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
                    FileFormat = new mzIdentML110.Generated.FileFormatType
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType
                        {
                            accession = thermoRawFile ? "MS:1000563" : "MS:1000584",
                            name = thermoRawFile ? "Thermo RAW format" : "mzML format",
                            cvRef = "PSI-MS"
                        }
                    },
                    SpectrumIDFormat = new mzIdentML110.Generated.SpectrumIDFormatType
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType
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
            Dictionary<PeptideWithSetModifications, int> peptide_evidence_ids = new Dictionary<PeptideWithSetModifications, int>();
            Dictionary<string, Tuple<int, HashSet<string>>> peptide_ids = new Dictionary<string, Tuple<int, HashSet<string>>>(); //key is peptide sequence, value is <peptide id for that peptide, peptide evidences>, list of spectra id's
            Dictionary<Tuple<string, int>, Tuple<int, int>> psm_per_scan = new Dictionary<Tuple<string, int>, Tuple<int, int>>(); //key is <filename, scan numer> value is <scan result id, scan item id #'s (could be more than one ID per scan)>

            var unambiguousPsms = items.Where(psm => psm.FullSequence != null);

            foreach (Psm psm in unambiguousPsms)
            {
                foreach (PeptideWithSetModifications peptide in psm.CompactPeptides.SelectMany(c => c.Value.Item2).Distinct())
                {
                    //if first peptide on list hasn't been added, add peptide and peptide evidence
                    if (!peptide_ids.TryGetValue(peptide.Sequence, out Tuple<int, HashSet<string>> peptide_id))
                    {
                        peptide_id = new Tuple<int, HashSet<string>>(p_index, new HashSet<string>());
                        p_index++;
                        _mzid.SequenceCollection.Peptide[peptide_id.Item1] = new mzIdentML110.Generated.PeptideType
                        {
                            PeptideSequence = peptide.BaseSequence,
                            id = "P_" + peptide_id.Item1,
                            Modification = new mzIdentML110.Generated.ModificationType[peptide.NumMods]
                        };
                        int mod_id = 0;
                        foreach (KeyValuePair<int, ModificationWithMass> mod in peptide.allModsOneIsNterminus)
                        {
                            _mzid.SequenceCollection.Peptide[peptide_id.Item1].Modification[mod_id] = new mzIdentML110.Generated.ModificationType()
                            {
                                location = mod.Key - 1,
                                locationSpecified = true,
                                monoisotopicMassDelta = mod.Value.monoisotopicMass,
                                residues = new string[1] { peptide.BaseSequence[Math.Min(Math.Max(0, mod.Key - 2), peptide.Length - 1)].ToString() },
                                monoisotopicMassDeltaSpecified = true,
                                cvParam = new mzIdentML110.Generated.CVParamType[1]
                                {
                                    GetUnimodCvParam(mod.Value)
                                }
                            };
                            mod_id++;
                        }
                        peptide_ids.Add(peptide.Sequence, peptide_id);
                    }

                    if (!peptide_evidence_ids.ContainsKey(peptide))
                    {
                        _mzid.SequenceCollection.PeptideEvidence[pe_index] = new mzIdentML110.Generated.PeptideEvidenceType()
                        {
                            id = "PE_" + pe_index,
                            peptide_ref = "P_" + peptide_id.Item1,
                            dBSequence_ref = "DBS_" + peptide.Protein.Accession,
                            isDecoy = peptide.Protein.IsDecoy,
                            startSpecified = true,
                            start = peptide.OneBasedStartResidueInProtein,
                            endSpecified = true,
                            end = peptide.OneBasedEndResidueInProtein,
                            pre = peptide.PreviousAminoAcid.ToString(),
                            post = (peptide.OneBasedEndResidueInProtein < peptide.Protein.BaseSequence.Length) ? peptide.Protein[peptide.OneBasedEndResidueInProtein].ToString() : "-",
                        };
                        peptide_evidence_ids.Add(peptide, pe_index);
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
                        spectrumID = "scan=" + psm.ScanNumber.ToString(),
                        SpectrumIdentificationItem = new mzIdentML110.Generated.SpectrumIdentificationItemType[500],
                        cvParam = new mzIdentML110.Generated.CVParamType[1]
                        {
                            new mzIdentML110.Generated.CVParamType
                            {
                                name = "scan start time",
                                cvRef = "PSI-MS",
                                accession = "MS:1000016",
                                value = psm.ScanRetentionTime.ToString()
                            }
                        }
                    };
                    psm_per_scan.Add(new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber), scan_result_scan_item);
                    sir_id++;
                }
                else
                {
                    psm_per_scan[new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber)] = new Tuple<int, int>(scan_result_scan_item.Item1, scan_result_scan_item.Item2 + 1);
                    scan_result_scan_item = psm_per_scan[new Tuple<string, int>(psm.FullFilePath, psm.ScanNumber)];
                }
                foreach (PeptideWithSetModifications p in psm.CompactPeptides.SelectMany(c => c.Value.Item2).Distinct())
                {
                    peptide_ids[p.Sequence].Item2.Add("SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2);
                }
                _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2] = new mzIdentML110.Generated.SpectrumIdentificationItemType()
                {
                    rank = 1,
                    chargeState = psm.ScanPrecursorCharge,
                    id = "SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2,
                    experimentalMassToCharge = Math.Round(psm.ScanPrecursorMonoisotopicPeakMz, 5),
                    passThreshold = psm.FdrInfo.QValue <= threshold,
                    //NOTE:ONLY CAN HAVE ONE PEPTIDE REF PER SPECTRUM IDENTIFICATION ITEM
                    peptide_ref = "P_" + peptide_ids[psm.FullSequence].Item1,
                    PeptideEvidenceRef = new mzIdentML110.Generated.PeptideEvidenceRefType[psm.CompactPeptides.SelectMany(c => c.Value.Item2).Distinct().Count()],
                    cvParam = new mzIdentML110.Generated.CVParamType[2]
                    {
                        new mzIdentML110.Generated.CVParamType
                        {
                            name = "MetaMorpheus:score",
                            cvRef = "PSI-MS",
                            accession = "MS:1002827",
                            value = psm.Score.ToString()
                        },
                        new mzIdentML110.Generated.CVParamType
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
                foreach (PeptideWithSetModifications p in psm.CompactPeptides.SelectMany(c => c.Value.Item2).Distinct())
                {
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2].PeptideEvidenceRef[pe]
                        = new mzIdentML110.Generated.PeptideEvidenceRefType
                        {
                            peptideEvidence_ref = "PE_" + peptide_evidence_ids[p]
                        };
                    pe++;
                }
            }

            _mzid.AnalysisProtocolCollection = new mzIdentML110.Generated.AnalysisProtocolCollectionType()
            {
                SpectrumIdentificationProtocol = new mzIdentML110.Generated.SpectrumIdentificationProtocolType[1]
                {
                    new mzIdentML110.Generated.SpectrumIdentificationProtocolType
                    {
                        id = "SIP",
                        analysisSoftware_ref = "AS_MetaMorpheus",
                        SearchType = new mzIdentML110.Generated.ParamType
                        {
                            Item = new mzIdentML110.Generated.CVParamType
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
                                new mzIdentML110.Generated.CVParamType
                                {
                                    accession = "MS:1001211",
                                    cvRef = "PSI-MS",
                                    name = "parent mass type mono"
                                },
                                new mzIdentML110.Generated.CVParamType
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
                            new mzIdentML110.Generated.CVParamType
                            {
                                accession = "MS:1001412",
                                name = "search tolerance plus value",
                                value = productTolerance.Value.ToString(),
                                cvRef = "PSI-MS",
                                unitAccession = productTolerance is PpmTolerance? "UO:0000169": "UO:0000221",
                                unitName = productTolerance is PpmTolerance? "parts per million" : "dalton" ,
                                unitCvRef = "UO"
                            },
                            new mzIdentML110.Generated.CVParamType
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
                       ParentTolerance = new mzIdentML110.Generated.CVParamType[2]
                        {
                              new mzIdentML110.Generated.CVParamType
                            {
                                accession = "MS:1001412",
                                name = "search tolerance plus value",
                                value = parentTolerance.Value.ToString(),
                                cvRef = "PSI-MS",
                                unitAccession = parentTolerance is PpmTolerance? "UO:0000169": "UO:0000221",
                                unitName = parentTolerance is PpmTolerance? "parts per million" : "dalton" ,
                                unitCvRef = "UO"
                            },
                            new mzIdentML110.Generated.CVParamType
                            {
                                accession = "MS:1001413",
                                name = "search tolerance minus value",
                                value = parentTolerance.Value.ToString(),
                                cvRef = "PSI-MS",
                                unitAccession = parentTolerance is PpmTolerance? "UO:0000169": "UO:0000221",
                                unitName = parentTolerance is PpmTolerance? "parts per million" : "dalton" ,
                                unitCvRef = "UO"
                            }
                        },
                        Threshold = new mzIdentML110.Generated.ParamListType()
                        {
                            Items = new mzIdentML110.Generated.CVParamType[1]
                            {
                                new mzIdentML110.Generated.CVParamType
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
                    missedCleavagesSpecified = true,
                    missedCleavages = missedCleavages,
                    SiteRegexp = protease.SiteRegexp,
                    EnzymeName = new mzIdentML110.Generated.ParamListType()
                    {
                        Items = new mzIdentML110.Generated.AbstractParamType[1]
                        {
                            new mzIdentML110.Generated.CVParamType
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
                    cvParam = new mzIdentML110.Generated.CVParamType[1]
                    {
                        GetUnimodCvParam(mod)
                    }
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
                    cvParam = new mzIdentML110.Generated.CVParamType[1]
                    {
                        GetUnimodCvParam(mod)
                    }
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
                        new mzIdentML110.Generated.CVParamType
                        {
                            accession = "MS:1001447",
                            name = "prot:FDR threshold",
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
                foreach (EngineLayer.ProteinGroup proteinGroup in groups)
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
                            new mzIdentML110.Generated.CVParamType
                            {
                                accession = "MS:1002828",
                                name = "MetaMorpheus:protein score",
                                cvRef = "PSI-MS",
                                value = proteinGroup.ProteinGroupScore.ToString()
                            },
                            new mzIdentML110.Generated.CVParamType
                            {
                                accession = "MS:1002373",
                                name = "protein group-level q-value",
                                cvRef = "PSI-MS",
                                value = proteinGroup.QValue.ToString()
                            },
                            new mzIdentML110.Generated.CVParamType
                            {
                                accession = "MS:1001093",
                                name = "sequence coverage",
                                cvRef = "PSI-MS",
                                value = proteinGroup.SequenceCoveragePercent.First().ToString()
                            },
                            new mzIdentML110.Generated.CVParamType
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
                            if (peptide_evidence_ids.ContainsKey(peptide))
                            {
                                if (peptide.Protein == protein)
                                {
                                    _mzid.DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[pag_protein_index].PeptideHypothesis[peptide_id] = new mzIdentML110.Generated.PeptideHypothesisType()
                                    {
                                        peptideEvidence_ref = "PE_" + peptide_evidence_ids[peptide],
                                        SpectrumIdentificationItemRef = new mzIdentML110.Generated.SpectrumIdentificationItemRefType[peptide_ids[peptide.Sequence].Item2.Count],
                                    };

                                    int i = 0;
                                    foreach (string sii in peptide_ids[peptide.Sequence].Item2)
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

        public static MassDiffAcceptor GetMassDiffAcceptor(Tolerance precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    if (precursorMassTolerance is PpmTolerance)
                        return new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value);
                    else
                        return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value);

                case MassDiffAcceptorType.OneMM:
                    return new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, precursorMassTolerance);

                case MassDiffAcceptorType.TwoMM:
                    return new DotMassDiffAcceptor("2mm", new List<double> { 0, 1.0029, 2.0052 }, precursorMassTolerance);

                case MassDiffAcceptorType.ThreeMM:
                    return new DotMassDiffAcceptor("3mm", new List<double> { 0, 1.0029, 2.0052, 3.0077 }, precursorMassTolerance);

                case MassDiffAcceptorType.ModOpen:
                    return new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });

                case MassDiffAcceptorType.Open:
                    return new OpenSearchMode();

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac);

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
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
                        + "\t" + bin.Combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.MedianLength).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.UniprotID);
                }
            }
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this);
            List<Psm> allPsms = new List<Psm>();
            FlashLFQResults flashLfqResults = null;

            Status("Loading modifications...", taskId);

            #region Load modifications

            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = CommonParameters.LocalizeAll ? GlobalVariables.AllModTypesKnown.ToList() : CommonParameters.ListOfModTypesLocalize.ToList();

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

            Status("Loading proteins...", new List<string> { taskId });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, SearchParameters.SearchTarget, SearchParameters.DecoyType, localizeableModificationTypes, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();
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
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + ". ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");
            proseCreatedWhileRunning.Append("report all ambiguity = " + CommonParameters.ReportAllAmbiguity + "; ");

            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxParallelFilesToAnalyze.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxParallelFilesToAnalyze.Value;
            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);

            HashSet<IDigestionParams> ListOfDigestionParams = GetListOfDistinctDigestionParams(CommonParameters, fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b)));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);
            Status("Searching files...", new List<string> { taskId, "Individual Spectra Files" });

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                myFileManager.DoneWithFile(origDataFile);

                var fileSpecificPsms = new Psm[arrayOfMs2ScansSortedByMass.Length];

                if (SearchParameters.SearchType == SearchType.Modern)
                {
                    for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                    {
                        List<CompactPeptide> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                        #region Generate indices for modern search

                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, ionTypes, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, SearchParameters.MaxFragmentSize, new List<string> { taskId });
                        List<int>[] fragmentIndex = null;
                        lock (indexLock)
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, taskId);

                        #endregion Generate indices for modern search

                        Status("Searching files...", taskId);

                        new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, ionTypes, currentPartition, combinedParams, SearchParameters.AddCompIons, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                    }
                }
                else if (SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    List<List<ProductType>> terminusSeparatedIons = ProductTypeMethod.SeparateIonsByTerminus(ionTypes);
                    foreach (List<ProductType> terminusSpecificIons in terminusSeparatedIons)
                    {
                        for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                        {
                            List<CompactPeptide> peptideIndex = null;
                            List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                            List<int>[] fragmentIndex = new List<int>[1];

                            #region Generate indices for nonspecifc search

                            Status("Getting fragment dictionary...", new List<string> { taskId });
                            var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, terminusSpecificIons, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, SearchParameters.MaxFragmentSize, new List<string> { taskId });
                            lock (indexLock)
                                GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, taskId);

                            #endregion Generate indices for nonspecifc search

                            #region Generate indices for nonspecific search

                            Status("Getting precursor dictionary...", new List<string> { taskId });
                            List<CompactPeptide> peptideIndexPrecursor = null;
                            List<Protein> proteinListSubsetPrecursor = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));
                            List<int>[] fragmentIndexPrecursor = new List<int>[1];
                            var indexEnginePrecursor = new PrecursorIndexingEngine(proteinListSubsetPrecursor, variableModifications, fixedModifications, terminusSpecificIons, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, 0, new List<string> { taskId });
                            lock (indexLock)
                                GenerateIndexes(indexEnginePrecursor, dbFilenameList, ref peptideIndexPrecursor, ref fragmentIndexPrecursor, taskId);

                            if (peptideIndex.Count != peptideIndexPrecursor.Count)
                                throw new MetaMorpheusException("peptideIndex not identical between indexing engines");

                            #endregion Generate indices for nonspecific search

                            Status("Searching files...", taskId);

                            new NonSpecificEnzymeSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, fragmentIndexPrecursor, terminusSpecificIons, currentPartition, combinedParams, SearchParameters.AddCompIons, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                            ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                        }
                    }
                }
                else //If classic search
                {
                    Status("Starting search...", thisId);
                    new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ionTypes, massDiffAcceptor, SearchParameters.AddCompIons, combinedParams, combinedParams.ProductMassTolerance, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }
                lock (psmLock)
                {
                    allPsms.AddRange(fileSpecificPsms);
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
                List<List<ProductType>> terminusSeparatedIons = ProductTypeMethod.SeparateIonsByTerminus(ionTypes);
                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
                foreach (List<ProductType> terminusSpecificIons in terminusSeparatedIons)
                    new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, allPsms, proteinList, fixedModifications, variableModifications, terminusSpecificIons, ListOfDigestionParams, massDiffAcceptor, CommonParameters.ReportAllAmbiguity, new List<string> { taskId }).Run();
            }
            else
            {
                SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, ionTypes, ListOfDigestionParams, CommonParameters.ReportAllAmbiguity, new List<string> { taskId });
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;
            }

            ProteinParsimonyResults proteinAnalysisResults = null;
            if (SearchParameters.DoParsimony)
                proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, SearchParameters.ModPeptidesAreDifferent, new List<string> { taskId }).Run());

            Status("Resolving most probable peptide...", new List<string> { taskId });
            foreach (var huh in allPsms)
            {
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            }

            Status("Ordering and grouping psms...", taskId);

            allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            Status("Running FDR analysis...", taskId);
            int massDiffAcceptorNumNotches = GetNumNotches(SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            var fdrAnalysisResults = new FdrAnalysisEngine(allPsms, massDiffAcceptorNumNotches, CommonParameters.CalculateEValue, new List<string> { taskId }).Run();

            List<EngineLayer.ProteinGroup> proteinGroups = null;

            if (SearchParameters.DoParsimony)
            {
                var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, allPsms, SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreDifferent, true, new List<string> { taskId }).Run();
                proteinGroups = proteinScoringAndFdrResults.sortedAndScoredProteinGroups;
            }

            if (SearchParameters.DoLocalizationAnalysis)
            {
                Status("Running localization analysis...", taskId);
                Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    Status("Running localization analysis...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                    new LocalizationEngine(allPsms.Where(b => b.FullFilePath.Equals(origDataFile)).ToList(), ionTypes, myMsDataFile, combinedParams.ProductMassTolerance, new List<string> { taskId, "Individual Spectra Files", origDataFile }, SearchParameters.AddCompIons).Run();
                    myFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                });
            }

            new ModificationAnalysisEngine(allPsms, new List<string> { taskId }).Run();

            if (SearchParameters.DoQuantification)
            {
                // pass quantification parameters to FlashLFQ
                Status("Quantifying...", taskId);

                // construct file info for FlashLFQ
                var rawfileinfos = new List<RawFileInfo>();
                foreach (var file in currentRawFileList)
                {
                    ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[currentRawFileList.IndexOf(file)]);
                    if (myFileManager.SeeIfOpen(file))
                        rawfileinfos.Add(new RawFileInfo(file, myFileManager.LoadFile(file, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks)));
                    else
                        rawfileinfos.Add(new RawFileInfo(file));
                }

                // get PSMs to pass to FlashLFQ
                var unambiguousPsmsBelowOnePercentFdr = allPsms.Where(p => p.FdrInfo.QValue < 0.01 && !p.IsDecoy && p.FullSequence != null);
                var mypsmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

                // pass protein group info for each PSM
                Dictionary<Psm, List<string>> psmToProteinGroupNames = new Dictionary<Psm, List<string>>();
                if (proteinGroups != null)
                {
                    EngineLayer.ProteinGroup.FilesForQuantification = rawfileinfos.Select(p => p.fullFilePathWithExtension).ToArray();

                    foreach (var proteinGroup in proteinGroups)
                    {
                        foreach (var psm in proteinGroup.AllPsmsBelowOnePercentFDR)
                        {
                            if (psmToProteinGroupNames.TryGetValue(psm, out List<string> proteinGroupNames))
                                proteinGroupNames.Add(proteinGroup.ProteinGroupName);
                            else
                                psmToProteinGroupNames.Add(psm, new List<string> { proteinGroup.ProteinGroupName });
                        }
                    }
                }
                else
                {
                    // if protein groups were not constructed, just use accession numbers
                    foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                    {
                        var proteinsAccessionString = psm.CompactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.Accession).Distinct();
                        psmToProteinGroupNames.Add(psm, proteinsAccessionString.ToList());
                    }
                }

                // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
                // the peptides should still be quantified but not considered for protein quantification
                foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                    if (!psmToProteinGroupNames.ContainsKey(psm))
                        psmToProteinGroupNames.Add(psm, new List<string>() { "" });

                // pass PSM info to FlashLFQ
                var flashLFQIdentifications = new List<Identification>();
                foreach (var file in mypsmsGroupedByFile)
                {
                    var rawfileinfo = rawfileinfos.Where(p => p.fullFilePathWithExtension.Equals(file.Key)).First();
                    foreach (var psm in file)
                        flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence, (double)psm.PeptideMonisotopicMass, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroupNames[psm]));
                }

                // run FlashLFQ
                var FlashLfqEngine = new FlashLFQEngine(flashLFQIdentifications, SearchParameters.QuantifyPpmTol, 5.0, SearchParameters.MatchBetweenRuns, 5.0, false, 2, false, true, true, GlobalVariables.ElementsLocation);
                flashLfqResults = FlashLfqEngine.Run();

                // get protein intensity back from FlashLFQ
                if (proteinGroups != null)
                {
                    Dictionary<string, EngineLayer.ProteinGroup> proteinGroupNameToProteinGroup = new Dictionary<string, EngineLayer.ProteinGroup>();
                    foreach (var proteinGroup in proteinGroups)
                    {
                        proteinGroup.IntensitiesByFile = new double[EngineLayer.ProteinGroup.FilesForQuantification.Length];
                        if (!proteinGroupNameToProteinGroup.ContainsKey(proteinGroup.ProteinGroupName))
                            proteinGroupNameToProteinGroup.Add(proteinGroup.ProteinGroupName, proteinGroup);
                    }

                    foreach (var flashLfqProteinGroup in flashLfqResults.proteinGroups)
                    {
                        if (proteinGroupNameToProteinGroup.TryGetValue(flashLfqProteinGroup.Key, out EngineLayer.ProteinGroup metamorpheusProteinGroup))
                        {
                            for (int i = 0; i < EngineLayer.ProteinGroup.FilesForQuantification.Length; i++)
                                metamorpheusProteinGroup.IntensitiesByFile[i] = flashLfqProteinGroup.Value.intensities[rawfileinfos[i]];
                        }
                    }
                }
            }

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            if (SearchParameters.DoHistogramAnalysis)
            {
                var limitedpsms_with_fdr = allPsms.Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                {
                    Status("Running histogram analysis...", new List<string> { taskId });
                    var myTreeStructure = new BinTreeStructure();
                    myTreeStructure.GenerateBins(limitedpsms_with_fdr, SearchParameters.HistogramBinTolInDaltons);
                    var writtenFile = Path.Combine(OutputFolder, "aggregate.mytsv");
                    WriteTree(myTreeStructure, writtenFile);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });
                }
            }

            Status("Writing results...", taskId);
            {
                if (currentRawFileList.Count > 1)
                {
                    var writtenFile = Path.Combine(OutputFolder, "aggregatePSMs.psmtsv");
                    WritePsmsToTsv(allPsms, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });

                    var writtenFileForPercolator = Path.Combine(OutputFolder, "forPercolator.tsv");
                    WritePsmsForPercolator(allPsms, writtenFileForPercolator);
                    SucessfullyFinishedWritingFile(writtenFileForPercolator, new List<string> { taskId });
                }
                myTaskResults.AddNiceText("All target PSMS within 1% FDR: " + allPsms.Count(a => a.FdrInfo.QValue < .01 && !a.IsDecoy));
            }

            var peptides = allPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
            {
                if (currentRawFileList.Count > 1)
                {
                    var writtenFile = Path.Combine(OutputFolder, "aggregatePeptides.psmtsv");
                    WritePsmsToTsv(peptides, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });
                }
                myTaskResults.AddNiceText("Target peptides within 1% FDR: " + peptides.Count(a => a.FdrInfo.QValue < .01 && !a.IsDecoy));
            }

            if (SearchParameters.DoParsimony)
                myTaskResults.AddNiceText("Target protein groups within 1% FDR: " + proteinGroups.Count(b => b.QValue < 0.01 && !b.isDecoy) + Environment.NewLine);

            var psmsGroupedByFile = allPsms.GroupBy(p => p.FullFilePath);

            // individual psm files (with global psm fdr, global parsimony)
            foreach (var group in psmsGroupedByFile)
            {
                var psmsForThisFile = group.ToList();

                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                {
                    var writtenFile = Path.Combine(OutputFolder, strippedFileName + "_PSMs.psmtsv");
                    WritePsmsToTsv(psmsForThisFile, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                    myTaskResults.AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue < .01 && a.IsDecoy == false));
                }

                var writtenFileForPercolator = Path.Combine(OutputFolder, strippedFileName + "_forPercolator.tsv");
                WritePsmsForPercolator(psmsForThisFile, writtenFileForPercolator);
                SucessfullyFinishedWritingFile(writtenFileForPercolator, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
            }
            foreach (var group in psmsGroupedByFile)
            {
                var psmsForThisFile = group.ToList();

                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                {
                    var peptidesForFile = psmsForThisFile.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                    var writtenFile = Path.Combine(OutputFolder, strippedFileName + "_Peptides.psmtsv");
                    WritePsmsToTsv(peptidesForFile, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                    myTaskResults.AddNiceText("Target peptides within 1% FDR in " + strippedFileName + ": " + peptidesForFile.Count(a => a.FdrInfo.QValue < .01 && a.IsDecoy == false));
                }
            }

            if (SearchParameters.DoParsimony)
            {
                if (currentRawFileList.Count > 1)
                    WriteProteinGroupsToTsv(proteinGroups, OutputFolder, "aggregateProteinGroups", new List<string> { taskId }, psmsGroupedByFile.Select(b => b.Key).ToList());

                // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                foreach (var fullFilePath in currentRawFileList)
                {
                    List<Psm> psmsForThisFile = psmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();

                    var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                    var subsetProteinGroupsForThisFile = new List<EngineLayer.ProteinGroup>();
                    foreach (var pg in proteinGroups)
                        subsetProteinGroupsForThisFile.Add(pg.ConstructSubsetProteinGroup(fullFilePath));

                    var subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile, SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreDifferent, false, new List<string> { taskId, "Individual Spectra Files", fullFilePath }).Run();
                    subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.sortedAndScoredProteinGroups;

                    myTaskResults.AddNiceText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.Count(b => b.QValue < 0.01 && !b.isDecoy));

                    WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, OutputFolder, strippedFileName + "_ProteinGroups", new List<string> { taskId, "Individual Spectra Files", fullFilePath }, new List<string> { fullFilePath });

                    Status("Writing mzid...", new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                    var mzidFilePath = Path.Combine(OutputFolder, strippedFileName + ".mzid");
                    WriteMzidentml(psmsForThisFile, subsetProteinGroupsForThisFile, variableModifications, fixedModifications, new List<Protease> { CommonParameters.DigestionParams.Protease }, 0.01, CommonParameters.ProductMassTolerance, CommonParameters.PrecursorMassTolerance, CommonParameters.DigestionParams.MaxMissedCleavages, mzidFilePath);
                    SucessfullyFinishedWritingFile(mzidFilePath, new List<string> { taskId, "Individual Spectra Files", fullFilePath });

                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", fullFilePath }));
                }
            }

            if (SearchParameters.DoQuantification)
            {
                foreach (var file in flashLfqResults.peaks)
                    WritePeakQuantificationResultsToTsv(file.Value, OutputFolder, file.Key.filenameWithoutExtension + "_QuantifiedPeaks", new List<string> { taskId, "Individual Spectra Files", file.Key.fullFilePathWithExtension });

                if (currentRawFileList.Count > 1)
                    WritePeakQuantificationResultsToTsv(flashLfqResults.peaks.SelectMany(p => p.Value).ToList(), OutputFolder, "aggregateQuantifiedPeaks", new List<string> { taskId });

                WritePeptideQuantificationResultsToTsv(flashLfqResults.peptideBaseSequences.Select(p => p.Value).OrderBy(p => p.Sequence).ToList(), OutputFolder, "aggregateQuantifiedPeptidesByBaseSeq", new List<string> { taskId });
                WritePeptideQuantificationResultsToTsv(flashLfqResults.peptideModifiedSequences.Select(p => p.Value).OrderBy(p => p.Sequence).ToList(), OutputFolder, "aggregateQuantifiedPeptidesByFullSeq", new List<string> { taskId });
            }

            if (SearchParameters.WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { taskId });
                List<Modification> modificationsToWriteIfBoth = new List<Modification>();
                List<Modification> modificationsToWriteIfInDatabase = new List<Modification>();
                List<Modification> modificationsToWriteIfObserved = new List<Modification>();

                var goodPsmsForEachProtein = allPsms.Where(b => b.FdrInfo.QValueNotch < 0.01 && !b.IsDecoy && b.FullSequence != null && b.ProteinAccesion != null).GroupBy(b => b.CompactPeptides.First().Value.Item2.First().Protein).ToDictionary(b => b.Key);
                // Add user mod selection behavours to Pruned DB
                foreach (var modType in SearchParameters.ModsToWriteSelection)
                {
                    if (modType.Value == 1) // Write if observed and in database
                        modificationsToWriteIfBoth.AddRange(GlobalVariables.AllModsKnown.Where(b => b.modificationType.Equals(modType.Key)));
                    if (modType.Value == 2) // Write if in database
                        modificationsToWriteIfInDatabase.AddRange(GlobalVariables.AllModsKnown.Where(b => b.modificationType.Equals(modType.Key)));
                    if (modType.Value == 3) // Write if observed
                        modificationsToWriteIfObserved.AddRange(GlobalVariables.AllModsKnown.Where(b => b.modificationType.Equals(modType.Key)));
                }

                foreach (var protein in proteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        HashSet<Tuple<int, ModificationWithMass>> modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>();
                        if (goodPsmsForEachProtein.ContainsKey(protein))
                            modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>(goodPsmsForEachProtein[protein].SelectMany(b => b.CompactPeptides.First().Value.Item2.First().allModsOneIsNterminus.Select(c => new Tuple<int, ModificationWithMass>(GetOneBasedIndexInProtein(c.Key, b.CompactPeptides.First().Value.Item2.First()), c.Value))));

                        IDictionary<int, List<Modification>> modsToWrite = new Dictionary<int, List<Modification>>();

                        foreach (var observedMod in modsObservedOnThisProtein)
                        {
                            //Add if observed (regardless if in database)
                            var tempMod = observedMod.Item2;

                            if (modificationsToWriteIfObserved.Contains(tempMod as Modification))
                            {
                                if (!modsToWrite.ContainsKey(observedMod.Item1))
                                    modsToWrite.Add(observedMod.Item1, new List<Modification> { observedMod.Item2 as Modification });
                                else
                                    modsToWrite[observedMod.Item1].Add(observedMod.Item2 as Modification);
                                continue;
                            }
                        }

                        // Add if in database (two cases: always or if observed)
                        foreach (var modd in protein.OneBasedPossibleLocalizedModifications)
                            foreach (var mod in modd.Value)
                            {
                                //Add if always In Database or if was observed and in database and not set to not include
                                if (modificationsToWriteIfInDatabase.Contains(mod as Modification) ||
                                (modsObservedOnThisProtein.Contains(new Tuple<int, ModificationWithMass>(modd.Key, mod as ModificationWithMass)) && modificationsToWriteIfBoth.Contains(mod as Modification)))
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

        private static mzIdentML110.Generated.CVParamType GetUnimodCvParam(ModificationWithMass mod)
        {
            if (mod.linksToOtherDbs.ContainsKey("Unimod"))
            {
                return new mzIdentML110.Generated.CVParamType()
                {
                    accession = "UNIMOD:" + mod.linksToOtherDbs["Unimod"].First(),
                    name = mod.id,
                    cvRef = "PSI-MS",
                };
            }
            else
                return new mzIdentML110.Generated.CVParamType()
                {
                    accession = "MS:1001460",
                    name = "unknown modification",
                    cvRef = "UNIMOD",
                    value = mod.id,
                };
        }

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor ye = null;

            var split = text.Split(' ');

            switch (split[1])
            {
                case "dot":

                    var massShifts = Array.ConvertAll(split[4].Split(','), Double.Parse);
                    var newString = split[2].Replace("�", "");
                    var toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                    if (split[3].ToUpperInvariant().Equals("PPM"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                    else if (split[3].ToUpperInvariant().Equals("DA"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
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
                    throw new MetaMorpheusException("Could not parse search mode string");
            }
            return ye;
        }

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

        private static void WriteFragmentIndexNetSerializer(List<int>[] fragmentIndex, string fragmentIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
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

        private void WritePsmsForPercolator(List<Psm> psmList, string writtenFileForPercolator)
        {
            using (StreamWriter output = new StreamWriter(writtenFileForPercolator))
            {
                output.WriteLine("SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins");
                output.WriteLine("DefaultDirection\t-\t-\t1\t1\t\t");
                for (int i = 0; i < psmList.Count; i++)
                {
                    var heh = psmList[i];

                    output.Write(i.ToString());
                    output.Write('\t' + (heh.IsDecoy ? -1 : 1).ToString());
                    output.Write('\t' + heh.ScanNumber.ToString());

                    // Features
                    {
                        output.Write('\t' + string.Join("\t", heh.Features));
                    }

                    // HACKY: Ignores all ambiguity
                    var pwsm = heh.CompactPeptides.First().Value.Item2.First();

                    output.Write('\t' + (pwsm.PreviousAminoAcid + "." + pwsm.Sequence + "." + pwsm.NextAminoAcid).ToString());
                    output.Write('\t' + (pwsm.Protein.Accession).ToString());
                    output.WriteLine();
                }
            }
        }

        private int GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    return 1;

                case MassDiffAcceptorType.OneMM:
                    return 2;

                case MassDiffAcceptorType.TwoMM:
                    return 3;

                case MassDiffAcceptorType.ThreeMM:
                    return 4;

                case MassDiffAcceptorType.ModOpen:
                    return 1;

                case MassDiffAcceptorType.Open:
                    return 1;

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac).NumNotches;

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        private void GenerateIndexes(IndexingEngine indexEngine, List<DbForTask> dbFilenameList, ref List<CompactPeptide> peptideIndex, ref List<int>[] fragmentIndex, string taskId)
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
                fragmentIndex = indexResults.FragmentIndex;

                Status("Writing peptide index...", new List<string> { taskId });
                var peptideIndexFile = Path.Combine(output_folderForIndices, "peptideIndex.ind");
                WritePeptideIndex(peptideIndex, peptideIndexFile);
                SucessfullyFinishedWritingFile(peptideIndexFile, new List<string> { taskId });

                Status("Writing fragment index...", new List<string> { taskId });
                var fragmentIndexFile = Path.Combine(output_folderForIndices, "fragmentIndex.ind");
                WriteFragmentIndexNetSerializer(fragmentIndex, fragmentIndexFile);
                SucessfullyFinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
            }
            else
            {
                Status("Reading peptide index...", new List<string> { taskId });
                var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
                var ser = new NetSerializer.Serializer(messageTypes);
                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
                {
                    peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);
                }

                Status("Reading fragment index...", new List<string> { taskId });
                messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
                ser = new NetSerializer.Serializer(messageTypes);
                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
                    fragmentIndex = (List<int>[])ser.Deserialize(file);
            }
        }

        private void WriteProteinGroupsToTsv(List<EngineLayer.ProteinGroup> items, string outputFolder, string strippedFileName, List<string> nestedIds, List<string> FileNames)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, strippedFileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(EngineLayer.ProteinGroup.GetTabSeparatedHeader(FileNames.Count == 1));
                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        private void WritePeptideQuantificationResultsToTsv(List<FlashLFQ.Peptide> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(FlashLFQ.Peptide.TabSeparatedHeader);

                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        private void WritePeakQuantificationResultsToTsv(List<FlashLFQ.ChromatographicPeak> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(FlashLFQ.ChromatographicPeak.TabSeparatedHeader);

                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }
                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        #endregion Private Methods
    }
}