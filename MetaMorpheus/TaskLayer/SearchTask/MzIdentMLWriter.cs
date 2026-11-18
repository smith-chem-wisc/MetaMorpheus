using Chemistry;
using EngineLayer;
using MzLibUtil;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Xml;
using System.Xml.Serialization;
using Omics.Digestion;
using Omics.Modifications;

namespace TaskLayer
{
    public static class MzIdentMLWriter
    {
        public static void WriteMzIdentMl(IEnumerable<SpectralMatch> psms, List<EngineLayer.ProteinGroup> groups, List<Modification> variableMods, 
            List<Modification> fixedMods, List<SilacLabel> silacLabels, List<DigestionAgent> proteases, Tolerance productTolerance, 
            Tolerance parentTolerance, int missedCleavages, string outputPath, bool appendMotifToModNames, bool? semiSpecific = null)
        {

            //if SILAC, remove the silac labels, because the base/full sequences reported for output are not the same as the peptides in the best peptides list for the psm
            if (silacLabels != null)
            {
                List<string> labelsToSearch = new List<string>();
                foreach (SilacLabel label in silacLabels)
                {
                    labelsToSearch.Add(label.MassDifference);
                    if (label.AdditionalLabels != null)
                    {
                        foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                        {
                            labelsToSearch.Add(additionalLabel.MassDifference);
                        }
                    }
                }
                // SpectralMatch.FullSequence is assigned from PsmTsvWriter.Resolve(...).ResolvedValue,
                // which is null whenever the candidate sequences disagree -- the pipe-joined form is
                // ResolvedString and is discarded. So FullSequence is either a clean sequence or null,
                // never "|"-delimited, and the old !FullSequence.Contains("|") conjunct could only ever
                // be true or throw a NullReferenceException on an ambiguous PSM. Dropping it fixes the
                // crash without narrowing psms, which also feeds peptides/proteins/filenames below;
                // ambiguous PSMs are excluded from the identification data by unambiguousPsms instead.
                psms = psms.Where(p => p.BaseSequence != null && !labelsToSearch.Any(x => p.BaseSequence.Contains(x)));
            }

            List<PeptideWithSetModifications> peptides = psms.SelectMany(i => i.BestMatchingBioPolymersWithSetMods.Select(v => v.SpecificBioPolymer as PeptideWithSetModifications)).Distinct().ToList();
            List<Protein> proteins = peptides.Select(p => p.Protein).Distinct().ToList();
            List<string> filenames = psms.Select(i => i.FullFilePath).Distinct().ToList();
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
                fullName = "UNIMOD",
                uri = "http://www.unimod.org/obo/unimod.obo"
            },
              new mzIdentML110.Generated.cvType()
            {
                id = "UO",
                fullName = "UNIT-ONTOLOGY",
                uri = "http://purl.obolibrary.org/obo/uo.obo"
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
                            SpectrumIdentificationResult = new mzIdentML110.Generated.SpectrumIdentificationResultType[psms.Count()]
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
            foreach (string database in databases.Where(p=>p!=null))
            {
                _mzid.DataCollection.Inputs.SearchDatabase[database_index] = new mzIdentML110.Generated.SearchDatabaseType()
                {
                    id = "SDB_" + database_index,
                    // The file name, not the local path: a path names the searcher's own folders (often their user
                    // name), and a repository such as PRIDE matches the file by its name.
                    location = Path.GetFileName(database),
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
            foreach (Protein protein in proteins.Where(p=>p.DatabaseFilePath!=null))
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
            Dictionary<string, SpectraFormat> spectra_formats = new Dictionary<string, SpectraFormat>(); //key is datafile
            int spectra_data_id = 0;
            foreach (string data_filepath in filenames)
            {
                SpectraFormat format = SpectraFormat.Of(data_filepath);
                spectra_formats.Add(data_filepath, format);
                string spectral_data_id = "SD_" + spectra_data_id;
                spectral_ids.Add(data_filepath, spectra_data_id);
                _mzid.AnalysisCollection.SpectrumIdentification[0].InputSpectra[spectra_data_id] = new mzIdentML110.Generated.InputSpectraType()
                {
                    spectraData_ref = spectral_data_id
                };
                // a Bruker .d is a folder, which a trailing separator would leave with no file name
                string fileName = Path.GetFileName(Path.TrimEndingDirectorySeparator(data_filepath));
                _mzid.DataCollection.Inputs.SpectraData[spectra_data_id] = new mzIdentML110.Generated.SpectraDataType()
                {
                    id = spectral_data_id,
                    name = Path.GetFileNameWithoutExtension(fileName),
                    location = fileName, // the file name, as for SearchDatabase above
                    FileFormat = new mzIdentML110.Generated.FileFormatType
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType
                        {
                            accession = format.FileFormatAccession,
                            name = format.FileFormatName,
                            cvRef = "PSI-MS"
                        }
                    },
                    SpectrumIDFormat = new mzIdentML110.Generated.SpectrumIDFormatType
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType
                        {
                            accession = format.IdFormatAccession,
                            name = format.IdFormatName,
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

            var unambiguousPsms = psms.Where(psm => psm.FullSequence != null);

            foreach (SpectralMatch psm in unambiguousPsms)
            {
                foreach (PeptideWithSetModifications peptide in psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer).Distinct())
                {
                    //if first peptide on list hasn't been added, add peptide and peptide evidence
                    if (!peptide_ids.TryGetValue(peptide.FullSequence, out Tuple<int, HashSet<string>> peptide_id))
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
                        // in position order: the dictionary's own order depends on how the peptide was built (digested
                        // here, or read back from the Modern search's peptide-index cache), and the file should not
                        foreach (KeyValuePair<int, Modification> mod in peptide.AllModsOneIsNterminus.OrderBy(kv => kv.Key))
                        {
                            _mzid.SequenceCollection.Peptide[peptide_id.Item1].Modification[mod_id] = new mzIdentML110.Generated.ModificationType()
                            {
                                location = mod.Key - 1,
                                locationSpecified = true,
                                monoisotopicMassDelta = mod.Value.MonoisotopicMass.Value,
                                residues = new string[1] { peptide.BaseSequence[Math.Min(Math.Max(0, mod.Key - 2), peptide.Length - 1)].ToString() },
                                monoisotopicMassDeltaSpecified = true,
                                cvParam = new mzIdentML110.Generated.CVParamType[1]
                                {
                                    GetUnimodCvParam(mod.Value, appendMotifToModNames)
                                }
                            };
                            mod_id++;
                        }
                        peptide_ids.Add(peptide.FullSequence, peptide_id);
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
                            start = peptide.OneBasedStartResidue,
                            endSpecified = true,
                            end = peptide.OneBasedEndResidue,
                            pre = peptide.PreviousResidue.ToString(),
                            post = (peptide.OneBasedEndResidue < peptide.Protein.BaseSequence.Length) ? peptide.Protein[peptide.OneBasedEndResidue].ToString() : "-",
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
                        spectrumID = spectra_formats[psm.FullFilePath].SpectrumIdOf(psm),
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
                foreach (PeptideWithSetModifications p in psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer).Distinct())
                {
                    peptide_ids[p.FullSequence].Item2.Add("SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2);
                }
                _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2] = new mzIdentML110.Generated.SpectrumIdentificationItemType()
                {
                    rank = 1,
                    chargeState = psm.ScanPrecursorCharge,
                    id = "SII_" + scan_result_scan_item.Item1 + "_" + scan_result_scan_item.Item2,
                    experimentalMassToCharge = Math.Round(psm.ScanPrecursorMonoisotopicPeakMz, 5),
                    passThreshold = psm.FdrInfo.QValue <= 0.01,
                    //NOTE:ONLY CAN HAVE ONE PEPTIDE REF PER SPECTRUM IDENTIFICATION ITEM
                    peptide_ref = "P_" + peptide_ids[psm.FullSequence].Item1,
                    PeptideEvidenceRef = new mzIdentML110.Generated.PeptideEvidenceRefType[psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer).Distinct().Count()],
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
                if (psm.BioPolymerWithSetModsMonoisotopicMass.HasValue)
                {
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2].calculatedMassToCharge = Math.Round(psm.BioPolymerWithSetModsMonoisotopicMass.Value.ToMz(psm.ScanPrecursorCharge), 5);
                    _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item.Item1].SpectrumIdentificationItem[scan_result_scan_item.Item2].calculatedMassToChargeSpecified = true;
                }

                int pe = 0;
                foreach (PeptideWithSetModifications p in psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer).Distinct())
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
                                    value = "0.01"
                                }
                            }
                        }
                    }
                }
            };

            int protease_index = 0;
            foreach (DigestionAgent protease in proteases)
            {
                // The search's own specificity when the caller gives it: trypsin is Full even in a semi-specific search,
                // whose semi-specificity comes from SearchModeType. Otherwise what the protease itself says.
                bool isSemiSpecific = semiSpecific ?? protease.CleavageSpecificity == CleavageSpecificity.Semi;
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].Enzymes.Enzyme[protease_index] = new mzIdentML110.Generated.EnzymeType()
                {
                    id = "E_" + protease_index,
                    name = protease.Name,
                    semiSpecific = isSemiSpecific,
                    // Written only for a semi-specific search. The attribute is optional, and every other search's file stays
                    // as it always was, with no semiSpecific attribute at all.
                    semiSpecificSpecified = isSemiSpecific,
                    missedCleavagesSpecified = true,
                    missedCleavages = missedCleavages,
                    EnzymeName = new mzIdentML110.Generated.ParamListType()
                    {
                        Items = new mzIdentML110.Generated.AbstractParamType[1]
                        {
                            new mzIdentML110.Generated.CVParamType
                            {
                                accession = protease is Protease prot ? prot.PsiMsAccessionNumber : protease.Name,
                                name = protease is Protease prot2 ? prot2.PsiMsName: protease.Name,
                                cvRef = "PSI-MS"
                            }
                        }
                    }
                };
                protease_index++;
            }

            int mod_index = 0;
            foreach (Modification mod in fixedMods)
            {
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = new mzIdentML110.Generated.SearchModificationType()
                {
                    fixedMod = true,
                    massDelta = (float)mod.MonoisotopicMass,
                    residues = mod.Target.ToString(),
                    cvParam = new mzIdentML110.Generated.CVParamType[1]
                    {
                        GetUnimodCvParam(mod, appendMotifToModNames)
                    }
                };
                mod_index++;
            }

            foreach (Modification mod in variableMods)
            {
                _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = new mzIdentML110.Generated.SearchModificationType()
                {
                    fixedMod = false,
                    massDelta = (float)mod.MonoisotopicMass,
                    residues = mod.Target.ToString(),
                    cvParam = new mzIdentML110.Generated.CVParamType[1]
                    {
                        GetUnimodCvParam(mod, appendMotifToModNames)
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
                            value = "0.01"
                        }
                    }
                }
            };

            if (groups != null)
            {
                // Built up in lists rather than fixed-size arrays because a ProteinDetectionHypothesis
                // with no PeptideHypothesis, and a ProteinAmbiguityGroup with no
                // ProteinDetectionHypothesis, are both schema-invalid and have to be left out.
                var ambiguityGroups = new List<mzIdentML110.Generated.ProteinAmbiguityGroupType>();

                int group_id = 0;
                int protein_id = 0;
                foreach (EngineLayer.ProteinGroup proteinGroup in groups)
                {
                    var detectionHypotheses = new List<mzIdentML110.Generated.ProteinDetectionHypothesisType>();
                    foreach (Protein protein in proteinGroup.Proteins)
                    {
                        var peptideHypotheses = new List<mzIdentML110.Generated.PeptideHypothesisType>();
                        foreach (PeptideWithSetModifications peptide in proteinGroup.AllPeptides)
                        {
                            if (peptide_evidence_ids.ContainsKey(peptide))
                            {
                                if (peptide.Protein == protein)
                                {
                                    var peptideHypothesis = new mzIdentML110.Generated.PeptideHypothesisType()
                                    {
                                        peptideEvidence_ref = "PE_" + peptide_evidence_ids[peptide],
                                        SpectrumIdentificationItemRef = new mzIdentML110.Generated.SpectrumIdentificationItemRefType[peptide_ids[peptide.FullSequence].Item2.Count],
                                    };

                                    int i = 0;
                                    foreach (string sii in peptide_ids[peptide.FullSequence].Item2)
                                    {
                                        peptideHypothesis.SpectrumIdentificationItemRef[i] = new mzIdentML110.Generated.SpectrumIdentificationItemRefType()
                                        {
                                            spectrumIdentificationItem_ref = sii
                                        };
                                        i++;
                                    }
                                    peptideHypotheses.Add(peptideHypothesis);
                                }
                            }
                        }

                        int this_protein_id = protein_id;
                        protein_id++;

                        // No peptide evidence written for this protein, e.g. a group supported only
                        // by PSMs whose sequence could not be resolved.
                        if (peptideHypotheses.Count == 0)
                        {
                            continue;
                        }

                        detectionHypotheses.Add(new mzIdentML110.Generated.ProteinDetectionHypothesisType()
                        {
                            id = "PDH_" + this_protein_id,
                            dBSequence_ref = "DBS_" + protein.Accession,
                            passThreshold = proteinGroup.QValue <= 0.01, // hardcoded as 1% FDR but we could change this to the provided threshold
                            PeptideHypothesis = peptideHypotheses.ToArray(),
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
                                value = proteinGroup.SequenceCoverageFraction.First().ToString()
                            },
                            new mzIdentML110.Generated.CVParamType
                            {
                                accession = "MS:1001097",
                                name = "distinct peptide sequences",
                                cvRef = "PSI-MS",
                                value = proteinGroup.NumberOfUniquePeptides.ToString()
                            }
                            }
                        });
                    }

                    int this_group_id = group_id;
                    group_id++;

                    if (detectionHypotheses.Count == 0)
                    {
                        continue;
                    }

                    ambiguityGroups.Add(new mzIdentML110.Generated.ProteinAmbiguityGroupType()
                    {
                        id = "PAG_" + this_group_id,
                        ProteinDetectionHypothesis = detectionHypotheses.ToArray()
                    });
                }

                _mzid.DataCollection.AnalysisData.ProteinDetectionList = new mzIdentML110.Generated.ProteinDetectionListType()
                {
                    id = "PDL",
                    ProteinAmbiguityGroup = ambiguityGroups.ToArray()
                };
            }
            XmlWriter writer = XmlWriter.Create(outputPath, settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        // mzIdentML 1.1 takes a modification's term from UNIMOD or PSI-MOD, or else MS:1001460 "unknown modification".
        // Each term carries its own vocabulary's cvRef, declared in the cvList above. RESID is not one of mzIdentML's
        // vocabularies, so a RESID reference is never written: UniProt's modifications that have one also have a PSI-MOD
        // reference, and anything else is an unknown modification under its own name.
        private static mzIdentML110.Generated.CVParamType GetUnimodCvParam(Modification mod, bool appendMotifToModNames)
        {
            string modName = appendMotifToModNames ? mod.IdWithMotif : mod.OriginalId;

            if (mod.DatabaseReference != null && mod.DatabaseReference.TryGetValue("Unimod", out var unimod) && unimod.Count > 0)
            {
                return new mzIdentML110.Generated.CVParamType()
                {
                    accession = "UNIMOD:" + unimod.First(),
                    name = modName,
                    cvRef = "UNIMOD",
                };
            }

            if (mod.DatabaseReference != null && mod.DatabaseReference.TryGetValue("PSI-MOD", out var psiMod) && psiMod.Count > 0)
            {
                // UniProt's ptmlist already writes the accession as "MOD:00046"
                string accession = psiMod.First();
                return new mzIdentML110.Generated.CVParamType()
                {
                    accession = accession.StartsWith("MOD:", StringComparison.Ordinal) ? accession : "MOD:" + accession,
                    name = modName,
                    cvRef = "PSI-MOD",
                };
            }

            return new mzIdentML110.Generated.CVParamType()
            {
                accession = "MS:1001460",
                name = "unknown modification",
                cvRef = "PSI-MS",
                value = modName,
            };
        }

        /// <summary>
        /// A spectra file's format and the format of the spectrum IDs written for it, in PSI-MS terms. mzLib's readers
        /// name their formats too (MsDataFile.GetSourceFile), but not always in PSI-MS's words: its BAF reader says
        /// "mzML format" and its timsTOF reader ".D format".
        /// </summary>
        private sealed class SpectraFormat
        {
            public string FileFormatAccession { get; }
            public string FileFormatName { get; }
            public string IdFormatAccession { get; }
            public string IdFormatName { get; }

            // mzLib gives an MGF spectrum no native ID, so it is named by its scan number instead
            private readonly bool _idIsScanNumber;

            private SpectraFormat(string fileFormatAccession, string fileFormatName, string idFormatAccession, string idFormatName,
                bool idIsScanNumber = false)
            {
                FileFormatAccession = fileFormatAccession;
                FileFormatName = fileFormatName;
                IdFormatAccession = idFormatAccession;
                IdFormatName = idFormatName;
                _idIsScanNumber = idIsScanNumber;
            }

            public string SpectrumIdOf(SpectralMatch psm) => _idIsScanNumber ? "scan=" + psm.ScanNumber : psm.NativeId;

            public static SpectraFormat Of(string spectraFilePath)
            {
                switch (Path.GetExtension(Path.TrimEndingDirectorySeparator(spectraFilePath)).ToLowerInvariant())
                {
                    case ".raw":
                        return new SpectraFormat("MS:1000563", "Thermo RAW format", "MS:1000768", "Thermo nativeID format");

                    // The scan number MetaMorpheus reports for an MGF spectrum: its SCANS value, or else its one-based
                    // position in the file. PSI's own MGF IDs, index=<zero-based position>, need the position, which only
                    // mzLib's reader knows.
                    case ".mgf":
                        return new SpectraFormat("MS:1001062", "Mascot MGF format", "MS:1000776", "scan number only nativeID format",
                            idIsScanNumber: true);

                    // mzLib names an msalign spectrum scan=<its SCANS value>
                    case ".msalign":
                        return new SpectraFormat("MS:1002899", "msalign format", "MS:1000776", "scan number only nativeID format");

                    case ".d":
                        return OfBrukerFolder(spectraFilePath);

                    // mzML, including the -calib and -averaged files MetaMorpheus writes itself
                    default:
                        return new SpectraFormat("MS:1000584", "mzML format", "MS:1001530", "mzML unique identifier");
                }
            }

            // Which Bruker format a .d holds is decided by the files inside it, in the order mzLib's readers decide it:
            // analysis.baf (BrukerFileReader), then analysis.tdf before analysis.tsf (TimsTofFileReader).
            private static SpectraFormat OfBrukerFolder(string folder)
            {
                if (File.Exists(Path.Combine(folder, "analysis.baf")))
                {
                    // mzLib's IDs: scan=N
                    return new SpectraFormat("MS:1000815", "Bruker BAF format", "MS:1000772", "Bruker BAF nativeID format");
                }

                if (!File.Exists(Path.Combine(folder, "analysis.tdf")) && File.Exists(Path.Combine(folder, "analysis.tsf")))
                {
                    // mzLib's IDs: frame=N
                    return new SpectraFormat("MS:1003282", "Bruker TSF format", "MS:1003283", "Bruker TSF nativeID format");
                }

                // timsTOF. mzLib's spectra here merge many: a PASEF precursor is frames=A-B;scans=C-D. PSI's TDF IDs are
                // frame=F scan=S, and its combined-spectra form (MS:1002819) is a comma-separated list of those, so neither
                // describes mzLib's IDs exactly; this is the term for IDs made of TDF frame and scan numbers.
                return new SpectraFormat("MS:1002817", "Bruker TDF format", "MS:1002818", "Bruker TDF nativeID format");
            }
        }
    }
}