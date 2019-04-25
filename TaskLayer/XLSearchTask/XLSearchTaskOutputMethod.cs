using EngineLayer;
using EngineLayer.CrosslinkSearch;
using Proteomics;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml.Serialization;

namespace TaskLayer
{
    public partial class XLSearchTask : MetaMorpheusTask
    {
        public static void WritePsmCrossToTsv(List<CrosslinkSpectralMatch> items, string filePath, int writeType)
        {
            if (items.Count == 0)
            {
                return;
            }

            using (StreamWriter output = new StreamWriter(filePath))
            {
                string header = "";
                switch (writeType)
                {
                    case 1:
                        header = CrosslinkSpectralMatch.GetTabSepHeaderSingle();
                        break;
                    case 2:
                        header = CrosslinkSpectralMatch.GetTabSepHeaderCross();
                        break;
                    case 3:
                        header = CrosslinkSpectralMatch.GetTabSepHeaderGlyco();
                        break;
                    default:
                        break;
                }
                output.WriteLine(header);
                foreach (var heh in items)
                {
                    output.WriteLine(heh.ToString());
                }
            }
        }
        
        public void WriteCrosslinkToTxtForPercolator(List<CrosslinkSpectralMatch> items, string outputFolder, string fileName, Crosslinker crosslinker, List<string> nestedIds)
        {
            if (items.Count == 0)
            { return; }
            var writtenFile = Path.Combine(outputFolder, fileName + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum" +
                    "\tPeptide\tProtein");
                foreach (var item in items)
                {
                    if (item.BaseSequence != null && item.BetaPeptide.BaseSequence != null && item.ProteinAccession != null && item.BetaPeptide.ProteinAccession != null)
                    {
                        string x = "T"; int label = 1;
                        if (item.IsDecoy || item.BetaPeptide.IsDecoy)
                        {
                            x = "D"; label = -1;
                        }
                        output.WriteLine(
                            x + "-" + item.ScanNumber.ToString(CultureInfo.InvariantCulture) + "-" + item.ScanRetentionTime.ToString(CultureInfo.InvariantCulture)
                            + "\t" + label.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.XLTotalScore.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.DeltaScore.ToString(CultureInfo.InvariantCulture)
                            + "\t" + (item.XlRank[0] + item.XlRank[1]).ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                            + "\t" + ((item.PeptideMonisotopicMass.HasValue && item.BetaPeptide.PeptideMonisotopicMass.HasValue) ? ((item.ScanPrecursorMass - item.BetaPeptide.PeptideMonisotopicMass.Value - item.PeptideMonisotopicMass.Value - crosslinker.TotalMass) / item.ScanPrecursorMass * 1E6).ToString(CultureInfo.InvariantCulture) : "---")
                            + "\t" + item.BetaPeptide.BaseSequence.Length.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.BaseSequence.Length.ToString(CultureInfo.InvariantCulture)
                            + "\t" + (item.BetaPeptide.BaseSequence.Length + item.BaseSequence.Length).ToString(CultureInfo.InvariantCulture)
                            + "\t" + "-." + item.FullSequence + item.LinkPositions.First().ToString(CultureInfo.InvariantCulture) + "--" + item.BetaPeptide.FullSequence + item.BetaPeptide.LinkPositions.First().ToString(CultureInfo.InvariantCulture) + ".-"
                            + "\t" + item.BestMatchingPeptides.First().Peptide.Protein.Accession.ToString(CultureInfo.InvariantCulture)
                                   + "(" + item.XlProteinPos.ToString(CultureInfo.InvariantCulture) + ")"
                            + "\t" + item.BetaPeptide.BestMatchingPeptides.First().Peptide.Protein.Accession.ToString(CultureInfo.InvariantCulture)
                                   + "(" + item.BetaPeptide.XlProteinPos.ToString(CultureInfo.InvariantCulture) + ")"
                            );
                    }
                }
            }
            FinishedWritingFile(writtenFile, nestedIds);
        }
        
        public void WritePepXML_xl(List<CrosslinkSpectralMatch> items, List<Protein> proteinList, string databasePath, List<Modification> variableModifications, List<Modification> fixedModifications, List<string> localizeableModificationTypes, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (!items.Any())
            {
                return;
            }

            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(pepXML.Generated.msms_pipeline_analysis));
            var _pepxml = new pepXML.Generated.msms_pipeline_analysis();

            _pepxml.date = DateTime.Now;
            _pepxml.summary_xml = items[0].FullFilePath + ".pep.XML";

            string proteaseC = ""; string proteaseNC = "";
            foreach (var x in CommonParameters.DigestionParams.Protease.DigestionMotifs.Select(m => m.InducingCleavage)) { proteaseC += x; }
            foreach (var x in CommonParameters.DigestionParams.Protease.DigestionMotifs.Select(m => m.PreventingCleavage)) { proteaseNC += x; }

            Crosslinker crosslinker = XlSearchParameters.Crosslinker;

            string fileNameNoExtension = Path.GetFileNameWithoutExtension(items[0].FullFilePath);
            string filePathNoExtension = Path.ChangeExtension(items[0].FullFilePath, null);
            string modSites = crosslinker.CrosslinkerModSites.ToCharArray().Concat(crosslinker.CrosslinkerModSites2.ToCharArray()).Distinct().ToString();

            var para = new List<pepXML.Generated.nameValueType>();
            {
                para.Add(new pepXML.Generated.nameValueType { name = "threads", value = CommonParameters.MaxThreadsToUsePerFile.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "database", value = databasePath });
                para.Add(new pepXML.Generated.nameValueType { name = "MS_data_file", value = items[0].FullFilePath });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-link precursor Mass Tolerance", value = CommonParameters.PrecursorMassTolerance.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker type", value = crosslinker.CrosslinkerName });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker mass", value = crosslinker.TotalMass.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker cleavable", value = crosslinker.Cleavable.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker cleavable long mass", value = crosslinker.CleaveMassLong.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker cleavable short mass", value = crosslinker.CleaveMassShort.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker xl site", value = modSites });

                para.Add(new pepXML.Generated.nameValueType { name = "Generate decoy proteins", value = XlSearchParameters.DecoyType.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "MaxMissed Cleavages", value = CommonParameters.DigestionParams.MaxMissedCleavages.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Protease", value = CommonParameters.DigestionParams.Protease.Name });
                para.Add(new pepXML.Generated.nameValueType { name = "Initiator Methionine", value = CommonParameters.DigestionParams.InitiatorMethionineBehavior.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Modification Isoforms", value = CommonParameters.DigestionParams.MaxModificationIsoforms.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Min Peptide Len", value = CommonParameters.DigestionParams.MinPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Peptide Len", value = CommonParameters.DigestionParams.MaxPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Product Mass Tolerance", value = CommonParameters.ProductMassTolerance.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Ions to search", value = String.Join(", ", DissociationTypeCollection.ProductsFromDissociationType[CommonParameters.DissociationType]) });

                foreach (var fixedMod in fixedModifications)
                {
                    para.Add(new pepXML.Generated.nameValueType { name = "Fixed Modifications: " + fixedMod.IdWithMotif, value = fixedMod.MonoisotopicMass.ToString() });
                }
                foreach (var variableMod in variableModifications)
                {
                    para.Add(new pepXML.Generated.nameValueType { name = "Variable Modifications: " + variableMod.IdWithMotif, value = variableMod.MonoisotopicMass.ToString() });
                }

                para.Add(new pepXML.Generated.nameValueType { name = "Localize All Modifications", value = "true" });
            }

            _pepxml.msms_run_summary = new pepXML.Generated.msms_pipeline_analysisMsms_run_summary[1]
             {
                 new pepXML.Generated.msms_pipeline_analysisMsms_run_summary
                 {
                 base_name = filePathNoExtension,
                 raw_data_type = "raw",
                 raw_data = ".mzML",
                 sample_enzyme = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySample_enzyme()
                 {
                     name = CommonParameters.DigestionParams.Protease.Name,
                     specificity = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificity[1]
                     {
                         new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificity
                         {
                             cut = proteaseC,
                             no_cut = proteaseNC,
                         }
                     }
                 },

                 search_summary = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summary[1]
                 {
                     new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summary
                     {
                         base_name = filePathNoExtension,
                         search_engine_version = GlobalVariables.MetaMorpheusVersion,
                         precursor_mass_type = pepXML.Generated.massType.monoisotopic,
                         fragment_mass_type = pepXML.Generated.massType.monoisotopic,
                         search_id = 1,
                         search_database = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summarySearch_database
                         {
                             local_path = databasePath,
                             type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summarySearch_databaseType.AA,
                         },
                         enzymatic_search_constraint = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summaryEnzymatic_search_constraint
                         {
                             enzyme = CommonParameters.DigestionParams.Protease.Name,
                             max_num_internal_cleavages = CommonParameters.DigestionParams.MaxMissedCleavages.ToString(),
                             //min_number_termini = "2"
                         },

                         parameter = para.ToArray()
                     }
                 },
                 }
             };

            _pepxml.msms_run_summary[0].spectrum_query = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_query[items.Count];

            var searchHits = new List<pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit>();
            for (int i = 0; i < items.Count; i++)
            {
                var mods = new List<pepXML.Generated.modInfoDataTypeMod_aminoacid_mass>();
                var alphaPeptide = items[i].BestMatchingPeptides.First().Peptide;

                foreach (var modification in alphaPeptide.AllModsOneIsNterminus)
                {
                    var mod = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass
                    {
                        mass = modification.Value.MonoisotopicMass.Value,

                        // convert from one-based to zero-based (N-term is zero in the pepXML output)
                        position = (modification.Key - 1).ToString()
                    };
                    mods.Add(mod);
                }

                if (items[i].CrossType == PsmCrossType.Single)
                {
                    var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                    {

                        hit_rank = 1,
                        peptide = alphaPeptide.BaseSequence,
                        peptide_prev_aa = alphaPeptide.PreviousAminoAcid.ToString(),
                        peptide_next_aa = alphaPeptide.NextAminoAcid.ToString(),
                        protein = alphaPeptide.Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMass,
                        massdiff = (items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value).ToString(),
                        xlink_typeSpecified = true,
                        xlink_type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type.na,
                        modification_info = new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = mods.ToArray() },
                        search_score = new pepXML.Generated.nameValueType[]
                                    {
                                        new pepXML.Generated.nameValueType{ name = "xlTotalScore", value = items[i].XLTotalScore.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = items[i].FdrInfo.QValue.ToString() }
                                    },
                    };
                    searchHits.Add(searchHit);
                }
                else if (items[i].CrossType == PsmCrossType.DeadEnd || items[i].CrossType == PsmCrossType.DeadEndH2O || items[i].CrossType == PsmCrossType.DeadEndNH2 || items[i].CrossType == PsmCrossType.DeadEndTris)
                {
                    double crosslinkerDeadEndMass = 0;
                    switch (items[i].CrossType)
                    {
                        case PsmCrossType.DeadEndNH2:
                            crosslinkerDeadEndMass = crosslinker.DeadendMassNH2;
                            break;

                        case PsmCrossType.DeadEndTris:
                            crosslinkerDeadEndMass = crosslinker.DeadendMassTris;
                            break;

                        default:
                            crosslinkerDeadEndMass = crosslinker.DeadendMassH2O;
                            break;
                    }
                    var mod = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass { mass = crosslinkerDeadEndMass, position = items[i].LinkPositions.First().ToString() };
                    mods.Add(mod);
                    var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                    {
                        hit_rank = 1,
                        peptide = alphaPeptide.BaseSequence,
                        peptide_prev_aa = alphaPeptide.PreviousAminoAcid.ToString(),
                        peptide_next_aa = alphaPeptide.NextAminoAcid.ToString(),
                        protein = alphaPeptide.Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMass,
                        massdiff = (items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value - crosslinkerDeadEndMass).ToString(),
                        xlink_typeSpecified = true,
                        xlink_type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type.na,
                        modification_info = new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = mods.ToArray() },
                        search_score = new pepXML.Generated.nameValueType[]
                                    {
                                        new pepXML.Generated.nameValueType{ name = "xlTotalScore", value = items[i].XLTotalScore.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = items[i].FdrInfo.QValue.ToString() }
                                    },
                    };
                    searchHits.Add(searchHit);
                }
                else if (items[i].CrossType == PsmCrossType.Inter || items[i].CrossType == PsmCrossType.Intra || items[i].CrossType == PsmCrossType.Cross)
                {
                    var betaPeptide = items[i].BetaPeptide.BestMatchingPeptides.First().Peptide;
                    var modsBeta = new List<pepXML.Generated.modInfoDataTypeMod_aminoacid_mass>();

                    foreach (var mod in betaPeptide.AllModsOneIsNterminus)
                    {
                        var modBeta = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass
                        {
                            mass = mod.Value.MonoisotopicMass.Value,

                            // convert from one-based to zero-based (N-term is zero in the pepXML output)
                            position = (mod.Key - 1).ToString()
                        };
                        modsBeta.Add(modBeta);
                    }

                    var alpha = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                    {
                        peptide = alphaPeptide.BaseSequence,
                        peptide_prev_aa = alphaPeptide.PreviousAminoAcid.ToString(),
                        peptide_next_aa = alphaPeptide.NextAminoAcid.ToString(),
                        protein = alphaPeptide.Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].PeptideMonisotopicMass.Value,
                        complement_mass = (float)(items[i].ScanPrecursorMass - alphaPeptide.MonoisotopicMass),
                        designation = "alpha",
                        modification_info = new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = mods.ToArray() },
                        xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].XLTotalScore.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].LinkPositions.First().ToString() },
                                                }
                    };
                    var beta = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                    {
                        peptide = betaPeptide.BaseSequence,
                        peptide_prev_aa = betaPeptide.PreviousAminoAcid.ToString(),
                        peptide_next_aa = betaPeptide.NextAminoAcid.ToString(),
                        protein = betaPeptide.Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)betaPeptide.MonoisotopicMass,
                        complement_mass = (float)(items[i].ScanPrecursorMass - betaPeptide.MonoisotopicMass),
                        designation = "beta",
                        modification_info = new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = modsBeta.ToArray() },
                        xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].BetaPeptide.Score.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].BetaPeptide.LinkPositions.First().ToString() },
                                                }
                    };
                    var cross = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide[2] { alpha, beta };
                    var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                    {
                        hit_rank = 1,
                        peptide = "-",
                        peptide_prev_aa = "-",
                        peptide_next_aa = "-",
                        protein = "-",
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMass,
                        massdiff = (items[i].ScanPrecursorMass - betaPeptide.MonoisotopicMass - alphaPeptide.MonoisotopicMass - crosslinker.TotalMass).ToString(),
                        xlink_typeSpecified = true,
                        xlink_type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type.xl,
                        xlink = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink
                        {
                            identifier = crosslinker.CrosslinkerName,
                            mass = (float)crosslinker.TotalMass,
                            linked_peptide = cross
                        },
                        search_score = new pepXML.Generated.nameValueType[]
                                    {
                                        new pepXML.Generated.nameValueType{ name = "xlTotalScore", value = items[i].XLTotalScore.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = items[i].FdrInfo.QValue.ToString() }
                                    }
                    };
                    searchHits.Add(searchHit);
                }
                else if (items[i].CrossType == PsmCrossType.Loop)
                {
                    var thePeptide = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                    {
                        xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "link", value = items[i].LinkPositions.First().ToString() },
                                                    new pepXML.Generated.nameValueType{ name = "link", value = items[i].LinkPositions[1].ToString() }
                                                }
                    };
                    var cross = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide[1] { thePeptide };
                    var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                    {
                        hit_rank = 1,
                        peptide = alphaPeptide.BaseSequence,
                        peptide_prev_aa = alphaPeptide.PreviousAminoAcid.ToString(),
                        peptide_next_aa = alphaPeptide.NextAminoAcid.ToString(),
                        protein = alphaPeptide.Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMass,
                        massdiff = (items[i].ScanPrecursorMass - alphaPeptide.MonoisotopicMass - crosslinker.LoopMass).ToString(),
                        xlink_typeSpecified = true,
                        xlink_type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type.loop,
                        modification_info = new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = mods.ToArray() },
                        xlink = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink
                        {
                            identifier = crosslinker.CrosslinkerName,
                            mass = (float)crosslinker.TotalMass,
                            linked_peptide = cross
                        },
                        search_score = new pepXML.Generated.nameValueType[]
                                    {
                                        new pepXML.Generated.nameValueType{ name = "xlTotalScore", value = items[i].XLTotalScore.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = items[i].FdrInfo.QValue.ToString() }
                                    }
                    };
                    searchHits.Add(searchHit);
                }
            }

            for (int i = 0; i < items.Count; i++)
            {
                _pepxml.msms_run_summary[0].spectrum_query[i] = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_query()
                {
                    spectrum = fileNameNoExtension + "." + items[i].ScanNumber.ToString(),
                    start_scan = Convert.ToUInt32(items[i].ScanNumber),
                    end_scan = Convert.ToUInt32(items[i].ScanNumber),
                    precursor_neutral_mass = (float)items[i].ScanPrecursorMass,
                    assumed_charge = items[i].ScanPrecursorCharge.ToString(),
                    index = Convert.ToUInt32(i + 1),
                    retention_time_sec = (float)(items[i].ScanRetentionTime * 60),
                    search_result = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result[1]
                    {
                        new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result
                        {
                            search_hit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit[1]
                            {
                                searchHits[i]
                            }
                        }
                    }
                };
            }

            TextWriter writer = new StreamWriter(Path.Combine(outputFolder, fileName + ".pep.XML"));
            _indexedSerializer.Serialize(writer, _pepxml);
            writer.Close();
            FinishedWritingFile(Path.Combine(outputFolder, fileName + ".pep.XML"), nestedIds);
        }

        public void WriteOxoniumIons(Tuple<int, double[]>[] items, string filePath)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine("ScanNum\t109\t115\t126\t127\t138\t144\t163\t168\t186\t204\t274\t290\t292\t308\t366\t657\t673");
                foreach (var item in items)
                {
                    output.Write(item.Item1.ToString());
                    foreach (var d in item.Item2)
                    {
                        output.Write('\t' + (d/item.Item2.Max()).ToString());
                    }
                    output.Write("\r\n");
                }
            }
        }
    }
}