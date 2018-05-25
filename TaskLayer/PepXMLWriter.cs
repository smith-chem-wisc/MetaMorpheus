﻿using EngineLayer;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Serialization;

namespace TaskLayer
{
    public static class PepXMLWriter
    {
        public static void WritePepXML(List<PeptideSpectralMatch> items, List<DbForTask> database, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, CommonParameters CommonParameters,string outputPath)
        {
            if (!items.Any())
                return;

            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(pepXML.Generated.msms_pipeline_analysis));
            var _pepxml = new pepXML.Generated.msms_pipeline_analysis();

            #region Add element to pepXML

            _pepxml.date = DateTime.Now;
            _pepxml.summary_xml = items[0].FullFilePath + ".pep.xml";

            string proteaseNC = string.Join("", CommonParameters.DigestionParams.Protease.SequencesPreventingCleavage);
            string proteaseC = string.Join("", CommonParameters.DigestionParams.Protease.SequencesInducingCleavage);

            string fileNameNoExtension = Path.GetFileNameWithoutExtension(items[0].FullFilePath);
            string filePathNoExtension = Path.ChangeExtension(items[0].FullFilePath, null);

            var para = new List<pepXML.Generated.nameValueType>();
            {
                para.Add(new pepXML.Generated.nameValueType { name = "threads", value = CommonParameters.MaxThreadsToUsePerFile.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "database", value = database.First().FilePath });
                para.Add(new pepXML.Generated.nameValueType { name = "MS_data_file", value = items[0].FullFilePath });

                para.Add(new pepXML.Generated.nameValueType { name = "MaxMissed Cleavages", value = CommonParameters.DigestionParams.MaxMissedCleavages.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Protease", value = CommonParameters.DigestionParams.Protease.Name });
                para.Add(new pepXML.Generated.nameValueType { name = "Initiator Methionine", value = CommonParameters.DigestionParams.InitiatorMethionineBehavior.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Modification Isoforms", value = CommonParameters.DigestionParams.MaxModificationIsoforms.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Min Peptide Len", value = CommonParameters.DigestionParams.MinPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Peptide Len", value = CommonParameters.DigestionParams.MaxPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Product Mass Tolerance", value = CommonParameters.ProductMassTolerance.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Ions to search", value = "B " + CommonParameters.BIons.ToString() + " Y " + CommonParameters.YIons.ToString() + " C " + CommonParameters.CIons.ToString() + " Z " + CommonParameters.ZdotIons.ToString() });
                foreach (var item in fixedModifications)
                {
                    para.Add(new pepXML.Generated.nameValueType { name = "Fixed Modifications: " + item.id, value = item.monoisotopicMass.ToString() });
                }
                foreach (var item in variableModifications)
                {
                    para.Add(new pepXML.Generated.nameValueType { name = "Variable Modifications: " + item.id, value = item.monoisotopicMass.ToString() });
                }

                para.Add(new pepXML.Generated.nameValueType { name = "Localize All Modifications", value = CommonParameters.LocalizeAll.ToString() });
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
                         //Contact for support
                         //search_engine = pepXML.Generated.engineType.MetaMorpheus 
                         search_engine_version = GlobalVariables.MetaMorpheusVersion,
                         precursor_mass_type = pepXML.Generated.massType.monoisotopic,
                         fragment_mass_type = pepXML.Generated.massType.monoisotopic,
                         search_id = 1,
                         //generate database information
                         search_database = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summarySearch_database
                         {
                             local_path = database.First().FilePath,
                             type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summarySearch_databaseType.AA,
                         },
                         enzymatic_search_constraint = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summaryEnzymatic_search_constraint
                         {
                             enzyme = CommonParameters.DigestionParams.Protease.Name,
                             max_num_internal_cleavages = CommonParameters.DigestionParams.MaxMissedCleavages.ToString(),
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
                var pepWithSetMods = items[i].CompactPeptides.First().Value.Item2.First();
                var mods = new List<pepXML.Generated.modInfoDataTypeMod_aminoacid_mass>();
                foreach (var mod in pepWithSetMods.allModsOneIsNterminus)
                {
                    var pepXmlMod = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass
                    {
                        mass = mod.Value.monoisotopicMass,
                        position = (mod.Key - 1).ToString()
                    };
                    mods.Add(pepXmlMod);
                }

                var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                {
                    hit_rank = 1,
                    peptide = ((items[i].BaseSequence != null) ? items[i].BaseSequence : "-"),
                    peptide_prev_aa = pepWithSetMods.PreviousAminoAcid.ToString(),
                    peptide_next_aa = pepWithSetMods.NextAminoAcid.ToString(),
                    protein = ((pepWithSetMods.Protein.Accession != null) ? pepWithSetMods.Protein.Accession : "-"),
                    num_tot_proteins = 1,
                    calc_neutral_pep_mass = (float)items[i].ScanPrecursorMass,
                    massdiff = ((items[i].PeptideMonisotopicMass != null) ? (items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value).ToString() : "-"),
                    modification_info = (mods.Count == 0 ? new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = mods.ToArray() } : null),
                    search_score = new pepXML.Generated.nameValueType[]
                    {
                                        new pepXML.Generated.nameValueType{ name = "Score", value = items[i].Score.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = items[i].FdrInfo.QValue.ToString() }
                    },
                };
                searchHits.Add(searchHit);
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

            #endregion Add element to pepXML

            TextWriter writer = new StreamWriter(Path.Combine(outputPath));
            _indexedSerializer.Serialize(writer, _pepxml);
            writer.Close();
        }
    }
}
