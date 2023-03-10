using EngineLayer;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Serialization;

namespace TaskLayer
{
    public static class PepXMLWriter
    {
        public static void WritePepXml(List<PeptideSpectralMatch> psms, List<DbForTask> database, List<Modification> variableModifications, List<Modification> fixedModifications, CommonParameters CommonParameters, string outputPath, double qValueFilter)
        {
            // TODO: needs a unit test
            psms = psms.Where(p => p.FdrInfo.QValue <= qValueFilter && p.FdrInfo.QValueNotch < qValueFilter).ToList();

            if (!psms.Any())
            {
                return;
            }

            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(pepXML.Generated.msms_pipeline_analysis));
            var _pepxml = new pepXML.Generated.msms_pipeline_analysis();

            _pepxml.date = DateTime.Now;
            _pepxml.summary_xml = psms[0].FullFilePath + ".pep.XML";

            string proteaseNC = string.Join(string.Empty, CommonParameters.DigestionParams.Protease.DigestionMotifs.Select(m => m.InducingCleavage));
            string proteaseC = string.Join(string.Empty, CommonParameters.DigestionParams.Protease.DigestionMotifs.Select(m => m.InducingCleavage));

            string fileNameNoExtension = Path.GetFileNameWithoutExtension(psms[0].FullFilePath);
            string filePathNoExtension = Path.ChangeExtension(psms[0].FullFilePath, null);

            var para = new List<pepXML.Generated.nameValueType>();
            {
                para.Add(new pepXML.Generated.nameValueType { name = "threads", value = CommonParameters.MaxThreadsToUsePerFile.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "database", value = database.First().FilePath });
                para.Add(new pepXML.Generated.nameValueType { name = "MS_data_file", value = psms[0].FullFilePath });

                para.Add(new pepXML.Generated.nameValueType { name = "MaxMissed Cleavages", value = CommonParameters.DigestionParams.MaxMissedCleavages.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Protease", value = CommonParameters.DigestionParams.Protease.Name });
                para.Add(new pepXML.Generated.nameValueType { name = "Initiator Methionine", value = CommonParameters.DigestionParams.InitiatorMethionineBehavior.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Modification Isoforms", value = CommonParameters.DigestionParams.MaxModificationIsoforms.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Min Peptide Len", value = CommonParameters.DigestionParams.MinPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Peptide Len", value = CommonParameters.DigestionParams.MaxPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Product Mass Tolerance", value = CommonParameters.ProductMassTolerance.ToString() });
                // TODO: check this
                para.Add(new pepXML.Generated.nameValueType { name = "Ions to search", value = string.Join(", ", DissociationTypeCollection.ProductsFromDissociationType[CommonParameters.DissociationType]) });
                para.Add(new pepXML.Generated.nameValueType { name = "Q-value Filter", value = CommonParameters.QValueOutputFilter.ToString() });
                foreach (var item in fixedModifications)
                {
                    para.Add(new pepXML.Generated.nameValueType { name = "Fixed Modifications: " + item.IdWithMotif, value = item.MonoisotopicMass.ToString() });
                }
                foreach (var item in variableModifications)
                {
                    para.Add(new pepXML.Generated.nameValueType { name = "Variable Modifications: " + item.IdWithMotif, value = item.MonoisotopicMass.ToString() });
                }

                para.Add(new pepXML.Generated.nameValueType { name = "Localize All Modifications", value = "true" });
            }

            _pepxml.msms_run_summary = new pepXML.Generated.msms_pipeline_analysisMsms_run_summary[1]
             {
                 new pepXML.Generated.msms_pipeline_analysisMsms_run_summary
                 {
                 base_name = filePathNoExtension,
                 raw_data_type = "raw",

                 raw_data = ".mzML", //TODO: use file format of spectra file used
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

                         // TODO: get MetaMorpheus recognized as a search engine type
                         //search_engine = pepXML.Generated.engineType.MetaMorpheus
                         search_engine_version = GlobalVariables.MetaMorpheusVersion,
                         precursor_mass_type = pepXML.Generated.massType.monoisotopic,
                         fragment_mass_type = pepXML.Generated.massType.monoisotopic,
                         search_id = 1,
                         //generate database information
                         //TODO: multiple databases
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

            _pepxml.msms_run_summary[0].spectrum_query = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_query[psms.Count];

            var searchHits = new List<pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit>();

            foreach (var psm in psms)
            {
                PeptideWithSetModifications peptide = psm.BestMatchingPeptides.First().Peptide;

                var mods = new List<pepXML.Generated.modInfoDataTypeMod_aminoacid_mass>();
                foreach (var mod in peptide.AllModsOneIsNterminus)
                {
                    var pepXmlMod = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass
                    {
                        mass = (double)mod.Value.MonoisotopicMass,
                        position = (mod.Key - 1).ToString()
                    };
                    mods.Add(pepXmlMod);
                }

                var proteinAccessions = psm.BestMatchingPeptides.Select(p => p.Peptide.Protein.Accession).Distinct();

                var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                {
                    // TODO: handle PSM ambiguity if pepXML supports it (base sequence, mod localization, protein)
                    // TODO: add target/decoy/contaminant designation for each PSM
                    // TODO: add amino acid substitution
                    hit_rank = 1,
                    peptide = ((psm.BaseSequence != null) ? psm.BaseSequence : "Ambiguous"),
                    peptide_prev_aa = peptide.PreviousAminoAcid.ToString(),
                    peptide_next_aa = peptide.NextAminoAcid.ToString(),
                    protein = ((peptide.Protein.Accession != null) ? peptide.Protein.Accession : string.Join("|", proteinAccessions)),
                    num_tot_proteins = (uint)proteinAccessions.Count(),
                    calc_neutral_pep_mass = (float)((psm.PeptideMonisotopicMass != null) ? psm.PeptideMonisotopicMass : float.NaN),
                    massdiff = ((psm.PeptideMonisotopicMass != null) ? (psm.ScanPrecursorMass - psm.PeptideMonisotopicMass.Value).ToString() : "Ambiguous"),
                    modification_info = (mods.Count == 0 ? new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = mods.ToArray() } : null),
                    search_score = new pepXML.Generated.nameValueType[]
                    {
                                        new pepXML.Generated.nameValueType{ name = "Score", value = psm.Score.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = psm.FdrInfo.QValue.ToString() }
                    },
                };
                searchHits.Add(searchHit);
            }

            for (int i = 0; i < psms.Count; i++)
            {
                _pepxml.msms_run_summary[0].spectrum_query[i] = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_query()
                {
                    spectrum = fileNameNoExtension + "." + psms[i].ScanNumber.ToString(),
                    start_scan = Convert.ToUInt32(psms[i].ScanNumber),
                    end_scan = Convert.ToUInt32(psms[i].ScanNumber),
                    precursor_neutral_mass = (float)psms[i].ScanPrecursorMass,
                    assumed_charge = psms[i].ScanPrecursorCharge.ToString(),
                    index = Convert.ToUInt32(i + 1),
                    retention_time_sec = (float)(psms[i].ScanRetentionTime * 60),
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

            TextWriter writer = new StreamWriter(Path.Combine(outputPath));
            _indexedSerializer.Serialize(writer, _pepxml);
            writer.Close();
        }
    }
}