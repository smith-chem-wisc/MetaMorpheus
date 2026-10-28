using Chemistry;
using EngineLayer;
using EngineLayer.SpectrumMatch;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Serialization;
using EngineLayer.DatabaseLoading;
using Omics.Fragmentation.Peptide;
using Omics.Modifications;

namespace TaskLayer
{
    public static class PepXMLWriter
    {
        // the unmodified termini of a peptide: H on the N terminus, OH on the C terminus
        private static readonly double UnmodifiedNTerminusMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double UnmodifiedCTerminusMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass
            + PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        public static void WritePepXml(List<SpectralMatch> psms, List<DbForTask> database, List<Modification> variableModifications, List<Modification> fixedModifications, CommonParameters CommonParameters, string outputPath,
            FilteredPsms identificationFilter = null)
        {
            if (!psms.Any())
            {
                return;
            }

            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(pepXML.Generated.msms_pipeline_analysis));
            var _pepxml = new pepXML.Generated.msms_pipeline_analysis();

            _pepxml.date = DateTime.Now;
            _pepxml.summary_xml = psms[0].FullFilePath + ".pep.XML";

            string proteaseNC = string.Join(string.Empty, CommonParameters.DigestionParams.DigestionAgent.DigestionMotifs.Select(m => m.InducingCleavage));
            string proteaseC = string.Join(string.Empty, CommonParameters.DigestionParams.DigestionAgent.DigestionMotifs.Select(m => m.InducingCleavage));

            string fileNameNoExtension = Path.GetFileNameWithoutExtension(psms[0].FullFilePath);
            string filePathNoExtension = Path.ChangeExtension(psms[0].FullFilePath, null);

            var para = new List<pepXML.Generated.nameValueType>();
            {
                para.Add(new pepXML.Generated.nameValueType { name = "threads", value = CommonParameters.MaxThreadsToUsePerFile.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "database", value = database.First().FilePath });
                para.Add(new pepXML.Generated.nameValueType { name = "MS_data_file", value = psms[0].FullFilePath });

                para.Add(new pepXML.Generated.nameValueType { name = "MaxMissed Cleavages", value = CommonParameters.DigestionParams.MaxMissedCleavages.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Protease", value = CommonParameters.DigestionParams.DigestionAgent.Name });
                para.Add(new pepXML.Generated.nameValueType { name = "Initiator Methionine", value = ((DigestionParams)CommonParameters.DigestionParams).InitiatorMethionineBehavior.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Modification Isoforms", value = CommonParameters.DigestionParams.MaxModificationIsoforms.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Min Peptide Len", value = CommonParameters.DigestionParams.MinLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Peptide Len", value = CommonParameters.DigestionParams.MaxLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Product Mass Tolerance", value = CommonParameters.ProductMassTolerance.ToString() });
                // TODO: check this
                para.Add(new pepXML.Generated.nameValueType { name = "Ions to search", value = string.Join(", ", DissociationTypeCollection.ProductsFromDissociationType[CommonParameters.DissociationType]) });
                para.Add(new pepXML.Generated.nameValueType { name = "Q-value Filter", value = CommonParameters.QValueThreshold.ToString() });
                if (identificationFilter != null && identificationFilter.IdentificationFilter != IdentificationFilterMode.Thresholds)
                {
                    // The filter a named identification filter applied to these matches, e.g. "pep q-value < 0.01".
                    para.Add(new pepXML.Generated.nameValueType { name = "Identification Filter", value = identificationFilter.GetFilterTypeString() + " " + identificationFilter.GetThresholdString() });
                }
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
                     name = CommonParameters.DigestionParams.DigestionAgent.Name,
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
                             enzyme = CommonParameters.DigestionParams.DigestionAgent.Name,
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
                PeptideWithSetModifications peptide = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer as PeptideWithSetModifications;

                // pepXML gives the mass of the modified residue or terminus, not the modification's mass shift, and keeps
                // the termini apart from the residues
                var mods = new List<pepXML.Generated.modInfoDataTypeMod_aminoacid_mass>();
                double? nTerminusMass = null;
                double? cTerminusMass = null;
                foreach (var mod in peptide.AllModsOneIsNterminus.OrderBy(kv => kv.Key))
                {
                    int position = mod.Key - 1; // one-based residue; 0 is the N terminus, Length + 1 the C terminus
                    double shift = (double)mod.Value.MonoisotopicMass;
                    if (position == 0)
                    {
                        nTerminusMass = UnmodifiedNTerminusMass + shift;
                    }
                    else if (position == peptide.Length + 1)
                    {
                        cTerminusMass = UnmodifiedCTerminusMass + shift;
                    }
                    else
                    {
                        mods.Add(new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass
                        {
                            mass = (Residue.TryGetResidue(peptide.BaseSequence[position - 1], out Residue residue) ? residue.MonoisotopicMass : double.NaN) + shift,
                            position = position.ToString()
                        });
                    }
                }
                bool isModified = mods.Count > 0 || nTerminusMass.HasValue || cTerminusMass.HasValue;

                var proteinAccessions = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer.Parent.Accession).Distinct();

                var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                {
                    // TODO: handle PSM ambiguity if pepXML supports it (base sequence, mod localization, protein)
                    // TODO: add target/decoy/contaminant designation for each PSM
                    // TODO: add amino acid substitution
                    hit_rank = 1,
                    peptide = ((psm.BaseSequence != null) ? psm.BaseSequence : "Ambiguous"),
                    peptide_prev_aa = peptide.PreviousResidue.ToString(),
                    peptide_next_aa = peptide.NextResidue.ToString(),
                    protein = ((peptide.Protein.Accession != null) ? peptide.Protein.Accession : string.Join("|", proteinAccessions)),
                    num_tot_proteins = (uint)proteinAccessions.Count(),
                    calc_neutral_pep_mass = (float)((psm.BioPolymerWithSetModsMonoisotopicMass != null) ? psm.BioPolymerWithSetModsMonoisotopicMass : float.NaN),
                    massdiff = ((psm.BioPolymerWithSetModsMonoisotopicMass != null) ? (psm.ScanPrecursorMass - psm.BioPolymerWithSetModsMonoisotopicMass.Value).ToString() : "Ambiguous"),
                    modification_info = isModified
                        ? new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = mods.Count > 0 ? mods.ToArray() : null,
                            mod_nterm_mass = nTerminusMass ?? 0,
                            mod_nterm_massSpecified = nTerminusMass.HasValue,
                            mod_cterm_mass = cTerminusMass ?? 0,
                            mod_cterm_massSpecified = cTerminusMass.HasValue,
                        }
                        : null,
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