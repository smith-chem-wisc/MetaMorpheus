using EngineLayer;
using EngineLayer.CrosslinkSearch;
using Proteomics;
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
        #region Private Methods

        private void WriteCrosslinkToTsv(List<Tuple<PsmCross, PsmCross>> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("File Name\tScan Numer\tPrecusor MZ\tPrecusor charge\tPrecusor mass" +
                    "\tPep1\tPep1 Protein Access\tPep1 Base sequence\tPep1 Full sequence(crosslink site)\tPep1 mass\tPep1 XLBestScore\tPep1 topPos" +
                    "\tPep2\tPep2 Protein Access\tPep2 Base sequence\tPep2 Full sequence(crosslink site)\tPep2 mass\tPep2 XLBestScore\tPep2 topPos\tTotalScore\tMass diff\tQValue");
                foreach (var item in items)
                {
                    var x = item.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession);
                    output.WriteLine(
                        item.Item1.FullFilePath
                        + "\t" + item.Item1.ScanNumber.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.ScanPrecursorMonoisotopicPeak.ToString() //CultureInfo.InvariantCulture
                        + "\t" + item.Item1.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                        + "\t"
                        + "\t" + item.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.BaseSequence
                        + "\t" + item.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence + "(" + item.Item1.Xlpos.ToString(CultureInfo.InvariantCulture) + ")"
                        + "\t" + item.Item1.MostProbableProteinInfo.PeptideMonoisotopicMass.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.Score.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.TopPosition[0].ToString(CultureInfo.InvariantCulture)
                        //+ "\t" + item.Item1.NScore.ToString(CultureInfo.InvariantCulture)
                        + "\t"
                        + "\t" + item.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item2.BaseSequence
                        + "\t" + item.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence + "(" + item.Item2.Xlpos.ToString(CultureInfo.InvariantCulture) + ")"
                        + "\t" + item.Item2.MostProbableProteinInfo.PeptideMonoisotopicMass.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item2.Score.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.TopPosition[1].ToString(CultureInfo.InvariantCulture)
                        //+ "\t" + item.Item2.NScore.ToString(CultureInfo.InvariantCulture)

                        + "\t" + item.Item1.XLTotalScore.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (item.Item1.ScanPrecursorMass - item.Item2.MostProbableProteinInfo.PeptideMonoisotopicMass - item.Item1.MostProbableProteinInfo.PeptideMonoisotopicMass)
                        + "\t" + item.Item1.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture));
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        private void WritePepXML(List<Tuple<PsmCross, PsmCross>> items, List<DbForTask> dbFilenameList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> localizeableModifications, string outputFolder, string fileName, List<string> nestedIds)
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(pepXML.Generated.msms_pipeline_analysis));
            var _pepxml = new pepXML.Generated.msms_pipeline_analysis()
            {
                #region Add element to pepXML

                date = DateTime.Now,
                summary_xml = items[0].Item1.FullFilePath + ".pep.xml"
            };

            string proteaseC = ""; string proteaseNC = "";
            foreach (var x in Protease.SequencesInducingCleavage) { proteaseC += x; }
            foreach (var x in Protease.SequencesPreventingCleavage) { proteaseNC += x; }

            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass().SelectCrosslinker(CrosslinkerType);
            string modsFixed = ""; string modsVar = "";
            foreach (var x in ListOfModsFixed) { modsFixed += x.Item2 + "."; }
            foreach (var x in ListOfModsVariable) { modsVar += x.Item2 + "."; }

            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, SearchDecoy, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();
            uint proteinTot = Convert.ToUInt32(proteinList.Count);

            string fileNameNoExtension = Path.GetFileNameWithoutExtension(items[0].Item1.FullFilePath);
            string filePathNoExtension = Path.ChangeExtension(items[0].Item1.FullFilePath, null);

            _pepxml.msms_run_summary = new pepXML.Generated.msms_pipeline_analysisMsms_run_summary[1]
             {
                 new pepXML.Generated.msms_pipeline_analysisMsms_run_summary
                 {
                 base_name = filePathNoExtension,
                 raw_data_type = "raw",
                 raw_data = ".mzML",
                 sample_enzyme = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySample_enzyme()
                 {
                     name = Protease.Name,
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
                         //search_engine = pepXML.Generated.engineType.Kojak,
                         search_engine_version = GlobalEngineLevelSettings.MetaMorpheusVersion,
                         precursor_mass_type = pepXML.Generated.massType.monoisotopic,
                         fragment_mass_type = pepXML.Generated.massType.monoisotopic,
                         search_id = 1,
                         search_database = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summarySearch_database
                         {
                             local_path = dbFilenameList[0].FilePath,
                             type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summarySearch_databaseType.AA,
                         },
                         enzymatic_search_constraint = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySearch_summaryEnzymatic_search_constraint
                         {
                             enzyme = Protease.Name,
                             max_num_internal_cleavages = MaxMissedCleavages.ToString(),
                             //min_number_termini = "2"
                         },
                         parameter = new pepXML.Generated.nameValueType[]
                         {
                             new pepXML.Generated.nameValueType{ name = "threads", value = "" },
                             new pepXML.Generated.nameValueType{ name = "database", value = dbFilenameList[0].FilePath },
                             new pepXML.Generated.nameValueType{ name = "MS_data_file", value = items[0].Item1.FullFilePath },

                             new pepXML.Generated.nameValueType{ name = "Search with All Possible Beta Peptides", value = CrosslinkSearchWithAllBeta.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-link Precusor Mass Tolence", value = XLprecusorMsTl.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker Type", value = crosslinker.CrosslinkerName },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker mass", value = crosslinker.TotalMass.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker cleavable", value = crosslinker.Cleavable.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker cleavable long mass", value = crosslinker.CleaveMassLong.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker cleavable short mass", value = crosslinker.CleaveMassShort.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker xl site", value = crosslinker.CrosslinkerModSite.ToString() },

                             new pepXML.Generated.nameValueType{ name = "Generate decoy proteins", value = SearchDecoy.ToString() },
                             new pepXML.Generated.nameValueType{ name = "MaxMissed Cleavages", value = MaxMissedCleavages.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Protease", value = Protease.Name },
                             new pepXML.Generated.nameValueType{ name = "Initiator Methionine", value = InitiatorMethionineBehavior.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Max Modification Isoforms", value = MaxModificationIsoforms.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Min Peptide Len", value = MinPeptideLength.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Max Peptide Len", value = MaxPeptideLength.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Product Mass Tolerance", value = ProductMassTolerance.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Ions to search", value = "B "+ BIons.ToString() + " Y " + YIons.ToString() + " C " +CIons.ToString() + " Z " + ZdotIons.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Allowed Beta Precusor Mass Difference", value =  XLBetaPrecusorMsTl.ToString()},

                             new pepXML.Generated.nameValueType{ name = "Fixed Modifications", value = modsFixed },
                             new pepXML.Generated.nameValueType{ name = "Variable Modificaions", value = modsVar },
                             new pepXML.Generated.nameValueType{ name = "Localize All Modifications", value = LocalizeAll.ToString() },
                         }
                     }
                 },
                 }
             };

            _pepxml.msms_run_summary[0].spectrum_query = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_query[items.Count];
            for (int i = 0; i < items.Count; i++)
            {
                _pepxml.msms_run_summary[0].spectrum_query[i] = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_query()
                {
                    spectrum = fileNameNoExtension + "." + items[i].Item1.ScanNumber.ToString(),
                    start_scan = Convert.ToUInt32(items[i].Item1.ScanNumber),
                    end_scan = Convert.ToUInt32(items[i].Item1.ScanNumber),
                    precursor_neutral_mass = (float)items[i].Item1.ScanPrecursorMonoisotopicPeak.Mz * items[i].Item1.ScanPrecursorCharge,
                    assumed_charge = items[i].Item1.ScanPrecursorCharge.ToString(),
                    index = Convert.ToUInt32(i + 1),
                    retention_time_sec = (float)items[i].Item1.ScanRetentionTime,
                    search_result = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result[1]
                    {
                        new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result
                        {
                            search_hit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit[1]
                            {
                                new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                                {
                                    hit_rank = 1,
                                    peptide = "-", peptide_prev_aa="-", peptide_next_aa="-", protein="-", num_tot_proteins = 1,
                                    calc_neutral_pep_mass = (float)items[i].Item1.ScanPrecursorMonoisotopicPeak.Mz * items[i].Item1.ScanPrecursorCharge,
                                    massdiff = (items[i].Item1.ScanPrecursorMass - items[i].Item2.MostProbableProteinInfo.PeptideMonoisotopicMass - items[i].Item1.MostProbableProteinInfo.PeptideMonoisotopicMass - crosslinker.TotalMass).ToString(),
                                    xlink_typeSpecified = true,
                                    xlink_type = pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink_type.xl,
                                    xlink = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlink
                                    {
                                        identifier = crosslinker.CrosslinkerName,
                                        mass = (float)crosslinker.TotalMass,
                                        linked_peptide = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide[2]
                                        {
                                            new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                                            {
                                                peptide = items[i].Item1.BaseSequence,
                                                peptide_prev_aa = items[i].Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().PreviousAminoAcid.ToString(),
                                                peptide_next_aa = items[i].Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().NextAminoAcid.ToString(),
                                                protein = items[i].Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().Protein.Accession,
                                                num_tot_proteins = proteinTot/2,
                                                calc_neutral_pep_mass = (float)items[i].Item1.MostProbableProteinInfo.PeptideMonoisotopicMass,
                                                complement_mass = (float)(items[i].Item1.ScanPrecursorMass - items[i].Item1.MostProbableProteinInfo.PeptideMonoisotopicMass),
                                                designation = "alpha",

                                                xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].Item1.XLBestScore.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].Item1.Xlpos.ToString() },
                                                }
                                            },
                                            new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                                            {
                                                peptide = items[i].Item2.BaseSequence,
                                                peptide_prev_aa = items[i].Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().PreviousAminoAcid.ToString(),
                                                peptide_next_aa = items[i].Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().NextAminoAcid.ToString(),
                                                protein = items[i].Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().Protein.Accession,
                                                num_tot_proteins = proteinTot,
                                                calc_neutral_pep_mass = (float)items[i].Item2.MostProbableProteinInfo.PeptideMonoisotopicMass,
                                                complement_mass = (float)(items[i].Item1.ScanPrecursorMass - items[i].Item2.MostProbableProteinInfo.PeptideMonoisotopicMass),
                                                designation = "beta",

                                                xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].Item2.XLBestScore.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].Item2.Xlpos.ToString() },
                                                }
                                            }
                                        }
                                    },
                                    search_score = new pepXML.Generated.nameValueType[]
                                    {
                                        new pepXML.Generated.nameValueType{ name = "xlTotalScore", value = items[i].Item1.XLTotalScore.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = items[i].Item1.FdrInfo.QValue.ToString() }
                                    }
                                }
                            }
                        }
                    }
                };

                #region mods infomation

                int modsFixedNum1 = items[i].Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
                int modsFixedNum2 = items[i].Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
                if (modsFixedNum1 != 0)
                {
                    modsFixedNum1 = items[i].Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
                    if (modsFixedNum1 == 1)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[1]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum1 == 2)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[2]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum1 == 3)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[3]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum1 == 4)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[4]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    for (int j = 0; j < modsFixedNum1; j++)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info.mod_aminoacid_mass[j].mass = items[i].Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Values.Select(p => p.monoisotopicMass).ToList()[j];
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info.mod_aminoacid_mass[j].position = items[i].Item1.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Keys.ToList()[j].ToString();
                    }
                }
                if (modsFixedNum2 != 0)
                {
                    modsFixedNum2 = items[i].Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
                    if (modsFixedNum2 == 1)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[1]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum2 == 2)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[2]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum2 == 3)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[3]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum2 == 4)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[4]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    for (int j = 0; j < modsFixedNum2; j++)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info.mod_aminoacid_mass[j].mass = items[i].Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Values.Select(p => p.monoisotopicMass).ToList()[j];
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info.mod_aminoacid_mass[j].position = items[i].Item2.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Keys.ToList()[j].ToString();
                    }
                }

                #endregion mods infomation
            }

            #endregion Add element to pepXML

            TextWriter writer = new StreamWriter(Path.Combine(outputFolder, fileName + ".pep.xml"));
            _indexedSerializer.Serialize(writer, _pepxml);
            writer.Close();
            SucessfullyFinishedWritingFile(Path.Combine(outputFolder, fileName + ".pep.xml"), nestedIds);
        }

        #endregion Private Methods
    }
}