﻿using EngineLayer;
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

        public void WriteCrosslinkToTsv(List<PsmCross> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("File Name\tScan Numer\tPrecusor MZ\tPrecusor charge\tPrecusor mass\tCrossType" +
                    "\tPep1\tPep1 Protein Access\tProtein link site\tPep1 Base sequence(crosslink site)\tPep1 Full sequence\tPep1 mass\tPep1 XLBestScore\tPep1 Rank" +
                    "\tPep2\tPep2 Protein Access\tProtein link site\tPep2 Base sequence(crosslink site)\tPep2 Full sequence\tPep2 mass\tPep2 XLBestScore\tPep2 Rank" +
                    "\tSummary\tQvalueTotalScore\tMass diff\tQValue\tParentIons\tParentIonsNum\tParentIonMaxIntensityRank\tCharge2Number\tLabel");
                foreach (var item in items)
                {
                    string label = "1";
                    if (item.IsDecoy || item.BetaPsmCross.IsDecoy)
                    {
                        label = "-1";
                    }

                    output.WriteLine(
                                            item.FullFilePath
                                            + "\t" + item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.ScanPrecursorMonoisotopicPeakMz.ToString() //CultureInfo.InvariantCulture
                                            + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.CrossType.ToString()
                                            + "\t"
                                            + "\t" + item.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.XlProteinPos.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.BaseSequence + "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                            + "\t" + item.FullSequence
                                            + "\t" + (item.PeptideMonisotopicMass.HasValue ? item.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                                            //"\t" + item.Score.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.XLBestScore.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.XlRank[0].ToString(CultureInfo.InvariantCulture)
                                            + "\t"
                                            + "\t" + item.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.BetaPsmCross.XlProteinPos.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.BetaPsmCross.BaseSequence + "(" + item.BetaPsmCross.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                            + "\t" + item.BetaPsmCross.FullSequence
                                            + "\t" + (item.BetaPsmCross.PeptideMonisotopicMass.HasValue ? item.BetaPsmCross.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                                            //"\t" + item.BetaPsmCross.Score.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.BetaPsmCross.XLBestScore.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.XlRank[1].ToString(CultureInfo.InvariantCulture)
                                            + "\t"
                                            + "\t" + item.XLQvalueTotalScore.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + ((item.PeptideMonisotopicMass.HasValue && item.BetaPsmCross.PeptideMonisotopicMass.HasValue) ? (item.BetaPsmCross.ScanPrecursorMass - item.BetaPsmCross.PeptideMonisotopicMass.Value - item.PeptideMonisotopicMass.Value).ToString(CultureInfo.InvariantCulture) : "---")
                                            + "\t" + (item.FdrInfo != null ? item.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture) : "-")
                                            + "\t" + item.ParentIonExist + "." + item.BetaPsmCross.ParentIonExist
                                            + "\t" + item.ParentIonExistNum.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + ((item.ParentIonMaxIntensityRanks != null) && (item.ParentIonMaxIntensityRanks.Any()) ? item.ParentIonMaxIntensityRanks.Min().ToString(CultureInfo.InvariantCulture) : "-")
                                            + "\t" + (item.Charge2IonExist + item.BetaPsmCross.Charge2IonExist).ToString(CultureInfo.InvariantCulture)
                                            + "\t" + label
                                            );
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        public void WriteAllToTsv(List<PsmCross> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("File Name\tScan Numer\tPrecusor MZ\tPrecusor charge\tPrecusor mass\tCross-link type" +
                    "\tPep1\tPep1 Protein Access(Protein link site)\tPep1 Base sequence(crosslink site)\tPep1 Full sequence\tPep1 mass\tPep1 Score\tPep1 XLBestScore\tPep1 Rank" +
                    "\tPep2\tPep2 Protein Access(Protein link site)\tPep2 Base sequence(crosslink site)\tPep2 Full sequence\tPep2 mass\tPep2 Score\tPep2 XLBestScore\tPep2 Rank" +
                    "\tSummary\tTotalScore\tMass diff\tQValue\tParentIons\tCharge2Number");
                foreach (var item in items)
                {
                    if (item.CrossType != PsmCrossType.Cross || item.CrossType != PsmCrossType.Inter || item.CrossType != PsmCrossType.Intra)
                    {
                        string position = "";
                        switch (item.CrossType)
                        {
                            case PsmCrossType.Singe:
                                break;
                            case PsmCrossType.Loop:
                                position = "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + "-" + item.XlPos2.ToString(CultureInfo.InvariantCulture) + ")";
                                break;
                            default:
                                position = "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + ")";
                                break;
                        }
                        output.WriteLine(
                            item.FullFilePath
                            + "\t" + item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorMonoisotopicPeakMz.ToString() //CultureInfo.InvariantCulture
                            + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.CrossType.ToString()
                            + "\t"
                            + "\t" + item.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.BaseSequence + position
                            + "\t" + item.FullSequence 
                            + "\t" + (item.PeptideMonisotopicMass.HasValue ? item.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                            + "\t" + item.Score.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.XLBestScore.ToString(CultureInfo.InvariantCulture)
                            + "\t" + (item.XlRank != null ? item.XlRank[0].ToString(CultureInfo.InvariantCulture) : "-")
                        );
                    }
                    else
                    {
                        output.WriteLine(
                                                item.FullFilePath
                                                + "\t" + item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.ScanPrecursorMonoisotopicPeakMz.ToString() //CultureInfo.InvariantCulture
                                                + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.CrossType.ToString()
                                                + "\t"
                                                + "\t" + item.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                                       + "(" + (item.XlProteinPos).ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.BaseSequence + "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.FullSequence
                                                + "\t" + (item.PeptideMonisotopicMass.HasValue ? item.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                                                + "\t" + item.Score.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.XLBestScore.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.XlRank[0].ToString(CultureInfo.InvariantCulture)
                                                + "\t"
                                                + "\t" + item.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                                       + "(" + (item.XlProteinPos).ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.BetaPsmCross.BaseSequence + "(" + item.BetaPsmCross.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.BetaPsmCross.FullSequence
                                                + "\t" + (item.BetaPsmCross.PeptideMonisotopicMass.HasValue ? item.BetaPsmCross.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                                                + "\t" + item.BetaPsmCross.Score.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.BetaPsmCross.XLBestScore.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.XlRank[1].ToString(CultureInfo.InvariantCulture)
                                                + "\t"
                                                + "\t" + item.XLTotalScore.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + ((item.PeptideMonisotopicMass.HasValue && item.BetaPsmCross.PeptideMonisotopicMass.HasValue) ? (item.BetaPsmCross.ScanPrecursorMass - item.BetaPsmCross.PeptideMonisotopicMass.Value - item.PeptideMonisotopicMass.Value).ToString(CultureInfo.InvariantCulture) : "---")
                                                + "\t" + (item.FdrInfo != null ? item.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture) : "-")
                                                + "\t" + item.ParentIonExist + "." + item.BetaPsmCross.ParentIonExist
                                                + "\t" + (item.Charge2IonExist + item.BetaPsmCross.Charge2IonExist).ToString(CultureInfo.InvariantCulture)
                                                );
                    }
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        public void WriteSingleToTsv(List<PsmCross> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("File Name\tScan Numer\tPrecusor MZ\tPrecusor charge\tPrecusor mass\tCross-link type" +
                    "\tPep1\tPep1 Protein Access(Protein link site)\tPep1 Base sequence(crosslink site)\tPep1 Full sequence\tPep1 mass\tPep1 Score\tPep1 XLTotalScore\tPep1 Rank" +
                    "\tQValue");
                foreach (var item in items)
                {
                    string position = "";
                    switch (item.CrossType)
                    {
                        case PsmCrossType.Singe:
                            break;
                        case PsmCrossType.Loop:
                            position = "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + "-" + item.XlPos2.ToString(CultureInfo.InvariantCulture) + ")";
                            break;
                        default:
                            position = "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + ")";
                            break;                           
                    }
                    output.WriteLine(
                        item.FullFilePath
                        + "\t" + item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.ScanPrecursorMonoisotopicPeakMz.ToString() //CultureInfo.InvariantCulture
                        + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.CrossType.ToString()
                        + "\t"
                        + "\t" + item.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.BaseSequence + position
                        + "\t" + item.FullSequence 
                        + "\t" + (item.PeptideMonisotopicMass.HasValue ? item.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                        + "\t" + item.Score.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.XLTotalScore.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (item.XlRank != null ? item.XlRank[0].ToString(CultureInfo.InvariantCulture) : "-")
                        + "\t" + (item.FdrInfo != null ? item.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture) : "-")
                    );
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        public void WriteCrosslinkToTxtForPercolator(List<PsmCross> items, string outputFolder, string fileName, CrosslinkerTypeClass crosslinker, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum" +
                    "\tPeptide\tProtein");
                foreach (var item in items)
                {
                    if (item.BaseSequence!=null && item.BetaPsmCross.BaseSequence != null && item.ProteinAccesion != null && item.BetaPsmCross.ProteinAccesion!= null)
                    {
                        string x = "T"; int label = 1;
                        if (item.IsDecoy || item.BetaPsmCross.IsDecoy)
                        {
                            x = "D"; label = -1;
                        }
                        output.WriteLine(
                            x + "-" + item.ScanNumber.ToString(CultureInfo.InvariantCulture) + "-" + item.ScanRetentionTime.ToString(CultureInfo.InvariantCulture)
                            + "\t" + label.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.XLTotalScore.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.DScore.ToString(CultureInfo.InvariantCulture)
                            + "\t" + (item.XlRank[0] + item.XlRank[1]).ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                            + "\t" + ((item.PeptideMonisotopicMass.HasValue && item.BetaPsmCross.PeptideMonisotopicMass.HasValue) ? ((item.ScanPrecursorMass - item.BetaPsmCross.PeptideMonisotopicMass.Value - item.PeptideMonisotopicMass.Value - crosslinker.TotalMass) / item.ScanPrecursorMass * 1E6).ToString(CultureInfo.InvariantCulture) : "---")
                            + "\t" + item.BetaPsmCross.BaseSequence.Length.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.BaseSequence.Length.ToString(CultureInfo.InvariantCulture)
                            + "\t" + (item.BetaPsmCross.BaseSequence.Length + item.BaseSequence.Length).ToString(CultureInfo.InvariantCulture)
                            + "\t" + "-." + item.BaseSequence + item.XlPos.ToString(CultureInfo.InvariantCulture) + "--" + item.BetaPsmCross.BaseSequence + item.BetaPsmCross.XlPos.ToString(CultureInfo.InvariantCulture) + ".-"
                            + "\t" + item.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                   + "(" + item.XlProteinPos.ToString(CultureInfo.InvariantCulture) + ")"
                            + "\t" + item.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                   + "(" + item.BetaPsmCross.XlProteinPos.ToString(CultureInfo.InvariantCulture) + ")"
                            );
                    }        
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        public void WritePepXML_xl(List<PsmCross> items, List<Protein> proteinList, string databasePath, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<string> localizeableModificationTypes, string outputFolder, string fileName, List<string> nestedIds)
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(pepXML.Generated.msms_pipeline_analysis));
            var _pepxml = new pepXML.Generated.msms_pipeline_analysis();

            #region Add element to pepXML

            _pepxml.date = DateTime.Now;
            _pepxml.summary_xml = items[0].FullFilePath + ".pep.xml";

            string proteaseC = ""; string proteaseNC = "";
            foreach (var x in CommonParameters.DigestionParams.Protease.SequencesInducingCleavage) { proteaseC += x; }
            foreach (var x in CommonParameters.DigestionParams.Protease.SequencesPreventingCleavage) { proteaseNC += x; }

            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass().SelectCrosslinker(XlSearchParameters.CrosslinkerType);

            uint proteinTot = Convert.ToUInt32(proteinList.Count);

            string fileNameNoExtension = Path.GetFileNameWithoutExtension(items[0].FullFilePath);
            string filePathNoExtension = Path.ChangeExtension(items[0].FullFilePath, null);

            var para = new List<pepXML.Generated.nameValueType>();
            {
                para.Add(new pepXML.Generated.nameValueType { name = "threads", value = "" });
                para.Add(new pepXML.Generated.nameValueType { name = "database", value = databasePath });
                para.Add(new pepXML.Generated.nameValueType { name = "MS_data_file", value = items[0].FullFilePath });

                para.Add(new pepXML.Generated.nameValueType { name = "Search with All Possible Beta Peptides", value = XlSearchParameters.CrosslinkSearchWithAllBeta.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-link Precusor Mass Tolence", value = XlSearchParameters.XlPrecusorMsTl.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker Type", value = crosslinker.CrosslinkerName });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker mass", value = crosslinker.TotalMass.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker cleavable", value = crosslinker.Cleavable.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker cleavable long mass", value = crosslinker.CleaveMassLong.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker cleavable short mass", value = crosslinker.CleaveMassShort.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Cross-linker xl site", value = crosslinker.CrosslinkerModSites.ToString() });

                para.Add(new pepXML.Generated.nameValueType { name = "Generate decoy proteins", value = XlSearchParameters.DecoyType.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "MaxMissed Cleavages", value = CommonParameters.DigestionParams.MaxMissedCleavages.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Protease", value = CommonParameters.DigestionParams.Protease.Name });
                para.Add(new pepXML.Generated.nameValueType { name = "Initiator Methionine", value = CommonParameters.DigestionParams.InitiatorMethionineBehavior.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Modification Isoforms", value = CommonParameters.DigestionParams.MaxModificationIsoforms.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Min Peptide Len", value = CommonParameters.DigestionParams.MinPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Max Peptide Len", value = CommonParameters.DigestionParams.MaxPeptideLength.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Product Mass Tolerance", value = CommonParameters.ProductMassTolerance.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Ions to search", value = "B " + CommonParameters.BIons.ToString() + " Y " + CommonParameters.YIons.ToString() + " C " + CommonParameters.CIons.ToString() + " Z " + CommonParameters.ZdotIons.ToString() });
                para.Add(new pepXML.Generated.nameValueType { name = "Allowed Beta Precusor Mass Difference", value = XlSearchParameters.XlBetaPrecusorMsTl.ToString() });
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
                         //search_engine = pepXML.Generated.engineType.Kojak,
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
                int modsFixedNum = items[i].CompactPeptides.First().Value.Item2.First().allModsOneIsNterminus.Count;
                var mods = new List<pepXML.Generated.modInfoDataTypeMod_aminoacid_mass> ();
                for (int j = 0; j < modsFixedNum; j++)
                {
                    var mod = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass
                    {
                        mass = items[i].CompactPeptides.First().Value.Item2.First().allModsOneIsNterminus.Values.Select(p => p.monoisotopicMass).ToList()[j],
                        position = (items[i].CompactPeptides.First().Value.Item2.First().allModsOneIsNterminus.Keys.ToList()[j]-1).ToString()
                    };
                    mods.Add(mod);
                }

                if (items[i].CrossType == PsmCrossType.Singe)
                {
                    var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                    {
                        hit_rank = 1,
                        peptide = items[i].BaseSequence,
                        peptide_prev_aa = items[i].CompactPeptides.First().Value.Item2.First().PreviousAminoAcid.ToString(),
                        peptide_next_aa = items[i].CompactPeptides.First().Value.Item2.First().NextAminoAcid.ToString(),
                        protein = items[i].CompactPeptides.First().Value.Item2.First().Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMonoisotopicPeakMz * items[i].ScanPrecursorCharge,
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
                if (items[i].CrossType == PsmCrossType.DeadEnd || items[i].CrossType == PsmCrossType.DeadEndH2O || items[i].CrossType == PsmCrossType.DeadEndNH2 || items[i].CrossType == PsmCrossType.DeadEndTris)
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
                    var mod = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass { mass = crosslinkerDeadEndMass, position = items[i].XlPos.ToString() };
                    mods.Add(mod);
                    var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                    {
                        hit_rank = 1,
                        peptide = items[i].BaseSequence,
                        peptide_prev_aa = items[i].CompactPeptides.First().Value.Item2.First().PreviousAminoAcid.ToString(),
                        peptide_next_aa = items[i].CompactPeptides.First().Value.Item2.First().NextAminoAcid.ToString(),
                        protein = items[i].CompactPeptides.First().Value.Item2.First().Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMonoisotopicPeakMz * items[i].ScanPrecursorCharge,
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
                if (items[i].CrossType == PsmCrossType.Inter || items[i].CrossType == PsmCrossType.Intra || items[i].CrossType == PsmCrossType.Cross)
                {
                    int modsFixedNumBeta = items[i].BetaPsmCross.CompactPeptides.First().Value.Item2.First().allModsOneIsNterminus.Count;
                    var modsBeta = new List<pepXML.Generated.modInfoDataTypeMod_aminoacid_mass>();
                    for (int j = 0; j < modsFixedNumBeta; j++)
                    {
                        var modBeta = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass
                        {
                            mass = items[i].BetaPsmCross.CompactPeptides.First().Value.Item2.First().allModsOneIsNterminus.Values.Select(p => p.monoisotopicMass).ToList()[j],
                            position = items[i].BetaPsmCross.CompactPeptides.First().Value.Item2.First().allModsOneIsNterminus.Keys.ToList()[j].ToString()
                        };
                        modsBeta.Add(modBeta);
                    }

                    var alpha = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                    {
                        peptide = items[i].BaseSequence,
                        peptide_prev_aa = items[i].CompactPeptides.First().Value.Item2.First().PreviousAminoAcid.ToString(),
                        peptide_next_aa = items[i].CompactPeptides.First().Value.Item2.First().NextAminoAcid.ToString(),
                        protein = items[i].CompactPeptides.First().Value.Item2.First().Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].PeptideMonisotopicMass.Value,
                        complement_mass = (float)(items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value),
                        designation = "alpha",
                        modification_info = new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = mods.ToArray() },
                        xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].XLBestScore.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].XlPos.ToString() },
                                                }
                    };
                    var beta = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                    {
                        peptide = items[i].BetaPsmCross.BaseSequence,
                        peptide_prev_aa = items[i].BetaPsmCross.CompactPeptides.First().Value.Item2.First().PreviousAminoAcid.ToString(),
                        peptide_next_aa = items[i].BetaPsmCross.CompactPeptides.First().Value.Item2.First().NextAminoAcid.ToString(),
                        protein = items[i].BetaPsmCross.CompactPeptides.First().Value.Item2.First().Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].BetaPsmCross.PeptideMonisotopicMass.Value,
                        complement_mass = (float)(items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value),
                        designation = "beta",
                        modification_info = new pepXML.Generated.modInfoDataType { mod_aminoacid_mass = modsBeta.ToArray() },
                        xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].BetaPsmCross.XLBestScore.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].BetaPsmCross.XlPos.ToString() },
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
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMonoisotopicPeakMz * items[i].ScanPrecursorCharge,
                        massdiff = (items[i].ScanPrecursorMass - items[i].BetaPsmCross.PeptideMonisotopicMass.Value - items[i].PeptideMonisotopicMass.Value - crosslinker.TotalMass).ToString(),
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
                if (items[i].CrossType == PsmCrossType.Loop)
                {
                    var thePeptide = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                    {
                        xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "link", value = items[i].XlPos.ToString() },
                                                    new pepXML.Generated.nameValueType{ name = "link", value = items[i].XlPos2.ToString() }
                                                }
                    };
                    var cross = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide[1] { thePeptide };
                    var searchHit = new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                    {
                        hit_rank = 1,
                        peptide = items[i].BaseSequence,
                        peptide_prev_aa = items[i].CompactPeptides.First().Value.Item2.First().PreviousAminoAcid.ToString(),
                        peptide_next_aa = items[i].CompactPeptides.First().Value.Item2.First().NextAminoAcid.ToString(),
                        protein = items[i].CompactPeptides.First().Value.Item2.First().Protein.Accession,
                        num_tot_proteins = 1,
                        calc_neutral_pep_mass = (float)items[i].ScanPrecursorMonoisotopicPeakMz * items[i].ScanPrecursorCharge,
                        massdiff = (items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value - crosslinker.LoopMass).ToString(),
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
                    precursor_neutral_mass = (float)items[i].ScanPrecursorMonoisotopicPeakMz * items[i].ScanPrecursorCharge,
                    assumed_charge = items[i].ScanPrecursorCharge.ToString(),
                    index = Convert.ToUInt32(i + 1),
                    retention_time_sec = (float)items[i].ScanRetentionTime,
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

            TextWriter writer = new StreamWriter(Path.Combine(outputFolder, fileName + ".pep.xml"));
            _indexedSerializer.Serialize(writer, _pepxml);
            writer.Close();
            SucessfullyFinishedWritingFile(Path.Combine(outputFolder, fileName + ".pep.xml"), nestedIds);
        }

        #endregion Private Methods
    }
}