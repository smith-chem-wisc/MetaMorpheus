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

        private void WriteCrosslinkToTsv(List<PsmCross> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("File Name\tScan Numer\tPrecusor MZ\tPrecusor charge\tPrecusor mass" +
                    "\tPep1\tPep1 Protein Access(Protein link site)\tPep1 Base sequence(crosslink site)\tPep1 Full sequence\tPep1 mass\tPep1 Score\tPep1 XLBestScore\tPep1 Rank" +
                    "\tPep2\tPep2 Protein Access(Protein link site)\tPep2 Base sequence(crosslink site)\tPep2 Full sequence\tPep2 mass\tPep2 Score\tPep2 XLBestScore\tPep2 Rank" +
                    "\tSummary\tTotalScore\tMass diff\tQValue\tParentIons\tCharge2Number\tLabel");
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
                                            + "\t"
                                            + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                                   + "(" + (item.MostProbableProteinInfo.PeptidesWithSetModifications.First().OneBasedStartResidueInProtein + item.XlPos - 1).ToString(CultureInfo.InvariantCulture) + ")"
                                            + "\t" + item.BaseSequence + "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                            + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence
                                            + "\t" + (item.PeptideMonisotopicMass.HasValue ? item.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                                            + "\t" + item.Score.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.XLBestScore.ToString(CultureInfo.InvariantCulture)
                                            + "\t" + item.XlRank[0].ToString(CultureInfo.InvariantCulture)
                                            + "\t"
                                            + "\t" + item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                                   + "(" + (item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().OneBasedStartResidueInProtein + item.BetaPsmCross.XlPos - 1).ToString(CultureInfo.InvariantCulture) + ")"
                                            + "\t" + item.BetaPsmCross.BaseSequence + "(" + item.BetaPsmCross.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                            + "\t" + item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence
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
                                            + "\t" + label
                                            );
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        private void WriteAllToTsv(List<PsmCross> items, string outputFolder, string fileName, List<string> nestedIds)
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
                    if (item.CrossType != PsmCrossType.Cross)
                    {
                        output.WriteLine(
                            item.FullFilePath
                            + "\t" + item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorMonoisotopicPeakMz.ToString() //CultureInfo.InvariantCulture
                            + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                            + "\t" + item.CrossType.ToString()
                            + "\t"
                            + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                            + "(" + (item.MostProbableProteinInfo.PeptidesWithSetModifications.First().OneBasedStartResidueInProtein + item.XlPos - 1).ToString(CultureInfo.InvariantCulture) + ")"
                            + "\t" + item.BaseSequence + "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                            + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence
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
                                                + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                                       + "(" + (item.MostProbableProteinInfo.PeptidesWithSetModifications.First().OneBasedStartResidueInProtein + item.XlPos - 1).ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.BaseSequence + "(" + item.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence
                                                + "\t" + (item.PeptideMonisotopicMass.HasValue ? item.PeptideMonisotopicMass.Value.ToString(CultureInfo.InvariantCulture) : "---")
                                                + "\t" + item.Score.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.XLBestScore.ToString(CultureInfo.InvariantCulture)
                                                + "\t" + item.XlRank[0].ToString(CultureInfo.InvariantCulture)
                                                + "\t"
                                                + "\t" + item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                                                       + "(" + (item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().OneBasedStartResidueInProtein + item.BetaPsmCross.XlPos - 1).ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.BetaPsmCross.BaseSequence + "(" + item.BetaPsmCross.XlPos.ToString(CultureInfo.InvariantCulture) + ")"
                                                + "\t" + item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().Sequence
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

        private void WriteCrosslinkToTxtForCLMSVault(List<PsmCross> items, string outputFolder, string fileName, CrosslinkerTypeClass crosslinker, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("Scan Number\tRet Time\tObs Mass\tCharge\tPSM Mass\tPPM Error\tScore\tdScore\tPep.Diff." +
                    "\tPeptide #1\tLink #1\tProtein #1" +
                    "\tPeptide #2\tLink #2\tProtein #2\tLinker Mass");
                foreach (var item in items)
                {
                    output.WriteLine(
                        item.ScanNumber.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.ScanRetentionTime.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.ScanPrecursorMass.ToString() //CultureInfo.InvariantCulture
                        + "\t" + item.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                        + "\t" + ((item.PeptideMonisotopicMass.HasValue && item.BetaPsmCross.PeptideMonisotopicMass.HasValue) ? (item.BetaPsmCross.PeptideMonisotopicMass.Value + item.PeptideMonisotopicMass.Value + crosslinker.TotalMass).ToString(CultureInfo.InvariantCulture) : "---")
                        + "\t" + ((item.PeptideMonisotopicMass.HasValue && item.BetaPsmCross.PeptideMonisotopicMass.HasValue) ? ((item.ScanPrecursorMass - item.BetaPsmCross.PeptideMonisotopicMass.Value - item.PeptideMonisotopicMass.Value - crosslinker.TotalMass) / item.ScanPrecursorMass * 10E6).ToString(CultureInfo.InvariantCulture) : "---")
                        + "\t" + item.XLTotalScore.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (item.Score + item.BetaPsmCross.Score).ToString(CultureInfo.InvariantCulture)
                        + "\t" + ((item.PeptideMonisotopicMass.HasValue && item.BetaPsmCross.PeptideMonisotopicMass.HasValue) ? (item.ScanPrecursorMass - item.BetaPsmCross.PeptideMonisotopicMass.Value - item.PeptideMonisotopicMass.Value - crosslinker.TotalMass).ToString(CultureInfo.InvariantCulture) : "---")

                        + "\t" + item.BaseSequence
                        + "\t" + item.XlPos.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)

                        + "\t" + item.BetaPsmCross.BaseSequence
                        + "\t" + item.BetaPsmCross.XlPos.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        + "\t" + crosslinker.TotalMass.ToString(CultureInfo.InvariantCulture)
                        );
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        private void WriteCrosslinkToTxtForPercolator(List<PsmCross> items, string outputFolder, string fileName, CrosslinkerTypeClass crosslinker, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum" +
                    "\tPeptide\tProtein");
                foreach (var item in items)
                {
                    string x = "T"; string label = "1";
                    if (item.IsDecoy || item.BetaPsmCross.IsDecoy)
                    {
                        x = "D"; label = "-1";
                    }
                    output.WriteLine(
                        x + "-" + item.ScanNumber.ToString(CultureInfo.InvariantCulture) + "-" + item.ScanRetentionTime.ToString(CultureInfo.InvariantCulture)
                        + "\t" + label
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
                        + "\t" + item.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        );
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        private void WritePepXML_xl(List<PsmCross> items, List<DbForTask> dbFilenameList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ModificationWithMass> localizeableModifications, string outputFolder, string fileName, List<string> nestedIds)
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
            string modsFixed = ""; string modsVar = "";
            foreach (var x in CommonParameters.ListOfModsFixed) { modsFixed += x.Item2 + "."; }
            foreach (var x in CommonParameters.ListOfModsVariable) { modsVar += x.Item2 + "."; }

            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, XlSearchParameters.DecoyType, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();

            uint proteinTot = Convert.ToUInt32(proteinList.Count);

            string fileNameNoExtension = Path.GetFileNameWithoutExtension(items[0].FullFilePath);
            string filePathNoExtension = Path.ChangeExtension(items[0].FullFilePath, null);

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
                             enzyme = CommonParameters.DigestionParams.Protease.Name,
                             max_num_internal_cleavages = CommonParameters.DigestionParams.MaxMissedCleavages.ToString(),
                             //min_number_termini = "2"
                         },
                         parameter = new pepXML.Generated.nameValueType[]
                         {
                             new pepXML.Generated.nameValueType{ name = "threads", value = "" },
                             new pepXML.Generated.nameValueType{ name = "database", value = dbFilenameList[0].FilePath },
                             new pepXML.Generated.nameValueType{ name = "MS_data_file", value = items[0].FullFilePath },

                             new pepXML.Generated.nameValueType{ name = "Search with All Possible Beta Peptides", value = XlSearchParameters.CrosslinkSearchWithAllBeta.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-link Precusor Mass Tolence", value = XlSearchParameters.XlPrecusorMsTl.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker Type", value = crosslinker.CrosslinkerName },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker mass", value = crosslinker.TotalMass.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker cleavable", value = crosslinker.Cleavable.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker cleavable long mass", value = crosslinker.CleaveMassLong.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker cleavable short mass", value = crosslinker.CleaveMassShort.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Cross-linker xl site", value = crosslinker.CrosslinkerModSite.ToString() },

                             new pepXML.Generated.nameValueType{ name = "Generate decoy proteins", value = XlSearchParameters.DecoyType.ToString() },
                             new pepXML.Generated.nameValueType{ name = "MaxMissed Cleavages", value = CommonParameters.DigestionParams.MaxMissedCleavages.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Protease", value = CommonParameters.DigestionParams.Protease.Name },
                             new pepXML.Generated.nameValueType{ name = "Initiator Methionine", value = CommonParameters.DigestionParams.InitiatorMethionineBehavior.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Max Modification Isoforms", value = CommonParameters.DigestionParams.MaxModificationIsoforms.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Min Peptide Len", value = CommonParameters.DigestionParams.MinPeptideLength.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Max Peptide Len", value = CommonParameters.DigestionParams.MaxPeptideLength.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Product Mass Tolerance", value = CommonParameters.ProductMassTolerance.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Ions to search", value = "B "+ CommonParameters.BIons.ToString() + " Y " + CommonParameters.YIons.ToString() + " C " + CommonParameters.CIons.ToString() + " Z " + CommonParameters.ZdotIons.ToString() },
                             new pepXML.Generated.nameValueType{ name = "Allowed Beta Precusor Mass Difference", value = XlSearchParameters.XlBetaPrecusorMsTl.ToString()},

                             new pepXML.Generated.nameValueType{ name = "Fixed Modifications", value = modsFixed },
                             new pepXML.Generated.nameValueType{ name = "Variable Modificaions", value = modsVar },
                             new pepXML.Generated.nameValueType{ name = "Localize All Modifications", value = CommonParameters.LocalizeAll.ToString() },
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
                                new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                                {
                                    hit_rank = 1,
                                    peptide = "-", peptide_prev_aa="-", peptide_next_aa="-", protein="-", num_tot_proteins = 1,
                                    calc_neutral_pep_mass = (float)items[i].ScanPrecursorMonoisotopicPeakMz * items[i].ScanPrecursorCharge,
                                    massdiff = (items[i].ScanPrecursorMass - items[i].BetaPsmCross.PeptideMonisotopicMass.Value - items[i].PeptideMonisotopicMass.Value - crosslinker.TotalMass).ToString(),
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
                                                peptide = items[i].BaseSequence,
                                                peptide_prev_aa = items[i].MostProbableProteinInfo.PeptidesWithSetModifications.First().PreviousAminoAcid.ToString(),
                                                peptide_next_aa = items[i].MostProbableProteinInfo.PeptidesWithSetModifications.First().NextAminoAcid.ToString(),
                                                protein = items[i].MostProbableProteinInfo.PeptidesWithSetModifications.First().Protein.Accession,
                                                num_tot_proteins = proteinTot/2,
                                                calc_neutral_pep_mass = (float)items[i].PeptideMonisotopicMass.Value,
                                                complement_mass = (float)(items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value),
                                                designation = "alpha",

                                                xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].XLBestScore.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].XlPos.ToString() },
                                                }
                                            },
                                            new pepXML.Generated.msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hitXlinkLinked_peptide
                                            {
                                                peptide = items[i].BetaPsmCross.BaseSequence,
                                                peptide_prev_aa = items[i].BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().PreviousAminoAcid.ToString(),
                                                peptide_next_aa = items[i].BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().NextAminoAcid.ToString(),
                                                protein = items[i].BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().Protein.Accession,
                                                num_tot_proteins = proteinTot,
                                                calc_neutral_pep_mass = (float)items[i].BetaPsmCross.PeptideMonisotopicMass.Value,
                                                complement_mass = (float)(items[i].ScanPrecursorMass - items[i].PeptideMonisotopicMass.Value),
                                                designation = "beta",

                                                xlink_score = new pepXML.Generated.nameValueType[]
                                                {
                                                    new pepXML.Generated.nameValueType{ name = "xlscore", value = items[i].BetaPsmCross.XLBestScore.ToString() },
                                                    new pepXML.Generated.nameValueType{name = "link", value = items[i].BetaPsmCross.XlPos.ToString() },
                                                }
                                            }
                                        }
                                    },
                                    search_score = new pepXML.Generated.nameValueType[]
                                    {
                                        new pepXML.Generated.nameValueType{ name = "xlTotalScore", value = items[i].XLTotalScore.ToString()},
                                        new pepXML.Generated.nameValueType{ name = "Qvalue", value = items[i].FdrInfo.QValue.ToString() }
                                    }
                                }
                            }
                        }
                    }
                };

                #region mods infomation

                int modsFixedNum1 = items[i].MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
                int modsFixedNum2 = items[i].BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
                if (modsFixedNum1 != 0)
                {
                    modsFixedNum1 = items[i].MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
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
                    if (modsFixedNum1 == 5)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[5]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum1 == 6)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[6]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    for (int j = 0; j < modsFixedNum1; j++)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info.mod_aminoacid_mass[j].mass = items[i].MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Values.Select(p => p.monoisotopicMass).ToList()[j];
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[0].modification_info.mod_aminoacid_mass[j].position = items[i].MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Keys.ToList()[j].ToString();
                    }
                }
                if (modsFixedNum2 != 0)
                {
                    modsFixedNum2 = items[i].BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Count;
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
                    if (modsFixedNum2 == 5)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[5]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    if (modsFixedNum2 == 6)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info = new pepXML.Generated.modInfoDataType
                        {
                            mod_aminoacid_mass = new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass[6]
                        {
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{},
                            new pepXML.Generated.modInfoDataTypeMod_aminoacid_mass{}
                        }
                        };
                    }
                    for (int j = 0; j < modsFixedNum2; j++)
                    {
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info.mod_aminoacid_mass[j].mass = items[i].BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Values.Select(p => p.monoisotopicMass).ToList()[j];
                        _pepxml.msms_run_summary[0].spectrum_query[i].search_result[0].search_hit[0].xlink.linked_peptide[1].modification_info.mod_aminoacid_mass[j].position = items[i].BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Keys.ToList()[j].ToString();
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