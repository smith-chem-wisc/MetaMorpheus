namespace OldInternalLogic
{
    public static class Exporters
    {
        //public static void WriteToTabDelimitedTextFile<T>(IEnumerable<T> items, string filepath)
        //{
        //    using (StreamWriter output = new StreamWriter(filepath))
        //    {
        //        if (typeof(T) == typeof(PeptideSpectrumMatch))
        //        {
        //            output.WriteLine(PeptideSpectrumMatch.Header);
        //        }
        //        else if (typeof(T) == typeof(IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>))
        //        {
        //            output.WriteLine(PeptideSpectrumMatch.Header + IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>.Header);
        //        }
        //        foreach (T item in items)
        //        {
        //            output.WriteLine(item.ToString());
        //        }
        //    }
        //}

        //public static void WriteMZIdentMLFile(string outputFilepath,
        //    IEnumerable<TandemMassSpectra> datas,
        //    int minimumAssumedPrecursorChargeState, int maximumAssumedPrecursorChargeState,
        //    int maximumNumberOfPeaks,
        //    bool assignChargeStates, bool deisotope, bool onTheFlyDecoys, int proteins,
        //    Protease protease, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior,
        //    int maximumVariableModificationIsoforms,
        //    MassTolerance precursorMassTolerance,
        //    IEnumerable<double> acceptedPrecursorMassErrors,
        //    MassTolerance productMassTolerance,
        //    double maximumFalseDiscoveryRate, string outputFolder,
        //    IEnumerable<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> psms,
        //    IEnumerable<Protein> allProteins)
        //{
        //    using (XmlTextWriter output = new XmlTextWriter(outputFilepath, Encoding.UTF8))
        //    {
        //        output.Formatting = Formatting.Indented;

        //        output.WriteStartDocument();

        //        output.WriteStartElement("MzIdentML");
        //        output.WriteAttributeString("xmlns", string.Empty, null, "http://psidev.info/psi/pi/mzIdentML/1.1.1");
        //        output.WriteAttributeString("xmlns", "xsi", null, "http://www.w3.org/2001/XMLSchema-instance");
        //        output.WriteAttributeString("xsi", "schemaLocation", "http://www.w3.org/2001/XMLSchema-instance", "http://psidev.info/psi/pi/mzIdentML/1.1.1 https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.1.1.xsd");
        //        output.WriteAttributeString("creationDate", DateTime.Now.ToString("s"));
        //        output.WriteAttributeString("id", Path.GetFileNameWithoutExtension(outputFilepath));
        //        output.WriteAttributeString("version", "1.1.1");

        //        output.WriteStartElement("cvList");
        //        output.WriteStartElement("cv");
        //        output.WriteAttributeString("id", "PSI-MS");
        //        output.WriteAttributeString("fullName", "Proteomics Standards Initiative Mass Spectrometry Vocabularies");
        //        output.WriteAttributeString("uri", "http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
        //        output.WriteAttributeString("version", "2.25.0");
        //        output.WriteEndElement();  // cv
        //        output.WriteStartElement("cv");
        //        output.WriteAttributeString("id", "UO");
        //        output.WriteAttributeString("fullName", "UNIT-ONTOLOGY");
        //        output.WriteAttributeString("uri", "http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
        //        output.WriteEndElement();  // cv
        //        output.WriteStartElement("cv");
        //        output.WriteAttributeString("id", "UNIMOD");
        //        output.WriteAttributeString("fullName", "UNIMOD");
        //        output.WriteAttributeString("uri", "http://www.unimod.org/obo/unimod.obo");
        //        output.WriteEndElement();  // cv
        //        output.WriteStartElement("cv");
        //        output.WriteAttributeString("id", "PSI-MOD");
        //        output.WriteAttributeString("fullName", "PSI-MOD");
        //        output.WriteAttributeString("uri", "http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/mod/data/PSI-MOD.obo");
        //        output.WriteEndElement();  // cv
        //        output.WriteEndElement();  // cvList

        //        output.WriteStartElement("AnalysisSoftwareList");
        //        output.WriteStartElement("AnalysisSoftware");
        //        output.WriteAttributeString("id", "AS_MetaMorpheus");
        //        output.WriteAttributeString("name", "MetaMorpheus");
        //        output.WriteStartElement("SoftwareName");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("name", "MetaMorpheus");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // SoftwareName
        //        output.WriteEndElement();  // AnalysisSoftware
        //        output.WriteEndElement();  // AnalysisSoftwareList

        //        output.WriteStartElement("SequenceCollection");
        //        foreach (Protein protein in allProteins)
        //        {
        //            output.WriteStartElement("DBSequence");
        //            string accession = protein.Accession;
        //            output.WriteAttributeString("id", "DBS_" + accession);
        //            output.WriteAttributeString("accession", accession);
        //            //output.WriteAttributeString("searchDatabase_ref", "SDB_" + Path.GetFileNameWithoutExtension(proteomeDatabaseFilepath));
        //            output.WriteAttributeString("length", protein.Length.ToString());
        //            output.WriteEndElement();  // DBSequence
        //        }
        //        int index = 1;
        //        Dictionary<PeptideSpectrumMatch, int> ids = new Dictionary<PeptideSpectrumMatch, int>();
        //        foreach (IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch> psm_with_fdr in psms)
        //        {
        //            PeptideSpectrumMatch psm = psm_with_fdr.Identification;
        //            PeptideWithSetModifications peptide = psm.Peptide;
        //            output.WriteStartElement("Peptide");
        //            output.WriteAttributeString("id", "P_" + index.ToString());
        //            output.WriteStartElement("PeptideSequence");
        //            output.WriteString(peptide.BaseSequence);
        //            output.WriteEndElement();  // PeptideSequence
        //            foreach (KeyValuePair<int, MorpheusModification> variable_mod_kvp in peptide.Modifications)
        //            {
        //                MorpheusModification variable_mod = variable_mod_kvp.Value;
        //                output.WriteStartElement("Modification");
        //                int location;
        //                int residue_index;
        //                ConvertModificationIndex(peptide, variable_mod_kvp.Key, out location, out residue_index);
        //                output.WriteAttributeString("location", location.ToString());
        //                output.WriteAttributeString("monoisotopicMassDelta", variable_mod.MonoisotopicMassShift.ToString());
        //                output.WriteAttributeString("residues", peptide[residue_index].ToString());
        //                output.WriteStartElement("cvParam");
        //                output.WriteAttributeString("accession", variable_mod.Database);
        //                output.WriteAttributeString("name", variable_mod.DatabaseName);
        //                output.WriteAttributeString("cvRef", variable_mod.Database);
        //                output.WriteEndElement();  // cvParam
        //                output.WriteEndElement();  // Modification
        //            }
        //            output.WriteEndElement();  // Peptide
        //            ids.Add(psm, index);
        //            index++;
        //        }
        //        index = 1;
        //        foreach (IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch> psm_with_fdr in psms)
        //        {
        //            PeptideSpectrumMatch psm = psm_with_fdr.Identification;
        //            PeptideWithSetModifications peptide = psm.Peptide;
        //            output.WriteStartElement("PeptideEvidence");
        //            output.WriteAttributeString("id", "PE_" + index.ToString());
        //            output.WriteAttributeString("peptide_ref", "P_" + index.ToString());
        //            output.WriteAttributeString("dBSequence_ref", "DBS_" + peptide.protein.Accession);
        //            output.WriteAttributeString("isDecoy", psm.isDecoy.ToString().ToLower());
        //            output.WriteAttributeString("start", peptide.OneBasedStartResidueInProtein.ToString());
        //            output.WriteAttributeString("end", peptide.OneBasedEndResidueInProtein.ToString());
        //            output.WriteAttributeString("pre", peptide.PreviousAminoAcid.ToString());
        //            output.WriteAttributeString("post", peptide.NextAminoAcid.ToString());
        //            output.WriteEndElement();  // PeptideEvidence
        //            index++;
        //        }
        //        output.WriteEndElement();  // SequenceCollection

        //        output.WriteStartElement("AnalysisCollection");
        //        output.WriteStartElement("SpectrumIdentification");
        //        output.WriteAttributeString("id", "SI");
        //        output.WriteAttributeString("spectrumIdentificationList_ref", "SIL");
        //        output.WriteAttributeString("spectrumIdentificationProtocol_ref", "SIP");
        //        Dictionary<string, string> spectral_data_ids = new Dictionary<string, string>();
        //        index = 1;
        //        foreach (var data in datas)
        //        {
        //            string spectral_data_id = "SD_" + index.ToString();
        //            spectral_data_ids.Add(data.filename, spectral_data_id);
        //            output.WriteStartElement("InputSpectra");
        //            output.WriteAttributeString("spectraData_ref", spectral_data_id);
        //            output.WriteEndElement();  // InputSpectra
        //            index++;
        //        }
        //        output.WriteStartElement("SearchDatabaseRef");
        //        //output.WriteAttributeString("searchDatabase_ref", "SDB_" + Path.GetFileNameWithoutExtension(proteomeDatabaseFilepath));
        //        output.WriteEndElement();  // SearchDatabaseRef
        //        output.WriteEndElement();  // SpectrumIdentification
        //        output.WriteEndElement();  // AnalysisCollection

        //        output.WriteStartElement("AnalysisProtocolCollection");
        //        output.WriteStartElement("SpectrumIdentificationProtocol");
        //        output.WriteAttributeString("id", "SIP");
        //        output.WriteAttributeString("analysisSoftware_ref", "AS_MetaMorpheus");
        //        output.WriteStartElement("SearchType");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001083");
        //        output.WriteAttributeString("name", "ms-ms search");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // SearchType
        //        output.WriteStartElement("AdditionalSearchParams");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001211");
        //        output.WriteAttributeString("name", "parent mass type mono");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001256");
        //        output.WriteAttributeString("name", "fragment mass type mono");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteStartElement("userParam");
        //        output.WriteAttributeString("name", "Minimum Assumed Precursor Charge State");
        //        output.WriteAttributeString("value", minimumAssumedPrecursorChargeState.ToString());
        //        output.WriteEndElement();  // userParam
        //        output.WriteStartElement("userParam");
        //        output.WriteAttributeString("name", "Maximum Assumed Precursor Charge State");
        //        output.WriteAttributeString("value", maximumAssumedPrecursorChargeState.ToString());
        //        output.WriteEndElement();  // userParam
        //        output.WriteStartElement("userParam");
        //        output.WriteAttributeString("name", "Maximum Number of MS/MS Peaks");
        //        output.WriteAttributeString("value", maximumNumberOfPeaks >= 0 ? maximumNumberOfPeaks.ToString() : "disabled");
        //        output.WriteEndElement();  // userParam
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1000034");
        //        output.WriteAttributeString("name", "charge deconvolution");
        //        output.WriteAttributeString("value", assignChargeStates.ToString().ToLower());
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1000033");
        //        output.WriteAttributeString("name", "deisotoping");
        //        output.WriteAttributeString("value", deisotope.ToString().ToLower());
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteStartElement("userParam");
        //        output.WriteAttributeString("name", "On-the-fly Decoy Protein Generation");
        //        output.WriteAttributeString("value", onTheFlyDecoys.ToString().ToLower());
        //        output.WriteEndElement();  // userParam
        //        output.WriteStartElement("userParam");
        //        output.WriteAttributeString("name", "Initiator Methionine Behavior");
        //        output.WriteAttributeString("value", initiatorMethionineBehavior.ToString().ToLower());
        //        output.WriteEndElement();  // userParam
        //        output.WriteStartElement("userParam");
        //        output.WriteAttributeString("name", "Maximum Number of Variable Modification Isoforms per Peptide");
        //        output.WriteAttributeString("value", maximumVariableModificationIsoforms.ToString());
        //        output.WriteEndElement();  // userParam
        //        output.WriteStartElement("userParam");
        //        output.WriteAttributeString("name", "Accepted Precursor Mass Errors");
        //        output.WriteAttributeString("value", string.Join("; ", acceptedPrecursorMassErrors));
        //        output.WriteAttributeString("unitAccession", "UO:0000221");
        //        output.WriteAttributeString("unitName", "dalton");
        //        output.WriteAttributeString("unitCvRef", "UO");
        //        output.WriteEndElement();  // userParam
        //        output.WriteEndElement();  // AdditionalSearchParams
        //        output.WriteStartElement("ModificationParams");
        //        //foreach (Modification fixed_mod in allMods)
        //        //{
        //        //    output.WriteStartElement("SearchModification");
        //        //    output.WriteAttributeString("massDelta", fixed_mod.MonoisotopicMassShift.ToString());
        //        //    output.WriteAttributeString("residues", fixed_mod.AminoAcid == char.MinValue ? "." : fixed_mod.AminoAcid.ToString());
        //        //    output.WriteStartElement("cvParam");
        //        //    output.WriteAttributeString("accession", fixed_mod.Database + ':' + fixed_mod.DatabaseAccessionNumber.ToString());
        //        //    output.WriteAttributeString("name", fixed_mod.DatabaseName);
        //        //    output.WriteAttributeString("cvRef", fixed_mod.Database);
        //        //    output.WriteEndElement();  // cvParam
        //        //    output.WriteEndElement();  // SearchModification
        //        //}
        //        output.WriteEndElement();  // ModificationParams
        //        output.WriteStartElement("Enzymes");
        //        output.WriteStartElement("Enzyme");
        //        output.WriteAttributeString("id", "E");
        //        output.WriteAttributeString("name", protease.Name);
        //        output.WriteAttributeString("semiSpecific", (protease.CleavageSpecificity == CleavageSpecificity.Semi || protease.CleavageSpecificity == CleavageSpecificity.SemiN || protease.CleavageSpecificity == CleavageSpecificity.SemiC).ToString().ToLower());
        //        output.WriteAttributeString("missedCleavages", maximumMissedCleavages.ToString());
        //        if (!string.IsNullOrWhiteSpace(protease.SiteRegexp))
        //        {
        //            output.WriteStartElement("SiteRegexp");
        //            output.WriteCData(protease.SiteRegexp);
        //            output.WriteEndElement();  // SiteRegexp
        //        }
        //        if (!string.IsNullOrWhiteSpace(protease.PsiMsAccessionNumber))
        //        {
        //            output.WriteStartElement("EnzymeName");
        //            output.WriteStartElement("cvParam");
        //            output.WriteAttributeString("accession", protease.PsiMsAccessionNumber);
        //            output.WriteAttributeString("name", protease.PsiMsName);
        //            output.WriteAttributeString("cvRef", "PSI-MS");
        //            output.WriteEndElement();  // cvParam
        //            output.WriteEndElement();  // EnzymeName
        //        }
        //        output.WriteEndElement();  // Enzyme
        //        output.WriteEndElement();  // Enzymes
        //        output.WriteStartElement("FragmentTolerance");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001412");
        //        output.WriteAttributeString("name", "search tolerance plus value");
        //        output.WriteAttributeString("value", productMassTolerance.Value.ToString());
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        if (productMassTolerance.Units == MassToleranceUnits.ppm)
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000169");
        //            output.WriteAttributeString("unitName", "parts per million");
        //        }
        //        else
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000221");
        //            output.WriteAttributeString("unitName", "dalton");
        //        }
        //        output.WriteAttributeString("unitCvRef", "UO");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001413");
        //        output.WriteAttributeString("name", "search tolerance minus value");
        //        output.WriteAttributeString("value", productMassTolerance.Value.ToString());
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        if (productMassTolerance.Units == MassToleranceUnits.ppm)
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000169");
        //            output.WriteAttributeString("unitName", "parts per million");
        //        }
        //        else
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000221");
        //            output.WriteAttributeString("unitName", "dalton");
        //        }
        //        output.WriteAttributeString("unitCvRef", "UO");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // FragmentTolerance
        //        output.WriteStartElement("ParentTolerance");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001412");
        //        output.WriteAttributeString("name", "search tolerance plus value");
        //        output.WriteAttributeString("value", precursorMassTolerance.Value.ToString());
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        if (precursorMassTolerance.Units == MassToleranceUnits.ppm)
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000169");
        //            output.WriteAttributeString("unitName", "parts per million");
        //        }
        //        else
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000221");
        //            output.WriteAttributeString("unitName", "dalton");
        //        }
        //        output.WriteAttributeString("unitCvRef", "UO");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001413");
        //        output.WriteAttributeString("name", "search tolerance minus value");
        //        output.WriteAttributeString("value", precursorMassTolerance.Value.ToString());
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        if (precursorMassTolerance.Units == MassToleranceUnits.ppm)
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000169");
        //            output.WriteAttributeString("unitName", "parts per million");
        //        }
        //        else
        //        {
        //            output.WriteAttributeString("unitAccession", "UO:0000221");
        //            output.WriteAttributeString("unitName", "dalton");
        //        }
        //        output.WriteAttributeString("unitCvRef", "UO");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // ParentTolerance
        //        output.WriteStartElement("Threshold");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001448");
        //        output.WriteAttributeString("name", "pep:FDR threshold");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteAttributeString("value", maximumFalseDiscoveryRate.ToString());
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // Threshold
        //        output.WriteEndElement();  // SpectrumIdentificationProtocol
        //        output.WriteStartElement("ProteinDetectionProtocol");
        //        output.WriteAttributeString("id", "PDP");
        //        output.WriteAttributeString("analysisSoftware_ref", "AS_MetaMorpheus");
        //        output.WriteStartElement("Threshold");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001447");
        //        output.WriteAttributeString("name", "prot:FDR threshold");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteAttributeString("value", maximumFalseDiscoveryRate.ToString());
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // Threshold
        //        output.WriteEndElement();  // ProteinDetectionProtocol
        //        output.WriteEndElement();  // AnalysisProtocolCollection

        //        output.WriteStartElement("DataCollection");
        //        output.WriteStartElement("Inputs");
        //        output.WriteStartElement("SearchDatabase");
        //        //output.WriteAttributeString("id", "SDB_" + Path.GetFileNameWithoutExtension(proteomeDatabaseFilepath));
        //        //output.WriteAttributeString("location", new Uri(proteomeDatabaseFilepath).AbsoluteUri);
        //        output.WriteAttributeString("numDatabaseSequences", proteins.ToString());
        //        output.WriteStartElement("FileFormat");
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1002660");
        //        output.WriteAttributeString("name", "UniProtKB XML sequence format");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // FileFormat
        //        //output.WriteStartElement("DatabaseName");
        //        ////output.WriteStartElement("userParam");
        //        ////output.WriteAttributeString("name", string.Join(",", proteomeDatabaseFilepaths));
        //        ////output.WriteEndElement();  // userParam
        //        //output.WriteEndElement();  // DatabaseName
        //        output.WriteStartElement("cvParam");
        //        output.WriteAttributeString("accession", "MS:1001073");
        //        output.WriteAttributeString("name", "database type amino acid");
        //        output.WriteAttributeString("cvRef", "PSI-MS");
        //        output.WriteEndElement();  // cvParam
        //        output.WriteEndElement();  // SearchDatabase
        //        foreach (var spectra in datas)
        //        {
        //            string data_filepath = spectra.filename;
        //            output.WriteStartElement("SpectraData");
        //            output.WriteAttributeString("id", spectral_data_ids[data_filepath]);
        //            output.WriteAttributeString("location", new Uri(data_filepath).AbsoluteUri);
        //            output.WriteStartElement("FileFormat");
        //            output.WriteStartElement("cvParam");
        //            if (Path.GetExtension(data_filepath).ToLower() == ".raw")
        //            {
        //                output.WriteAttributeString("accession", "MS:1000563");
        //                output.WriteAttributeString("name", "Thermo RAW format");
        //            }
        //            else if (Path.GetExtension(data_filepath).ToLower() == ".d")
        //            {
        //                output.WriteAttributeString("accession", "MS:1001509");
        //                output.WriteAttributeString("name", "Agilent MassHunter format");
        //            }
        //            else
        //            {
        //                output.WriteAttributeString("accession", "MS:1000584");
        //                output.WriteAttributeString("name", "mzML format");
        //            }
        //            output.WriteAttributeString("cvRef", "PSI-MS");
        //            output.WriteEndElement();  // cvParam
        //            output.WriteEndElement();  // FileFormat
        //            output.WriteStartElement("SpectrumIDFormat");
        //            output.WriteStartElement("cvParam");
        //            if (Path.GetExtension(data_filepath).ToLower() == ".raw")
        //            {
        //                output.WriteAttributeString("accession", "MS:1000768");
        //                output.WriteAttributeString("name", "Thermo nativeID format");
        //            }
        //            else if (Path.GetExtension(data_filepath).ToLower() == ".d")
        //            {
        //                output.WriteAttributeString("accession", "MS:1001508");
        //                output.WriteAttributeString("name", "Agilent MassHunter nativeID format");
        //            }
        //            else
        //            {
        //                output.WriteAttributeString("accession", "MS:1001530");
        //                output.WriteAttributeString("name", "mzML unique identifier");
        //            }
        //            output.WriteAttributeString("cvRef", "PSI-MS");
        //            output.WriteEndElement();  // cvParam
        //            output.WriteEndElement();  // SpectrumIDFormat
        //            output.WriteEndElement();  // SpectraData
        //        }
        //        output.WriteEndElement();  // Inputs
        //        output.WriteStartElement("AnalysisData");
        //        output.WriteStartElement("SpectrumIdentificationList");
        //        output.WriteAttributeString("id", "SIL");
        //        index = 1;
        //        foreach (IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch> psm_with_fdr in psms)
        //        {
        //            PeptideSpectrumMatch psm = psm_with_fdr.Identification;
        //            output.WriteStartElement("SpectrumIdentificationResult");
        //            output.WriteAttributeString("id", "SIR_" + index.ToString());
        //            output.WriteAttributeString("spectraData_ref", spectral_data_ids[psm.spectrumFilename]);
        //            output.WriteAttributeString("spectrumID", psm.spectrumID);
        //            output.WriteStartElement("SpectrumIdentificationItem");
        //            output.WriteAttributeString("id", "SII_" + index.ToString());
        //            output.WriteAttributeString("chargeState", psm.spectrumPrecursorCharge.ToString());
        //            output.WriteAttributeString("experimentalMassToCharge", psm.spectrumPrecursorMZ.ToString());
        //            output.WriteAttributeString("calculatedMassToCharge", MSPeak.MZFromMass(psm.Peptide.MonoisotopicMass, psm.spectrumPrecursorCharge).ToString());
        //            output.WriteAttributeString("passThreshold", (psm_with_fdr.QValue <= maximumFalseDiscoveryRate).ToString().ToLower());
        //            output.WriteAttributeString("rank", "1");
        //            output.WriteStartElement("PeptideEvidenceRef");
        //            output.WriteAttributeString("peptideEvidence_ref", "PE_" + index.ToString());
        //            output.WriteEndElement();  // PeptideEvidenceRef
        //            output.WriteStartElement("cvParam");
        //            output.WriteAttributeString("accession", "MS:1002662");
        //            output.WriteAttributeString("name", "MetaMorpheus:MetaMorpheus score");
        //            output.WriteAttributeString("cvRef", "PSI-MS");
        //            output.WriteAttributeString("value", psm.MetaMorpheusScore.ToString());
        //            output.WriteEndElement();  // cvParam
        //            output.WriteStartElement("cvParam");
        //            output.WriteAttributeString("accession", "MS:1002354");
        //            output.WriteAttributeString("name", "PSM-level q-value");
        //            output.WriteAttributeString("cvRef", "PSI-MS");
        //            output.WriteAttributeString("value", psm_with_fdr.QValue.ToString());
        //            output.WriteEndElement();  // cvParam
        //            output.WriteEndElement();  // SpectrumIdentificationItem
        //            output.WriteEndElement();  // SpectrumIdentificationResult
        //            index++;
        //        }
        //        output.WriteEndElement();  // SpectrumIdentificationList
        //        output.WriteEndElement();  // AnalysisData
        //        output.WriteEndElement();  // DataCollection

        //        output.WriteEndElement();  // MzIdentML

        //        output.WriteEndDocument();
        //    }
        //}

        //public static void MyAnalysis(IEnumerable<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> psms_with_fdr, string filepath)
        //{
        //    MyTreeStructure myTreeStructure = new MyTreeStructure();
        //    myTreeStructure.GenerateBins(psms_with_fdr.Where(b => (b.QValue <= 0.01 && !b.Identification.isDecoy)).Select(b => b.Identification), 0.003);
        //    myTreeStructure.AddToBins(psms_with_fdr.Where(b => b.QValue <= 0.01).Select(b => b.Identification));

        //    using (StreamWriter output = new StreamWriter(filepath))
        //    {
        //        output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tFracLocalizeableTarget");
        //        foreach (Bin bin in myTreeStructure.finalBins.OrderByDescending(b => b.Count))
        //        {
        //            output.WriteLine(bin.MassShift.ToString(CultureInfo.InvariantCulture)
        //                + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
        //                + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
        //                + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
        //                + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
        //                + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
        //                + "\t" + (bin.Count == 0 ? double.NaN : bin.CountDecoy / bin.Count).ToString(CultureInfo.InvariantCulture)
        //                + "\t" + (bin.CountTarget == 0 ? double.NaN : bin.LocalizeableTarget / bin.CountTarget).ToString(CultureInfo.InvariantCulture)
        //                + "\t" + bin.UnimodId);
        //        }
        //    }
        //}

        #region Private Methods

        private static void ConvertModificationIndex(PeptideWithSetModifications peptide, int modificationIndex, out int location, out int residueIndex)
        {
            if (modificationIndex <= 1)
            {
                location = 0;
                residueIndex = 0;
            }
            else if (modificationIndex >= peptide.Length + 2)
            {
                location = peptide.Length + 1;
                residueIndex = peptide.Length - 1;
            }
            else
            {
                location = modificationIndex - 1;
                residueIndex = modificationIndex - 2;
            }
        }

        #endregion Private Methods
    }
}