using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Omics.Modifications;
using ThermoFisher.CommonCore.Data;
using Omics;
using Transcriptomics.Digestion;
using Omics.BioPolymerGroup;
using Omics.SpectralMatch;
using System;

namespace EngineLayer
{
    public class ProteinGroup : IBioPolymerGroup
    {
        public ProteinGroup(HashSet<IBioPolymer> proteins, HashSet<IBioPolymerWithSetMods> peptides,
            HashSet<IBioPolymerWithSetMods> uniquePeptides)
        {
            Proteins = proteins;
            ListOfProteinsOrderedByAccession = Proteins.OrderBy(p => p.Accession).ToList();
            ProteinGroupName = string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.Accession));
            AllPeptides = peptides;
            UniquePeptides = uniquePeptides;
            AllPsmsBelowOnePercentFDR = new HashSet<SpectralMatch>();
            SequenceCoverageFraction = new List<double>();
            SequenceCoverageDisplayList = new List<string>();
            SequenceCoverageDisplayListWithMods = new List<string>();
            FragmentSequenceCoverageDisplayList = new List<string>();
            ProteinGroupScore = 0;
            BestPeptideScore = 0;
            QValue = 0;
            IsDecoy = false;
            IsContaminant = false;
            ModsInfo = new List<string>();

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.IsDecoy)
                {
                    IsDecoy = true;
                    break;
                }

                if (protein.IsContaminant)
                {
                    IsContaminant = true;
                    break;
                }
            }
        }

        public bool IsDecoy { get; }

        public bool IsContaminant { get; }

        public List<SpectraFileInfo> FilesForQuantification { get; set; }

        public HashSet<IBioPolymer> Proteins { get; set; }

        public string ProteinGroupName { get; private set; }

        public double ProteinGroupScore { get; set; }

        public HashSet<IBioPolymerWithSetMods> AllPeptides { get; set; }

        public HashSet<IBioPolymerWithSetMods> UniquePeptides { get; set; }
        /// <summary>
        /// Contains all PSMs associated with this protein group that pass the configured quality threshold.
        /// The specific filtering criteria depends on the <see cref="FilterType"/> and threshold passed to 
        /// <see cref="ProteinScoringAndFdrEngine"/> during protein scoring:
        /// <list type="bullet">
        ///   <item>
        ///     <term><see cref="FilterType.QValue"/></term>
        ///     <description>PSMs where QValue ≤ threshold AND QValueNotch ≤ threshold</description>
        ///   </item>
        ///   <item>
        ///     <term><see cref="FilterType.PepQValue"/></term>
        ///     <description>PSMs where PEP_QValue ≤ threshold</description>
        ///   </item>
        /// </list>
        /// The default threshold is 0.01 (1% FDR), but this can vary based on the filter configuration.
        /// This collection is populated during <see cref="ProteinScoringAndFdrEngine.ScoreProteinGroups"/> 
        /// and is used for:
        /// <list type="bullet">
        ///   <item>Calculating the <see cref="ProteinGroupScore"/> via the <see cref="Score"/> method</item>
        ///   <item>Determining <see cref="BestPeptideScore"/>, <see cref="BestPeptideQValue"/>, and <see cref="BestPeptidePEP"/></item>
        ///   <item>Computing sequence coverage in <see cref="CalculateSequenceCoverage"/></item>
        ///   <item>Reporting the number of PSMs in protein group output</item>
        /// </list>
        /// </summary>
        /// <remarks>
        /// Note: The property name "AllPsmsBelowOnePercentFDR" is a legacy name. The actual threshold 
        /// used is determined by the FilterThreshold parameter passed to ProteinScoringAndFdrEngine.
        /// </remarks>
        public HashSet<SpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        public List<double> SequenceCoverageFraction { get; private set; }

        public List<string> SequenceCoverageDisplayList { get; private set; }

        public List<string> SequenceCoverageDisplayListWithMods { get; private set; }

        public List<string> FragmentSequenceCoverageDisplayList { get; private set; }

        public double QValue { get; set; }

        public double BestPeptideQValue { get; set; }

        /// <summary>
        /// The minimum Posterior Error Probability (PEP) among all PSMs in <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// Lower values indicate higher confidence that the best peptide identification is correct.
        /// This value is populated during <see cref="ProteinScoringAndFdrEngine.DoProteinFdr"/> and is used
        /// for protein group ranking when using PEP-based filtering (<see cref="FilterType.PepQValue"/>).
        /// </summary>
        public double BestPeptidePEP { get; set; }

        public double BestPeptideScore { get; set; }

        public int CumulativeTarget { get; set; }

        public int CumulativeDecoy { get; set; }

        public bool DisplayModsOnPeptides { get; set; }

        public List<string> ModsInfo { get; private set; }

        public Dictionary<SpectraFileInfo, double> IntensitiesByFile { get; set; }

        /// <summary>
        /// Samples that carry quantification for this group, in output column order, as written back by
        /// mzLib's QuantificationEngine.
        /// </summary>
        /// <remarks>
        /// This is deliberately NOT <see cref="FilesForQuantification"/>. That pair
        /// (<see cref="FilesForQuantification"/> / <see cref="IntensitiesByFile"/>) is keyed by
        /// <see cref="SpectraFileInfo"/> and is the label-free path FlashLFQ fills; an isobaric run's
        /// samples are <see cref="IsobaricQuantSampleInfo"/> channels, several per file, which that
        /// dictionary cannot represent without losing the channel. Keeping the two side by side leaves
        /// the LFQ path untouched while the engine populates this one.
        /// </remarks>
        public List<ISampleInfo> SamplesForQuantification { get; set; }

        /// <summary>
        /// Per-sample intensities for this group, written back by mzLib's QuantificationEngine.
        /// A sample the group was not observed in is absent rather than present with a zero.
        /// </summary>
        public Dictionary<ISampleInfo, double> IntensitiesBySample { get; set; }

        /// <summary>
        /// Per-sample-group quantification and modification-occupancy results, built from
        /// <see cref="IntensitiesBySample"/>. Null until something populates it.
        /// </summary>
        public List<SampleGroupResult> SampleGroupResults { get; set; }

        /// <summary>
        /// True if any protein in this group is an entrapment sequence.
        /// </summary>
        /// <remarks>
        /// Computed on demand rather than in the constructor, unlike <see cref="IsDecoy"/> and
        /// <see cref="IsContaminant"/>, so that adding it changes no existing construction path.
        /// </remarks>
        public bool IsEntrapment => Proteins != null && Proteins.Any(p => p.IsEntrapment);

        private List<IBioPolymer> ListOfProteinsOrderedByAccession;

        private string UniquePeptidesOutput;
        private string SharedPeptidesOutput;

        //Get unique and identified peptides for output
        //Convert the output if it's a SILAC experiment
        public void GetIdentifiedPeptidesOutput(List<SilacLabel> labels)
        {
            var SharedPeptides = AllPeptides.Except(UniquePeptides);
            if (labels == null)
            {
                //TODO add unit test with displaymodsonpeptides
                if (!DisplayModsOnPeptides)
                {
                    UniquePeptidesOutput =
                        GlobalVariables.CheckLengthOfOutput(string.Join("|",
                            UniquePeptides.Select(p => p.BaseSequence).Distinct()));
                    SharedPeptidesOutput =
                        GlobalVariables.CheckLengthOfOutput(string.Join("|",
                            SharedPeptides.Select(p => p.BaseSequence).Distinct()));
                }
                else
                {
                    UniquePeptidesOutput =
                        GlobalVariables.CheckLengthOfOutput(string.Join("|",
                            UniquePeptides.Select(p => p.FullSequence).Distinct()));
                    SharedPeptidesOutput =
                        GlobalVariables.CheckLengthOfOutput(string.Join("|",
                            SharedPeptides.Select(p => p.FullSequence).Distinct()));
                }
            }
            else
            {
                if (!DisplayModsOnPeptides)
                {
                    UniquePeptidesOutput = GlobalVariables.CheckLengthOfOutput(string.Join("|",
                        UniquePeptides.Select(p =>
                            SilacConversions.GetAmbiguousLightSequence(p.BaseSequence, labels, true)).Distinct()));
                    SharedPeptidesOutput = GlobalVariables.CheckLengthOfOutput(string.Join("|",
                        SharedPeptides.Select(p =>
                            SilacConversions.GetAmbiguousLightSequence(p.BaseSequence, labels, true)).Distinct()));
                }
                else
                {
                    UniquePeptidesOutput = GlobalVariables.CheckLengthOfOutput(string.Join("|",
                        UniquePeptides.Select(p =>
                            SilacConversions.GetAmbiguousLightSequence(p.FullSequence, labels, false)).Distinct()));
                    SharedPeptidesOutput = GlobalVariables.CheckLengthOfOutput(string.Join("|",
                        SharedPeptides.Select(p =>
                            SilacConversions.GetAmbiguousLightSequence(p.FullSequence, labels, false)).Distinct()));
                }
            }
        }

        public string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("Protein Accession" + '\t');
            sb.Append("Gene" + '\t');
            sb.Append("Organism" + '\t');
            sb.Append("Protein Full Name" + '\t');
            sb.Append("Protein Unmodified Mass" + '\t');
            sb.Append("Number of Proteins in Group" + '\t');
            sb.Append("Unique Peptides" + '\t');
            sb.Append("Shared Peptides" + '\t');
            sb.Append("Number of Peptides" + '\t');
            sb.Append("Number of Unique Peptides" + '\t');
            sb.Append("Sequence Coverage Fraction" + '\t');
            sb.Append("Sequence Coverage" + '\t');
            sb.Append("Sequence Coverage with Mods" + '\t');
            sb.Append("Fragment Sequence Coverage" + '\t');
            sb.Append("Modification Info List" + "\t");
            if (FilesForQuantification != null)
            {
                bool unfractionated = FilesForQuantification.Select(p => p.Fraction).Distinct().Count() == 1;
                bool conditionsUndefined = FilesForQuantification.All(p => string.IsNullOrEmpty(p.Condition));

                // this is a hacky way to test for SILAC-labeled data...
                // Currently SILAC will report 1 column of intensities per label per spectra file, and is NOT summarized
                // into biorep-level intensity values. the SILAC code uses the "condition" field to organize this info,
                // even if the experimental design is not defined by the user. So the following bool is a way to distinguish
                // between experimental design being used in SILAC automatically vs. being defined by the user
                bool silacExperimentalDesign =
                    FilesForQuantification.Any(p => !File.Exists(p.FullFilePathWithExtension));

                foreach (var sampleGroup in FilesForQuantification.GroupBy(p => p.Condition))
                {
                    foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                    {
                        if ((conditionsUndefined && unfractionated) || silacExperimentalDesign)
                        {
                            // if the data is unfractionated and the conditions haven't been defined, just use the file name as the intensity header
                            sb.Append("Intensity_" + sample.First().FilenameWithoutExtension + "\t");
                        }
                        else
                        {
                            // if the data is fractionated and/or the conditions have been defined, label the header w/ the condition and biorep number
                            sb.Append("Intensity_" + sample.First().Condition + "_" +
                                      (sample.First().BiologicalReplicate + 1) + "\t");
                        }
                    }
                }
            }

            sb.Append("Number of PSMs" + '\t');
            sb.Append("Protein Decoy/Contaminant/Target" + '\t');
            sb.Append("Protein Cumulative Target" + '\t');
            sb.Append("Protein Cumulative Decoy" + '\t');
            sb.Append("Protein QValue" + '\t');
            sb.Append("Best Peptide Score" + '\t');
            sb.Append("Best Peptide Notch QValue" + '\t');
            sb.Append("Best Peptide PEP");
            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of protein accession numbers
            sb.Append(ProteinGroupName);
            sb.Append("\t");

            // genes
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|",
                ListOfProteinsOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault()))));
            sb.Append("\t");

            // organisms
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|",
                ListOfProteinsOrderedByAccession.Select(p => p.Organism).Distinct())));
            sb.Append("\t");

            // list of protein names
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|",
                ListOfProteinsOrderedByAccession.Select(p => p.FullName).Distinct())));
            sb.Append("\t");

            // list of masses
            var sequences = ListOfProteinsOrderedByAccession.Select(p => p.BaseSequence).Distinct();
            List<double> masses = new List<double>();
            foreach (var sequence in sequences)
            {
                try
                {
                    if (GlobalVariables.AnalyteType == AnalyteType.Oligo)
                        masses.Add(new OligoWithSetMods(sequence, GlobalVariables.AllRnaModsKnownDictionary).MonoisotopicMass);
                    else
                        masses.Add(new Proteomics.AminoAcidPolymer.Peptide(sequence).MonoisotopicMass);

                }
                catch (System.Exception)
                {
                    masses.Add(double.NaN);
                }
            }

            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", masses)));
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + Proteins.Count);
            sb.Append("\t");

            // list of unique peptides
            if (UniquePeptidesOutput != null)
            {
                sb.Append(GlobalVariables.CheckLengthOfOutput(UniquePeptidesOutput));
            }

            sb.Append("\t");

            // list of shared peptides
            if (SharedPeptidesOutput != null)
            {
                sb.Append(GlobalVariables.CheckLengthOfOutput(SharedPeptidesOutput));
            }

            sb.Append("\t");

            // number of peptides
            if (!DisplayModsOnPeptides)
            {
                sb.Append("" + AllPeptides.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + AllPeptides.Select(p => p.FullSequence).Distinct().Count());
            }

            sb.Append("\t");

            // number of unique peptides
            if (!DisplayModsOnPeptides)
            {
                sb.Append("" + UniquePeptides.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + UniquePeptides.Select(p => p.FullSequence).Distinct().Count());
            }

            sb.Append("\t");

            // sequence coverage percent
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|",
                SequenceCoverageFraction.Select(p => string.Format("{0:0.#####}", p)))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", SequenceCoverageDisplayList)));
            sb.Append("\t");

            // sequence coverage with mods
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", SequenceCoverageDisplayListWithMods)));
            sb.Append("\t");

            // fragment sequence coverage
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", FragmentSequenceCoverageDisplayList)));
            sb.Append("\t");

            //Detailed mods information list
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", ModsInfo)));
            sb.Append("\t");

            // MS1 intensity (retrieved from FlashLFQ in the SearchTask)
            if (IntensitiesByFile != null && FilesForQuantification != null)
            {
                foreach (var sampleGroup in FilesForQuantification.GroupBy(p => p.Condition))
                {
                    foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                    {
                        // if the samples are fractionated, the protein will only have 1 intensity in the first fraction
                        // and the other fractions will be zero. we could find the first/only fraction with an intensity,
                        // but simply summing the fractions is easier than finding the single non-zero value
                        double summedIntensity = sample.Sum(file => IntensitiesByFile[file]);

                        if (summedIntensity > 0)
                        {
                            sb.Append(summedIntensity);
                        }

                        sb.Append("\t");
                    }
                }
            }

            // number of PSMs for listed peptides
            sb.Append("" + AllPsmsBelowOnePercentFDR.Count);
            sb.Append("\t");

            // isDecoy
            if (IsDecoy)
            {
                sb.Append("D");
            }
            else if (IsContaminant)
            {
                sb.Append("C");
            }
            else
            {
                sb.Append("T");
            }

            sb.Append("\t");

            // cumulative target
            sb.Append(CumulativeTarget);
            sb.Append("\t");

            // cumulative decoy
            sb.Append(CumulativeDecoy);
            sb.Append("\t");

            // q value
            sb.Append(QValue);
            sb.Append("\t");

            // best peptide score
            sb.Append(BestPeptideScore);
            sb.Append("\t");

            // best peptide q value
            sb.Append(BestPeptideQValue);
            sb.Append("\t");

            // best peptide PEP
            sb.Append(BestPeptidePEP);

            return sb.ToString();
        }

        // Score() method is only used internally, to make protein grouping faster
        // this is NOT an output and is NOT used for protein FDR calculations
        public void Score()
        {
            // sum the scores of the best PSM per base sequence
            ProteinGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.BaseSequence)
                .Select(p => p.Select(x => x.Score).Max()).Sum();
        }

        public void CalculateSequenceCoverage()
        {
            var proteinsWithUnambigSeqPsms = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();
            var proteinsWithPsmsWithLocalizedMods = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            foreach (var protein in Proteins)
            {
                proteinsWithUnambigSeqPsms.Add(protein, new List<IBioPolymerWithSetMods>());
                proteinsWithPsmsWithLocalizedMods.Add(protein, new List<IBioPolymerWithSetMods>());
            }

            foreach (var psm in AllPsmsBelowOnePercentFDR)
            {
                // null BaseSequence means that the amino acid sequence is ambiguous; do not use these to calculate sequence coverage
                if (psm.BaseSequence != null)
                {
                    psm.GetAminoAcidCoverage();

                    foreach (var peptide in psm.BestMatchingBioPolymersWithSetMods.Select(psm => psm.SpecificBioPolymer).DistinctBy(pep => pep.FullSequence))
                    {
                        // might be unambiguous but also shared; make sure this protein group contains this peptide+protein combo
                        if (Proteins.Contains(peptide.Parent))
                        {
                            proteinsWithUnambigSeqPsms[peptide.Parent].Add(peptide);

                            // null FullSequence means that mods were not successfully localized; do not display them on the sequence coverage mods info
                            if (peptide.FullSequence != null)
                            {
                                proteinsWithPsmsWithLocalizedMods[peptide.Parent].Add(peptide);
                            }
                        }
                    }
                    
                }
            }

            //Calculate sequence coverage at the amino acid level by looking at fragment specific coverage
            //loop through proteins
            foreach (IBioPolymer protein in ListOfProteinsOrderedByAccession)
            {
                //create a hash set for storing covered one-based residue numbers of protein
                HashSet<int> coveredResiduesInProteinOneBased = new();

                //loop through PSMs
                foreach (SpectralMatch psm in AllPsmsBelowOnePercentFDR.Where(psm => psm.BaseSequence != null))
                {
                    //Calculate the covered bases within the psm. This is one based numbering for the peptide only
                    psm.GetAminoAcidCoverage();
                    if (psm.FragmentCoveragePositionInPeptide == null) continue;
                    //loop through each peptide within the psm
                    IEnumerable<IBioPolymerWithSetMods> pwsms = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer)
                        .Where(p => p.Parent.Accession == protein.Accession);
                    foreach (var pwsm in pwsms)
                    {
                        //create a hashset to store the covered residues for the peptide, converted to the corresponding indices of the protein
                        HashSet<int> coveredResiduesInPeptide = new();
                        //add the peptide start position within the protein to each covered index of the psm
                        foreach (var position in psm.FragmentCoveragePositionInPeptide)
                        {
                            coveredResiduesInPeptide.Add(position + pwsm.OneBasedStartResidue -
                                                         1); //subtract one because these are both one based
                        }

                        //Add the peptide specific positions, to the overall hashset for the protein
                        coveredResiduesInProteinOneBased.UnionWith(coveredResiduesInPeptide);
                    }
                }

                // create upper/lowercase string
                char[] fragmentCoverageArray = protein.BaseSequence.ToLower().ToCharArray();
                foreach (var residue in coveredResiduesInProteinOneBased)
                {
                    fragmentCoverageArray[residue - 1] = char.ToUpper(fragmentCoverageArray[residue - 1]);
                }

                FragmentSequenceCoverageDisplayList.Add(new string(fragmentCoverageArray));
            }

            //Calculates the coverage at the peptide level... if a peptide is present all of the AAs in the peptide are covered
            foreach (var protein in ListOfProteinsOrderedByAccession)
            {
                HashSet<int> coveredOneBasedResidues = new HashSet<int>();

                // get residue numbers of each peptide in the protein and identify them as observed if the sequence is unambiguous
                foreach (var peptide in proteinsWithUnambigSeqPsms[protein])
                {
                    for (int i = peptide.OneBasedStartResidue; i <= peptide.OneBasedEndResidue; i++)
                    {
                        coveredOneBasedResidues.Add(i);
                    }
                }

                // calculate sequence coverage percent
                double seqCoverageFract = (double)coveredOneBasedResidues.Count / protein.Length;

                // add the percent coverage
                SequenceCoverageFraction.Add(seqCoverageFract);

                // convert the observed amino acids to upper case if they are unambiguously observed
                string sequenceCoverageDisplay = protein.BaseSequence.ToLower();
                var coverageArray = sequenceCoverageDisplay.ToCharArray();
                foreach (var obsResidueLocation in coveredOneBasedResidues)
                {
                    coverageArray[obsResidueLocation - 1] = char.ToUpper(coverageArray[obsResidueLocation - 1]);
                }

                sequenceCoverageDisplay = new string(coverageArray);

                // add the coverage display
                SequenceCoverageDisplayList.Add(sequenceCoverageDisplay);

                // put mods in the sequence coverage display
                // get mods to display in sequence (only unambiguously identified mods)
                var modsOnThisProtein = new HashSet<KeyValuePair<int, Modification>>();
                foreach (var pep in proteinsWithPsmsWithLocalizedMods[protein])
                {
                    foreach (var mod in pep.AllModsOneIsNterminus)
                    {
                        if (!mod.Value.ModificationType.Contains("PeptideTermMod")
                            && !mod.Value.ModificationType.Contains("Common Variable")
                            && !mod.Value.ModificationType.Contains("Common Fixed"))
                        {
                            modsOnThisProtein.Add(
                                new KeyValuePair<int, Modification>(pep.OneBasedStartResidue + mod.Key - 2,
                                    mod.Value));
                        }
                    }
                }

                var tempMods = modsOnThisProtein.OrderBy(p => p.Key).ToList();
                foreach (var mod in tempMods)
                {
                    if (mod.Value.LocationRestriction.Equals("N-terminal."))
                    {
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(
                            0,
                            $"[{mod.Value.IdWithMotif}]-");
                    }
                    else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                    {
                        int modStringIndex = sequenceCoverageDisplay.Length - (protein.Length - mod.Key);
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(
                            modStringIndex,
                            $"[{mod.Value.IdWithMotif}]");
                    }
                    else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                    {
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(
                            sequenceCoverageDisplay.Length,
                            $"-[{mod.Value.IdWithMotif}]");
                    }
                }

                SequenceCoverageDisplayListWithMods.Add(sequenceCoverageDisplay);

                if (!modsOnThisProtein.Any())
                {
                    continue;
                }

                // calculate spectral count % of modified observations
                var pepModTotals = new List<int>(); // count of modified peptides for each mod/index
                var pepTotals = new List<int>(); // count of all peptides for each mod/index
                var modIndex = new List<(int index, string modName)>(); // index and name of the modified position

                foreach (var pep in proteinsWithPsmsWithLocalizedMods[protein])
                {
                    foreach (var mod in pep.AllModsOneIsNterminus)
                    {
                        int pepNumTotal = 0; //For one mod, The total Pep Num

                        if (mod.Value.ModificationType.Contains("Common Variable")
                            || mod.Value.ModificationType.Contains("Common Fixed")
                            || mod.Value.LocationRestriction.Equals(ModLocationOnPeptideOrProtein.PepC)
                            || mod.Value.LocationRestriction.Equals(ModLocationOnPeptideOrProtein.NPep))
                        {
                            continue;
                        }

                        int indexInProtein;
                        if (mod.Value.LocationRestriction.Equals("N-terminal."))
                        {
                            indexInProtein = 1;
                        }
                        else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                        {
                            indexInProtein = pep.OneBasedStartResidue + mod.Key - 2;
                        }
                        else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                        {
                            indexInProtein = protein.Length;
                        }
                        else
                        {
                            // In case it's a peptide terminal mod, skip!
                            // we don't want this annotated in the protein's modifications
                            continue;
                        }

                        var modKey = (indexInProtein, mod.Value.IdWithMotif);
                        if (modIndex.Contains(modKey))
                        {
                            pepModTotals[modIndex.IndexOf(modKey)] += 1;
                        }
                        else
                        {
                            modIndex.Add(modKey);
                            foreach (var pept in proteinsWithPsmsWithLocalizedMods[protein])
                            {
                                if (indexInProtein >= pept.OneBasedStartResidue - (indexInProtein == 1 ? 1 : 0)
                                    && indexInProtein <= pept.OneBasedEndResidue)
                                {
                                    pepNumTotal += 1;
                                }
                            }

                            pepTotals.Add(pepNumTotal);
                            pepModTotals.Add(1);
                        }
                    }
                }

                var modStrings = new List<(int aaNum, string part)>();
                for (int i = 0; i < pepModTotals.Count; i++)
                {
                    string aa = modIndex[i].index.ToString();
                    string modName = modIndex[i].modName.ToString();
                    string occupancy = ((double)pepModTotals[i] / (double)pepTotals[i]).ToString("F2");
                    string fractOccupancy = $"{pepModTotals[i].ToString()}/{pepTotals[i].ToString()}";
                    string tempString = ($"#aa{aa}[{modName},info:occupancy={occupancy}({fractOccupancy})]");
                    modStrings.Add((modIndex[i].index, tempString));
                }

                var modInfoString = string.Join(";", modStrings.OrderBy(x => x.aaNum).Select(x => x.part));

                if (!string.IsNullOrEmpty(modInfoString))
                {
                    ModsInfo.Add(modInfoString);
                }
            }
        }

        public void MergeProteinGroupWith(ProteinGroup other)
        {
            this.Proteins.UnionWith(other.Proteins);
            this.AllPeptides.UnionWith(other.AllPeptides);
            this.UniquePeptides.UnionWith(other.UniquePeptides);
            this.AllPsmsBelowOnePercentFDR.UnionWith(other.AllPsmsBelowOnePercentFDR);
            other.ProteinGroupScore = 0;

            ListOfProteinsOrderedByAccession = Proteins.OrderBy(p => p.Accession).ToList();

            ProteinGroupName = string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.Accession));
        }

        public ProteinGroup ConstructSubsetProteinGroup(string fullFilePath, List<SilacLabel> silacLabels = null)
        {
            var allPsmsForThisFile =
                new HashSet<SpectralMatch>(
                    AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allPeptidesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(
                    allPsmsForThisFile.SelectMany(p => p.BestMatchingBioPolymersWithSetMods.Select(v => v.SpecificBioPolymer)));
            var allUniquePeptidesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(UniquePeptides.Intersect(allPeptidesForThisFile));

            ProteinGroup subsetPg = new ProteinGroup(Proteins, allPeptidesForThisFile, allUniquePeptidesForThisFile)
            {
                AllPsmsBelowOnePercentFDR = allPsmsForThisFile,
                DisplayModsOnPeptides = DisplayModsOnPeptides
            };

            SpectraFileInfo spectraFileInfo = null;
            if (FilesForQuantification != null)
            {
                spectraFileInfo = FilesForQuantification.Where(p => p.FullFilePathWithExtension == fullFilePath)
                    .FirstOrDefault();
                //check that file name wasn't changed (can occur in SILAC searches)
                if (!silacLabels.IsNullOrEmpty() && spectraFileInfo == null)
                {
                    foreach (SilacLabel label in silacLabels)
                    {
                        string fakeFilePath = SilacConversions
                            .GetHeavyFileInfo(new SpectraFileInfo(fullFilePath, "", 0, 0, 0), label)
                            .FullFilePathWithExtension;
                        spectraFileInfo = FilesForQuantification.Where(p => p.FullFilePathWithExtension == fakeFilePath)
                            .FirstOrDefault();
                        if (spectraFileInfo != null)
                        {
                            break;
                        }
                    }

                    //if still no hits, might be SILAC turnover
                    if (spectraFileInfo == null)
                    {
                        string filepathWithoutExtension = Path.Combine(Path.GetDirectoryName(fullFilePath),
                            Path.GetFileNameWithoutExtension(fullFilePath));
                        string extension = Path.GetExtension(fullFilePath);
                        string fakeFilePath = filepathWithoutExtension + SilacConversions.ORIGINAL_TURNOVER_LABEL_NAME +
                                              extension;
                        spectraFileInfo = FilesForQuantification.Where(p => p.FullFilePathWithExtension == fakeFilePath)
                            .FirstOrDefault();
                    }
                }

                subsetPg.FilesForQuantification = new List<SpectraFileInfo> { spectraFileInfo };
            }

            if (IntensitiesByFile == null)
            {
                subsetPg.IntensitiesByFile = null;
            }
            else
            {
                subsetPg.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>
                    { { spectraFileInfo, IntensitiesByFile[spectraFileInfo] } };
            }

            return subsetPg;
        }

        //method only considers accessions, not peptides
        public bool Equals(ProteinGroup grp)
        {
            //Check for null and compare run-time types.
            if (grp == null) 
            {
                return false;
            }
            else if (!this.ListOfProteinsOrderedByAccession.Select(a=>a.Accession).ToList().SequenceEqual(grp.ListOfProteinsOrderedByAccession.Select(a => a.Accession).ToList()))
            {
                return false;
            }

            return true;
        }

        #region IBioPolymerGroup

        // mzLib's QuantificationEngine takes List<IBioPolymerGroup> in and writes its results back onto
        // the same objects. Implementing the interface here -- rather than wrapping ProteinGroup in an
        // adapter -- is what lets those results land on the real protein groups that MetaMorpheus's
        // writers already hold.
        //
        // Every member below is an EXPLICIT implementation, because each one collides with a
        // MetaMorpheus name that means the same thing ("Protein" where the interface says "BioPolymer").
        // Nothing here introduces new state or changes existing behaviour; the members that already
        // match by name and signature -- IsDecoy, IsContaminant, QValue, CumulativeTarget,
        // CumulativeDecoy, Score(), CalculateSequenceCoverage(), GetTabSeparatedHeader(), and the three
        // quantification properties above -- satisfy the interface as they stand.

        HashSet<IBioPolymer> IBioPolymerGroup.BioPolymers { get => Proteins; set => Proteins = value; }

        string IBioPolymerGroup.BioPolymerGroupName => ProteinGroupName;

        double IBioPolymerGroup.BioPolymerGroupScore { get => ProteinGroupScore; set => ProteinGroupScore = value; }

        HashSet<IBioPolymerWithSetMods> IBioPolymerGroup.AllBioPolymersWithSetMods { get => AllPeptides; set => AllPeptides = value; }

        HashSet<IBioPolymerWithSetMods> IBioPolymerGroup.UniqueBioPolymersWithSetMods { get => UniquePeptides; set => UniquePeptides = value; }

        double IBioPolymerGroup.BestBioPolymerWithSetModsQValue { get => BestPeptideQValue; set => BestPeptideQValue = value; }

        double IBioPolymerGroup.BestBioPolymerWithSetModsScore { get => BestPeptideScore; set => BestPeptideScore = value; }

        List<IBioPolymer> IBioPolymerGroup.ListOfBioPolymersOrderedByAccession => ListOfProteinsOrderedByAccession;

        /// <remarks>
        /// MetaMorpheus's ProteinGroup always represents parent sequences (proteins or nucleic acids),
        /// never digestion products, so occupancy is calculated at parent-level positions.
        /// </remarks>
        BioPolymerGroupType IBioPolymerGroup.GroupType => BioPolymerGroupType.Parent;

        /// <remarks>
        /// A PROJECTION, not the live set. <see cref="AllPsmsBelowOnePercentFDR"/> is a
        /// <c>HashSet&lt;SpectralMatch&gt;</c> and the interface asks for <c>HashSet&lt;ISpectralMatch&gt;</c>,
        /// which are unrelated types however compatible their elements are. The getter therefore copies,
        /// so mutating what it returns does nothing to the group -- assign a whole set instead. The
        /// quantification engine only reads peptide membership, never this, so the copy is not on any
        /// hot path.
        /// </remarks>
        HashSet<ISpectralMatch> IBioPolymerGroup.AllPsmsBelowOnePercentFDR
        {
            get => new HashSet<ISpectralMatch>(AllPsmsBelowOnePercentFDR);
            set => AllPsmsBelowOnePercentFDR = value == null
                ? new HashSet<SpectralMatch>()
                : new HashSet<SpectralMatch>(value.OfType<SpectralMatch>());
        }

        void IBioPolymerGroup.MergeWith(IBioPolymerGroup otherBioPolymerGroup)
        {
            if (otherBioPolymerGroup is not ProteinGroup other)
            {
                throw new ArgumentException(
                    "A ProteinGroup can only be merged with another ProteinGroup, not with " +
                    $"{otherBioPolymerGroup?.GetType().Name ?? "null"}.", nameof(otherBioPolymerGroup));
            }

            MergeProteinGroupWith(other);
        }

        IBioPolymerGroup IBioPolymerGroup.ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel> silacLabels)
            => ConstructSubsetProteinGroup(fullFilePath, silacLabels);

        /// <summary>
        /// Reference equality, deliberately, and NOT the accession comparison of
        /// <see cref="Equals(ProteinGroup)"/>.
        /// </summary>
        /// <remarks>
        /// ProteinGroup does not override <see cref="object.GetHashCode"/>, so its hash is the reference
        /// hash, and mzLib's QuantMatrix keys its rows in a
        /// <c>Dictionary&lt;IBioPolymerGroup, int&gt;</c> with the default comparer -- which pairs this
        /// method with that hash. An accession-based Equals here would give keys that compare equal but
        /// hash differently, and rows would go missing.
        ///
        /// Overriding GetHashCode to fix that is not an option either: the hash would have to come from
        /// <see cref="ProteinGroupName"/>, which <see cref="MergeProteinGroupWith"/> reassigns, and a
        /// mutable hash key loses entries in any dictionary that already holds the group. Within one
        /// search the distinction is moot -- protein groups are unique by accession and each exists
        /// once -- so the two definitions agree on everything the engine actually does.
        /// </remarks>
        bool IEquatable<IBioPolymerGroup>.Equals(IBioPolymerGroup other) => ReferenceEquals(this, other);

        #endregion
    }
}