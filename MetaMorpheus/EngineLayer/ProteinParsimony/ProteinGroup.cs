using Proteomics;
using FlashLFQ;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Omics.Modifications;
using Omics;
using Omics.BioPolymerGroup;
using Omics.SpectralMatch;
using Transcriptomics.Digestion;
using MzLibUtil;
using Easy.Common.Extensions;

namespace EngineLayer
{
    /// <summary>
    /// MetaMorpheus-specific protein group, extending the generic BioPolymerGroup from mzLib.
    /// Adds PEP-based scoring, SILAC-aware peptide output, and MetaMorpheus-specific TSV column names.
    /// Quantification and modification occupancy are handled by the base class's SampleGroupResult system.
    ///
    /// Backward-compatible alias properties (e.g., Proteins, AllPeptides, ProteinGroupScore) delegate
    /// to the corresponding BioPolymerGroup base properties. Consumers should gradually migrate to
    /// the base class names.
    ///
    /// Score() is handled entirely by the base class.
    /// CalculateSequenceCoverage() and GetTabSeparatedHeader() are hidden (via new) because
    /// CalculateSequenceCoverage() accesses MetaMorpheus-specific SpectralMatch members
    /// (GetAminoAcidCoverage, BestMatchingBioPolymersWithSetMods, FragmentCoveragePositionInPeptide)
    /// that are not on the ISpectralMatch interface, and GetTabSeparatedHeader() uses
    /// MetaMorpheus-specific column names and includes BestPeptidePEP.
    /// </summary>
    public class ProteinGroup : BioPolymerGroup
    {
        public ProteinGroup(HashSet<IBioPolymer> proteins, HashSet<IBioPolymerWithSetMods> peptides,
            HashSet<IBioPolymerWithSetMods> uniquePeptides)
            : base(proteins, peptides, uniquePeptides)
        {
            AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            SequenceCoverageFraction = new List<double>();
            SequenceCoverageDisplayList = new List<string>();
            SequenceCoverageDisplayListWithMods = new List<string>();
            FragmentSequenceCoverageDisplayList = new List<string>();
            BestPeptideScore = 0;
        }

        #region Backward-Compatible Property Aliases

        /// <summary>Maps to <see cref="BioPolymerGroup.BioPolymers"/>.</summary>
        public HashSet<IBioPolymer> Proteins
        {
            get => BioPolymers;
            set => BioPolymers = value;
        }

        /// <summary>Maps to <see cref="BioPolymerGroup.BioPolymerGroupName"/>.</summary>
        public string ProteinGroupName => BioPolymerGroupName;

        /// <summary>Maps to <see cref="BioPolymerGroup.BioPolymerGroupScore"/>.</summary>
        public double ProteinGroupScore
        {
            get => BioPolymerGroupScore;
            set => BioPolymerGroupScore = value;
        }

        /// <summary>Maps to <see cref="BioPolymerGroup.AllBioPolymersWithSetMods"/>.</summary>
        public HashSet<IBioPolymerWithSetMods> AllPeptides
        {
            get => AllBioPolymersWithSetMods;
            set => AllBioPolymersWithSetMods = value;
        }

        /// <summary>Maps to <see cref="BioPolymerGroup.UniqueBioPolymersWithSetMods"/>.</summary>
        public HashSet<IBioPolymerWithSetMods> UniquePeptides
        {
            get => UniqueBioPolymersWithSetMods;
            set => UniqueBioPolymersWithSetMods = value;
        }

        /// <summary>Maps to <see cref="BioPolymerGroup.BestBioPolymerWithSetModsScore"/>.</summary>
        public double BestPeptideScore
        {
            get => BestBioPolymerWithSetModsScore;
            set => BestBioPolymerWithSetModsScore = value;
        }

        /// <summary>Maps to <see cref="BioPolymerGroup.BestBioPolymerWithSetModsQValue"/>.</summary>
        public double BestPeptideQValue
        {
            get => BestBioPolymerWithSetModsQValue;
            set => BestBioPolymerWithSetModsQValue = value;
        }

        /// <summary>Maps to <see cref="BioPolymerGroup.ListOfBioPolymersOrderedByAccession"/>.</summary>
        public List<IBioPolymer> ListOfProteinsOrderedByAccession => ListOfBioPolymersOrderedByAccession;

        /// <summary>Maps to <see cref="BioPolymerGroup.SamplesForQuantification"/>. Filtered to SpectraFileInfo.</summary>
        public List<SpectraFileInfo> FilesForQuantification
        {
            get => SamplesForQuantification?.OfType<SpectraFileInfo>().ToList();
            set => SamplesForQuantification = value?.Cast<ISampleInfo>().ToList();
        }

        /// <summary>Maps to <see cref="BioPolymerGroup.IntensitiesBySample"/>. Keyed by SpectraFileInfo.</summary>
        public Dictionary<SpectraFileInfo, double> IntensitiesByFile
        {
            get => IntensitiesBySample?.ToDictionary(kvp => (SpectraFileInfo)kvp.Key, kvp => kvp.Value);
            set => IntensitiesBySample = value?.ToDictionary(kvp => (ISampleInfo)kvp.Key, kvp => kvp.Value);
        }

        #endregion

        #region MetaMorpheus-Specific Properties

        /// <summary>
        /// The minimum Posterior Error Probability (PEP) among all PSMs in <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// Lower values indicate higher confidence. Populated during protein FDR and used for PEP-based ranking.
        /// </summary>
        public double BestPeptidePEP { get; set; }

        // Sequence coverage stored as flat lists (MM-specific format).
        // BioPolymerGroup uses CoverageResult instead; these are kept for TSV output compatibility.
        public List<double> SequenceCoverageFraction { get; private set; }
        public List<string> SequenceCoverageDisplayList { get; private set; }
        public List<string> SequenceCoverageDisplayListWithMods { get; private set; }
        public List<string> FragmentSequenceCoverageDisplayList { get; private set; }

        private string UniquePeptidesOutput;
        private string SharedPeptidesOutput;

        #endregion

        #region Peptide Output

        /// <summary>
        /// Populates unique and shared peptide output strings, converting to light SILAC sequences if needed.
        /// </summary>
        public void GetIdentifiedPeptidesOutput(List<SilacLabel> labels)
        {
            var SharedPeptides = AllPeptides.Except(UniquePeptides);
            if (labels == null)
            {
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

        #endregion

        #region TSV Output (MetaMorpheus column names, base quantification format)

        /// <summary>
        /// MetaMorpheus TSV header with "Protein" column names and BestPeptidePEP.
        /// Quantification/occupancy columns use the base BioPolymerGroup SampleGroupResult format.
        /// </summary>
        public new string GetTabSeparatedHeader()
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

            // Quantification and occupancy columns from base SampleGroupResult system.
            // Dynamic columns appear only when SampleGroupResults has been explicitly populated
            // upstream (e.g., LFQ-success path). Workflows that return early from quantification
            // leave it null/empty and emit only the static columns.
            if (SampleGroupResults != null)
            {
                foreach (var group in SampleGroupResults)
                {
                    sb.Append($"SpectralCount_{group.Label}\t");
                    if (group.HasIntensityData)
                        sb.Append($"Intensity_{group.Label}\t");
                    sb.Append($"CountOccupancy_{group.Label}\t");
                    if (group.HasIntensityData)
                        sb.Append($"IntensityOccupancy_{group.Label}\t");
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
                sb.Append("" + AllPeptides.Select(p => p.BaseSequence).Distinct().Count());
            else
                sb.Append("" + AllPeptides.Select(p => p.FullSequence).Distinct().Count());
            sb.Append("\t");

            // number of unique peptides
            if (!DisplayModsOnPeptides)
                sb.Append("" + UniquePeptides.Select(p => p.BaseSequence).Distinct().Count());
            else
                sb.Append("" + UniquePeptides.Select(p => p.FullSequence).Distinct().Count());
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

            // Quantification and occupancy from base SampleGroupResult system.
            // Mirrors the header: dynamic columns appear only when SampleGroupResults has been
            // populated upstream so header and rows stay consistent within a file.
            if (SampleGroupResults != null)
            {
                bool isProteinLevel = GroupType == BioPolymerGroupType.Parent;
                IEnumerable<string> orderedKeys = isProteinLevel
                    ? ListOfProteinsOrderedByAccession.Select(p => p.Accession)
                    : AllPeptides.Select(p => p.BaseSequence).Distinct().OrderBy(s => s);

                foreach (var group in SampleGroupResults)
                {
                    sb.Append(group.SpectralCount);
                    sb.Append("\t");

                    if (group.HasIntensityData)
                    {
                        if (group.Intensity > 0)
                            sb.Append(group.Intensity);
                        sb.Append("\t");
                    }

                    sb.Append(GlobalVariables.CheckLengthOfOutput(group.FormatOccupancy(orderedKeys, isProteinLevel, intensityBased: false)));
                    sb.Append("\t");

                    if (group.HasIntensityData)
                    {
                        sb.Append(GlobalVariables.CheckLengthOfOutput(group.FormatOccupancy(orderedKeys, isProteinLevel, intensityBased: true)));
                        sb.Append("\t");
                    }
                }
            }

            // number of PSMs for listed peptides
            sb.Append("" + AllPsmsBelowOnePercentFDR.Count);
            sb.Append("\t");

            // isDecoy
            if (IsDecoy)
                sb.Append("D");
            else if (IsContaminant)
                sb.Append("C");
            else
                sb.Append("T");
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

        #endregion

        #region Coverage

        /// <summary>
        /// Computes sequence coverage using MetaMorpheus-specific SpectralMatch members
        /// (GetAminoAcidCoverage, BestMatchingBioPolymersWithSetMods, FragmentCoveragePositionInPeptide).
        /// Results are stored in the flat list properties (SequenceCoverageFraction, etc.)
        /// rather than in <see cref="BioPolymerGroup.CoverageResult"/>.
        /// </summary>
        public new void CalculateSequenceCoverage()
        {
            var proteinsWithUnambigSeqPsms = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();
            var proteinsWithPsmsWithLocalizedMods = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            foreach (var protein in Proteins)
            {
                proteinsWithUnambigSeqPsms.Add(protein, new List<IBioPolymerWithSetMods>());
                proteinsWithPsmsWithLocalizedMods.Add(protein, new List<IBioPolymerWithSetMods>());
            }

            foreach (var psm in AllPsmsBelowOnePercentFDR.OfType<SpectralMatch>())
            {
                // null BaseSequence means that the amino acid sequence is ambiguous; do not use these to calculate sequence coverage
                if (psm.BaseSequence != null)
                {
                    psm.GetAminoAcidCoverage();

                    foreach (var peptide in psm.BestMatchingBioPolymersWithSetMods
                        .Select(p => p.SpecificBioPolymer).DistinctBy(pep => pep.FullSequence))
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
                foreach (SpectralMatch psm in AllPsmsBelowOnePercentFDR.OfType<SpectralMatch>()
                    .Where(psm => psm.BaseSequence != null))
                {
                    //Calculate the covered bases within the psm. This is one based numbering for the peptide only
                    psm.GetAminoAcidCoverage();
                    if (psm.FragmentCoveragePositionInPeptide == null) continue;

                    //loop through each peptide within the psm
                    IEnumerable<IBioPolymerWithSetMods> pwsms = psm.BestMatchingBioPolymersWithSetMods
                        .Select(p => p.SpecificBioPolymer)
                        .Where(p => p.Parent.Accession == protein.Accession);

                    foreach (var pwsm in pwsms)
                    {
                        //create a hashset to store the covered residues for the peptide, converted to the corresponding indices of the protein
                        HashSet<int> coveredResiduesInPeptide = new();
                        //add the peptide start position within the protein to each covered index of the psm
                        foreach (var position in psm.FragmentCoveragePositionInPeptide)
                        {
                            coveredResiduesInPeptide.Add(position + pwsm.OneBasedStartResidue - 1); //subtract one because these are both one based
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
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(0, $"[{mod.Value.IdWithMotif}]-");
                    }
                    else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                    {
                        int modStringIndex = sequenceCoverageDisplay.Length - (protein.Length - mod.Key);
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(modStringIndex, $"[{mod.Value.IdWithMotif}]");
                    }
                    else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                    {
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(sequenceCoverageDisplay.Length, $"-[{mod.Value.IdWithMotif}]");
                    }
                }

                SequenceCoverageDisplayListWithMods.Add(sequenceCoverageDisplay);
            }
        }

        #endregion

        #region Merge and Subset

        /// <summary>
        /// Merges another ProteinGroup into this one. Delegates entirely to
        /// <see cref="BioPolymerGroup.MergeWith"/> which handles PSMs, biopolymers, peptides, and name.
        /// </summary>
        public void MergeProteinGroupWith(ProteinGroup other)
        {
            base.MergeWith(other);
        }

        /// <summary>
        /// Creates a ProteinGroup subset containing only data from the specified spectra file.
        /// </summary>
        public ProteinGroup ConstructSubsetProteinGroup(string fullFilePath, List<SilacLabel> silacLabels = null)
        {
            var allPsmsForThisFile =
                new HashSet<ISpectralMatch>(
                    AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allPeptidesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(
                    allPsmsForThisFile.SelectMany(p => p.GetIdentifiedBioPolymersWithSetMods()));
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
                if (!MzLibUtil.ClassExtensions.IsNullOrEmpty(silacLabels) && spectraFileInfo == null)
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

                if (spectraFileInfo != null)
                {
                    subsetPg.FilesForQuantification = new List<SpectraFileInfo> { spectraFileInfo };
                }
            }

            if (IntensitiesByFile == null || spectraFileInfo == null)
            {
                subsetPg.IntensitiesByFile = null;
            }
            else
            {
                subsetPg.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>
                    { { spectraFileInfo, IntensitiesByFile.GetValueOrDefault(spectraFileInfo, 0) } };
            }

            return subsetPg;
        }

        #endregion

        #region Equality

        /// <summary>
        /// Compares by ordered accession list.
        /// </summary>
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

        #endregion
    }
}
