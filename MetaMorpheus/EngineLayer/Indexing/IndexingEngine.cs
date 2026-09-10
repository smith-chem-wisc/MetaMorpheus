using Chemistry;
using EngineLayer.NonSpecificEnzymeSearch;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using Omics;
using Omics.Fragmentation.Peptide;
using Omics.Modifications;
using UsefulProteomicsDatabases;

namespace EngineLayer.Indexing
{
    public class IndexingEngine : MetaMorpheusEngine
    {
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private const int FragmentBinsPerDalton = 1000;
        private readonly List<Protein> ProteinList;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly List<SilacLabel> SilacLabels;
        private readonly (SilacLabel StartLabel, SilacLabel EndLabel)? TurnoverLabels;
        private readonly int CurrentPartition;
        private readonly DecoyType DecoyType;
        private readonly double MaxFragmentSize;
        public readonly bool GeneratePrecursorIndex;
        public readonly List<FileInfo> ProteinDatabases;
        public readonly TargetContaminantAmbiguity TcAmbiguity;

        public IndexingEngine(List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications,
            List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel, int currentPartition, DecoyType decoyType,
            CommonParameters commonParams, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, 
            double maxFragmentSize, bool generatePrecursorIndex, List<FileInfo> proteinDatabases, TargetContaminantAmbiguity tcAmbiguity, List<string> nestedIds)
            : base(commonParams, fileSpecificParameters, nestedIds)
        {
            ProteinList = proteinList;
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            SilacLabels = silacLabels;
            if (startLabel != null || endLabel != null) //else it's null
            {
                TurnoverLabels = (startLabel, endLabel);
            }

            CurrentPartition = currentPartition + 1;
            DecoyType = decoyType;
            MaxFragmentSize = maxFragmentSize;
            GeneratePrecursorIndex = generatePrecursorIndex;
            ProteinDatabases = proteinDatabases;
            TcAmbiguity = tcAmbiguity;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Databases: " + string.Join(",", ProteinDatabases.OrderBy(p => p.Name).Select(p => p.Name + ":" + p.CreationTime)));
            sb.AppendLine("Partitions: " + CurrentPartition + "/" + CommonParameters.TotalPartitions);
            sb.AppendLine("Precursor Index: " + GeneratePrecursorIndex);
            sb.AppendLine("Search Decoys: " + DecoyType);
            sb.AppendLine("Number of proteins: " + ProteinList.Count);
            sb.AppendLine("Fixed mods: " + DescribeModifications(FixedModifications));
            sb.AppendLine("Variable mods: " + DescribeModifications(VariableModifications));
            sb.AppendLine("Silac labels: " + DescribeSilacLabels(SilacLabels));
            sb.AppendLine("Turnover labels: " + (TurnoverLabels == null
                ? "none"
                : DescribeSilacLabel(TurnoverLabels.Value.StartLabel) + ">" + DescribeSilacLabel(TurnoverLabels.Value.EndLabel)));
            sb.AppendLine("Dissociation Type: " + CommonParameters.DissociationType);
            sb.AppendLine("Contaminant Handling: " + TcAmbiguity);

            sb.AppendLine("protease: " + CommonParameters.DigestionParams.DigestionAgent);
            if (CommonParameters.DigestionParams is DigestionParams digestionParams)
                sb.AppendLine("initiatorMethionineBehavior: " + digestionParams.InitiatorMethionineBehavior);
            sb.AppendLine("maximumMissedCleavages: " + CommonParameters.DigestionParams.MaxMissedCleavages);
            sb.AppendLine("minPeptideLength: " + CommonParameters.DigestionParams.MinLength);
            sb.AppendLine("maxPeptideLength: " + CommonParameters.DigestionParams.MaxLength);
            sb.AppendLine("maximumVariableModificationIsoforms: " + CommonParameters.DigestionParams.MaxModificationIsoforms);
            sb.AppendLine("digestionTerminus: " + CommonParameters.DigestionParams.FragmentationTerminus);
            sb.AppendLine("maxModsForEachPeptide: " + CommonParameters.DigestionParams.MaxMods);
            sb.AppendLine("cleavageSpecificity: " + CommonParameters.DigestionParams.SearchModeType);
            if (CommonParameters.DigestionParams is DigestionParams digestionParam)
                sb.AppendLine("specificProtease: " + digestionParam.SpecificProtease);
            sb.AppendLine("maximumFragmentSize" + (int)Math.Round(MaxFragmentSize));

            sb.Append("Localizeable mods: " + ProteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Sum());
            return sb.ToString();
        }

        /// <summary>
        /// The identity of a modification list, for the index cache key.
        /// </summary>
        /// <remarks>
        /// This names every modification instead of counting them. Counting was the defect: the cache key
        /// is compared as literal text by <c>MetaMorpheusTask.SameSettings</c>, so two searches whose
        /// modification lists differed but happened to be the same length produced the same key, and the
        /// second search silently reused a peptide index built for the first one's modifications. The
        /// index stores decorated <c>PeptideWithSetModifications</c>, so that index is wrong rather than
        /// merely incomplete: a variable-mod swap loses every peptide bearing the new modification, and a
        /// fixed-mod swap puts wrong masses in every entry.
        ///
        /// Each entry is the readable id followed by a short hash of the modification's complete definition
        /// (<see cref="Modification.ToString"/>, i.e. its mods.txt record: target, location restriction,
        /// chemical formula, monoisotopic mass, neutral losses, diagnostic ions). The id alone would miss
        /// a user who edited a custom modification file in place, keeping the name and changing the mass;
        /// the hash covers every field, and keeps covering fields added to Modification later.
        ///
        /// The list order is deliberately preserved rather than sorted. Modification enumeration is
        /// truncated at MaxModificationIsoforms by a yield break in
        /// <c>ProteolyticPeptide.GetModifiedPeptides</c>, so which isoforms survive can depend on the
        /// order the modifications arrive in. Sorting here would let two orders share one key and reuse
        /// each other's index. Preserving order can only cost a needless rebuild, never a wrong reuse.
        /// </remarks>
        internal static string DescribeModifications(IEnumerable<Modification> modifications)
        {
            if (modifications == null)
            {
                return "none";
            }

            return string.Join(",", modifications.Select(m =>
                (m?.IdWithMotif ?? "unnamed") + "[" + ShortHash(m?.ToString()) + "]"));
        }

        /// <summary>
        /// The identity of the SILAC label list, for the index cache key. Labels reach
        /// <c>Protein.Digest</c> alongside the modifications and change which peptides land in the index,
        /// but appeared nowhere in the key. <see cref="SilacLabel"/> does not override ToString, so the
        /// fields are named here. Order is preserved for the same reason as the modifications.
        /// </summary>
        internal static string DescribeSilacLabels(IEnumerable<SilacLabel> labels)
        {
            if (labels == null)
            {
                return "none";
            }

            return string.Join(",", labels.Select(DescribeSilacLabel));
        }

        internal static string DescribeSilacLabel(SilacLabel label)
        {
            if (label == null)
            {
                return "none";
            }

            var sb = new StringBuilder();
            sb.Append(label.OriginalAminoAcid).Append('>').Append(label.AminoAcidLabel)
                .Append('(').Append(label.LabelChemicalFormula).Append(',').Append(label.MassDifference).Append(')');

            if (label.AdditionalLabels != null)
            {
                foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                {
                    sb.Append('+').Append(DescribeSilacLabel(additionalLabel));
                }
            }

            return sb.ToString();
        }

        /// <summary>
        /// A short, stable hash of a definition string, for use inside the index cache key.
        /// </summary>
        /// <remarks>
        /// Line endings are normalised first. <see cref="Modification.ToString"/> builds its text with
        /// AppendLine, which emits the writing machine's line ending, so the same modification would
        /// otherwise hash differently on Windows and Linux and force a rebuild on nothing.
        ///
        /// Eight bytes separates the handful of modifications in a search many times over, and keeps the
        /// key readable next to the id it qualifies. This is a cache key, not a security boundary.
        /// </remarks>
        internal static string ShortHash(string definition)
        {
            if (string.IsNullOrEmpty(definition))
            {
                return "none";
            }

            string normalized = definition.Replace("\r\n", "\n").Replace('\r', '\n');
            byte[] hash = SHA256.HashData(Encoding.UTF8.GetBytes(normalized));
            return Convert.ToHexString(hash, 0, 8).ToLowerInvariant();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;

            // digest database
            List<PeptideWithSetModifications> peptides = new List<PeptideWithSetModifications>();

            if (CommonParameters.DigestionParams is not DigestionParams digestionParams)
                throw new MetaMorpheusException("Digestion parameters must be of type DigestionParams. Not yet implemented for Rna Digestion");
            
            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (i) =>
            {
                List<PeptideWithSetModifications> localPeptides = new List<PeptideWithSetModifications>();

                for (; i < ProteinList.Count; i += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    localPeptides.AddRange(ProteinList[i].Digest(digestionParams, FixedModifications, VariableModifications, SilacLabels, TurnoverLabels));

                    progress++;
                    var percentProgress = (int)((progress / ProteinList.Count) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins...", NestedIds));
                    }
                }

                lock (peptides)
                {
                    peptides.AddRange(localPeptides);
                }
            });

            // sort peptides by mass
            peptides.Sort((x, y) => x.MonoisotopicMass.CompareTo(y.MonoisotopicMass));

            //create precursor index (if specified)
            List<int>[] precursorIndex = null;
            if (GeneratePrecursorIndex)
            {
                precursorIndex = CreateNewPrecursorIndex(peptides);
            }
            bool addInteriorTerminalModsToPrecursorIndex = GeneratePrecursorIndex && CommonParameters.DigestionParams.DigestionAgent.Name.Contains("single");
            List<Modification> terminalModifications = addInteriorTerminalModsToPrecursorIndex ?
                NonSpecificEnzymeSearchEngine.GetVariableTerminalMods(CommonParameters.DigestionParams.FragmentationTerminus, VariableModifications) :
                null;

            // create fragment index
            List<int>[] fragmentIndex;

            try
            {
                fragmentIndex = new List<int>[(int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1];
            }
            catch (OutOfMemoryException)
            {
                throw new MetaMorpheusException("Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }

            // populate fragment index
            progress = 0;
            oldPercentProgress = 0;
            List<Product> fragments = new List<Product>();

            for (int peptideId = 0; peptideId < peptides.Count; peptideId++)
            {
                peptides[peptideId].Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, fragments, CommonParameters.FragmentationParameters);

                foreach (var theoreticalFragment in fragments)
                {
                    double theoreticalFragmentMass = theoreticalFragment.NeutralMass;

                    //if low res round
                    if (CommonParameters.DissociationType == MassSpectrometry.DissociationType.LowCID)
                    {
                        theoreticalFragmentMass = Math.Round(theoreticalFragmentMass / 1.0005079, 0) * 1.0005079;
                    }

                    if (theoreticalFragmentMass < MaxFragmentSize && theoreticalFragmentMass > 0)
                    {
                        int fragmentBin = (int)Math.Round(theoreticalFragmentMass * FragmentBinsPerDalton);

                        if (fragmentIndex[fragmentBin] == null)
                        {
                            fragmentIndex[fragmentBin] = new List<int> { peptideId };
                        }
                        else
                        {
                            fragmentIndex[fragmentBin].Add(peptideId);
                        }
                    }
                }

                //Add terminal mods if needed (do it here rather than earlier so that we don't have to fragment twice)
                if (addInteriorTerminalModsToPrecursorIndex)
                {
                    AddInteriorTerminalModsToPrecursorIndex(precursorIndex, fragments, peptides[peptideId], peptideId, terminalModifications);
                }

                progress++;
                var percentProgress = (int)((progress / peptides.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Fragmenting peptides...", NestedIds));
                }
            }

            return new IndexingResults(peptides, fragmentIndex, precursorIndex, this);
        }

        private List<int>[] CreateNewPrecursorIndex(List<PeptideWithSetModifications> peptidesSortedByMass)
        {
            // create precursor index
            List<int>[] precursorIndex = null;
            try
            {
                precursorIndex = new List<int>[(int)Math.Ceiling(MaxFragmentSize) * FragmentBinsPerDalton + 1];
            }
            catch (OutOfMemoryException)
            {
                throw new MetaMorpheusException("Max precursor mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }

            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(0, "Creating precursor index...", NestedIds));

            //Add all the precursors
            for (int i = 0; i < peptidesSortedByMass.Count; i++)
            {
                double mass = peptidesSortedByMass[i].MonoisotopicMass;
                if (!double.IsNaN(mass))
                {
                    if (mass > MaxFragmentSize) //if the precursor is larger than the index allows, then stop adding precursors
                    {
                        break;
                    }

                    int precursorBin = (int)Math.Round(mass * FragmentBinsPerDalton);

                    if (precursorIndex[precursorBin] == null)
                    {
                        precursorIndex[precursorBin] = new List<int> { i };
                    }
                    else
                    {
                        precursorIndex[precursorBin].Add(i);
                    }
                }
                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Creating precursor index...", NestedIds));
                }
            }
            return precursorIndex;
        }

        //add possible protein/peptide terminal modifications that aren't on the terminal amino acids
        //The purpose is for terminal mods that are contained WITHIN the Single peptide
        private void AddInteriorTerminalModsToPrecursorIndex(List<int>[] precursorIndex, List<Product> fragmentMasses, PeptideWithSetModifications peptide, int peptideId, List<Modification> variableModifications)
        {
            //Get database annotated mods
            Dictionary<int, List<Modification>> databaseAnnotatedMods = NonSpecificEnzymeSearchEngine.GetTerminalModPositions(peptide, CommonParameters.DigestionParams, variableModifications);
            foreach (KeyValuePair<int, List<Modification>> relevantDatabaseMod in databaseAnnotatedMods)
            {
                int fragmentNumber = relevantDatabaseMod.Key;
                Product fragmentAtIndex = fragmentMasses.FirstOrDefault(x => x.FragmentNumber == fragmentNumber);

                double basePrecursorMass;
                if (fragmentAtIndex is null)
                {
                    basePrecursorMass = peptide.MonoisotopicMass;
                }
                else
                {
                    basePrecursorMass = fragmentAtIndex.NeutralMass -
                                        DissociationTypeCollection.GetMassShiftFromProductType(fragmentAtIndex.ProductType) +
                                        WaterMonoisotopicMass;
                }

                foreach (Modification mod in relevantDatabaseMod.Value)
                {
                    double modifiedMass = basePrecursorMass + mod.MonoisotopicMass.Value;
                    if (modifiedMass <= MaxFragmentSize) //if the precursor is larger than the index allows, then don't add it
                    {
                        int precursorBin = (int)Math.Round(modifiedMass * FragmentBinsPerDalton);

                        if (precursorIndex[precursorBin] == null)
                        {
                            precursorIndex[precursorBin] = new List<int> { peptideId };
                        }
                        else
                        {
                            precursorIndex[precursorBin].Add(peptideId);
                        }
                    }
                }
            }
        }
    }
}