using Chemistry;
using EngineLayer.NonSpecificEnzymeSearch;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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

        public IndexingEngine(List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications,
            List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel, int currentPartition, DecoyType decoyType,
            CommonParameters commonParams, double maxFragmentSize, bool generatePrecursorIndex, List<FileInfo> proteinDatabases, List<string> nestedIds)
            : base(commonParams, nestedIds)
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
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Databases: " + string.Join(",", ProteinDatabases.OrderBy(p => p.Name).Select(p => p.Name + ":" + p.CreationTime)));
            sb.AppendLine("Partitions: " + CurrentPartition + "/" + CommonParameters.TotalPartitions);
            sb.AppendLine("Precursor Index: " + GeneratePrecursorIndex);
            sb.AppendLine("Search Decoys: " + DecoyType);
            sb.AppendLine("Number of proteins: " + ProteinList.Count);
            sb.AppendLine("Number of fixed mods: " + FixedModifications.Count);
            sb.AppendLine("Number of variable mods: " + VariableModifications.Count);
            sb.AppendLine("Dissociation Type: " + CommonParameters.DissociationType);

            sb.AppendLine("protease: " + CommonParameters.DigestionParams.Protease);
            sb.AppendLine("initiatorMethionineBehavior: " + CommonParameters.DigestionParams.InitiatorMethionineBehavior);
            sb.AppendLine("maximumMissedCleavages: " + CommonParameters.DigestionParams.MaxMissedCleavages);
            sb.AppendLine("minPeptideLength: " + CommonParameters.DigestionParams.MinPeptideLength);
            sb.AppendLine("maxPeptideLength: " + CommonParameters.DigestionParams.MaxPeptideLength);
            sb.AppendLine("maximumVariableModificationIsoforms: " + CommonParameters.DigestionParams.MaxModificationIsoforms);
            sb.AppendLine("digestionTerminus: " + CommonParameters.DigestionParams.FragmentationTerminus);
            sb.AppendLine("maxModsForEachPeptide: " + CommonParameters.DigestionParams.MaxModsForPeptide);
            sb.AppendLine("cleavageSpecificity: " + CommonParameters.DigestionParams.SearchModeType);
            sb.AppendLine("specificProtease: " + CommonParameters.DigestionParams.SpecificProtease);
            sb.AppendLine("maximumFragmentSize" + (int)Math.Round(MaxFragmentSize));

            sb.Append("Localizeable mods: " + ProteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Sum());
            return sb.ToString();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;

            // digest database
            List<PeptideWithSetModifications> globalPeptides = new List<PeptideWithSetModifications>();

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (i) =>
            {
                List<PeptideWithSetModifications> localPeptides = new List<PeptideWithSetModifications>();

                for (; i < ProteinList.Count; i += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    localPeptides.AddRange(ProteinList[i].Digest(CommonParameters.DigestionParams, FixedModifications, VariableModifications, SilacLabels, TurnoverLabels));

                    progress++;
                    var percentProgress = (int)((progress / ProteinList.Count) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Digesting proteins...", NestedIds));
                    }
                }

                lock (globalPeptides)
                {
                    globalPeptides.AddRange(localPeptides);
                }
            });

            // sort peptides by mass
            List<PeptideWithSetModifications> peptidesSortedByMass = globalPeptides.OrderBy(p => p.MonoisotopicMass).ToList();
            globalPeptides = null; //RAM

            //create precursor index (if specified)
            List<int>[] precursorIndex = null;
            if (GeneratePrecursorIndex)
            {
                precursorIndex = CreateNewPrecursorIndex(peptidesSortedByMass);
            }
            bool addInteriorTerminalModsToPrecursorIndex = GeneratePrecursorIndex && CommonParameters.DigestionParams.Protease.Name.Contains("single");
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
            for (int peptideId = 0; peptideId < peptidesSortedByMass.Count; peptideId++)
            {
                List<Product> fragments = peptidesSortedByMass[peptideId].Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus).ToList();
                List<double> fragmentMasses = fragments.Select(m => m.NeutralMass).ToList();

                foreach (double theoreticalFragmentMass in fragmentMasses)
                {
                    double tfm = theoreticalFragmentMass;
                    //if low res round
                    if (CommonParameters.DissociationType == MassSpectrometry.DissociationType.LowCID)
                    {
                        tfm = Math.Round(theoreticalFragmentMass / 1.0005079, 0) * 1.0005079;
                    }

                    if (tfm < MaxFragmentSize && tfm > 0)
                    {
                        int fragmentBin = (int)Math.Round(tfm * FragmentBinsPerDalton);

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

                //Add terminal mods if needed (do it here rather than earlier so that we don't have to fragment twice
                if (addInteriorTerminalModsToPrecursorIndex)
                {
                    AddInteriorTerminalModsToPrecursorIndex(precursorIndex, fragments, peptidesSortedByMass[peptideId], peptideId, terminalModifications);
                }

                progress++;
                var percentProgress = (int)((progress / peptidesSortedByMass.Count) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Fragmenting peptides...", NestedIds));
                }
            }

            return new IndexingResults(peptidesSortedByMass, fragmentIndex, precursorIndex, this);
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
            int progress = 0;
            int oldPercentProgress = 0;
            int maxProgress = peptidesSortedByMass.Count;
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
                var percentProgress = (int)((progress / maxProgress) * 100);

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
                Product fragmentAtIndex = fragmentMasses.Where(x => x.TerminusFragment.FragmentNumber == fragmentNumber).FirstOrDefault();
                double basePrecursorMass = fragmentAtIndex == null ? peptide.MonoisotopicMass : fragmentAtIndex.NeutralMass - DissociationTypeCollection.GetMassShiftFromProductType(fragmentAtIndex.ProductType) + WaterMonoisotopicMass;

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