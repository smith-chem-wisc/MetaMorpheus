using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.Indexing
{
    public class IndexingEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private const int decimalDigitsForFragmentMassRounding = 3;
        private readonly List<Protein> proteinList;

        private readonly List<ModificationWithMass> fixedModifications;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly List<ProductType> lp;
        private readonly int currentPartition;
        private readonly bool searchDecoys;
        private readonly CommonParameters CommonParameters;

        #endregion Private Fields

        #region Public Constructors

        public IndexingEngine(List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<ProductType> lp, int currentPartition, bool searchDecoys, CommonParameters CommonParameters, List<string> nestedIds) : base(nestedIds)
        {
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.lp = lp;
            this.currentPartition = currentPartition + 1;
            this.searchDecoys = searchDecoys;
            this.CommonParameters = CommonParameters;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Partitions: " + currentPartition + "/" + CommonParameters.TotalPartitions);
            sb.AppendLine("Search Decoys: " + searchDecoys);
            sb.AppendLine("Number of proteins: " + proteinList.Count);
            sb.AppendLine("Number of fixed mods: " + fixedModifications.Count);
            sb.AppendLine("Number of variable mods: " + variableModifications.Count);
            sb.AppendLine("lp: " + string.Join(",", lp));
            sb.AppendLine("protease: " + CommonParameters.Protease);
            sb.AppendLine("initiatorMethionineBehavior: " + CommonParameters.InitiatorMethionineBehavior);
            sb.AppendLine("maximumMissedCleavages: " + CommonParameters.MaxMissedCleavages);
            sb.AppendLine("minPeptideLength: " + CommonParameters.MinPeptideLength);
            sb.AppendLine("maxPeptideLength: " + CommonParameters.MaxPeptideLength);
            sb.AppendLine("maximumVariableModificationIsoforms: " + CommonParameters.MaxModificationIsoforms);
            sb.Append("Localizeable mods: " + proteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Sum());
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var myDictionary = new List<CompactPeptide>();
            var myFragmentDictionary = new Dictionary<float, List<int>>(100000);
            int totalProteins = proteinList.Count;
            var observed_sequences = new HashSet<CompactPeptide>();
            int proteinsSeen = 0;
            int old_progress = 0;
            TerminusType terminusType = ProductTypeToTerminusType.IdentifyTerminusType(lp);
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                var myInnerDictionary = new Dictionary<float, List<int>>(100000);
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(CommonParameters.Protease, (int)CommonParameters.MaxMissedCleavages, CommonParameters.MinPeptideLength, CommonParameters.MaxPeptideLength, CommonParameters.InitiatorMethionineBehavior, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, (int)CommonParameters.MaxModificationIsoforms, (int)CommonParameters.Max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            var correspondingCompactPeptide = yyy.CompactPeptide(terminusType);
                            var observed = observed_sequences.Contains(correspondingCompactPeptide);
                            if (observed)
                                continue;
                            lock (observed_sequences)
                            {
                                observed = observed_sequences.Contains(correspondingCompactPeptide);
                                if (observed)
                                    continue;
                                observed_sequences.Add(correspondingCompactPeptide);
                            }

                            int index;
                            lock (myDictionary)
                            {
                                index = myDictionary.Count;
                                myDictionary.Add(correspondingCompactPeptide);
                            }

                            foreach (var huhu in correspondingCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp))
                            {
                                if (!double.IsNaN(huhu))
                                {
                                    var rounded = (float)Math.Round(huhu, decimalDigitsForFragmentMassRounding);
                                    if (myInnerDictionary.TryGetValue(rounded, out List<int> value))
                                    {
                                        if (!value.Contains(index))
                                            value.Add(index);
                                    }
                                    else
                                        myInnerDictionary.Add(rounded, new List<int> { index });
                                }
                            }
                        }
                    }
                }
                lock (myFragmentDictionary)
                {
                    foreach (var huhu in myInnerDictionary)
                    {
                        foreach (var hhhh in huhu.Value)
                        {
                            if (myFragmentDictionary.TryGetValue(huhu.Key, out List<int> value))
                                value.Add(hhhh);
                            else
                                myFragmentDictionary.Add(huhu.Key, new List<int> { hhhh });
                        }
                    }
                    proteinsSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In indexing loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });

            return new IndexingResults(myDictionary, myFragmentDictionary, this);
        }

        #endregion Protected Methods
    }
}