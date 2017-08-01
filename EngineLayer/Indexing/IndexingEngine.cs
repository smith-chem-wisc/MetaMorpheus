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

        private const int max_mods_for_peptide = 3;
        private const int decimalDigitsForFragmentMassRounding = 3;
        private readonly int maximumMissedCleavages;
        private readonly int? minPeptideLength;
        private readonly int? maxPeptideLength;
        private readonly int maximumVariableModificationIsoforms;
        private readonly List<Protein> proteinList;

        private readonly Protease protease;

        private readonly List<ModificationWithMass> fixedModifications;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly InitiatorMethionineBehavior initiatorMethionineBehavior;

        private readonly List<ProductType> lp;
        private readonly bool addCompIons;

        #endregion Private Fields

        #region Public Constructors

        public IndexingEngine(List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, Protease protease, InitiatorMethionineBehavior initiatorMethionineBehavior, int maximumMissedCleavages, int? minPeptideLength, int? maxPeptideLength, int maximumVariableModificationIsoforms, List<ProductType> lp, List<string> nestedIds, bool addCompIons) : base(nestedIds)
        {
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.protease = protease;
            this.initiatorMethionineBehavior = initiatorMethionineBehavior;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.minPeptideLength = minPeptideLength;
            this.maxPeptideLength = maxPeptideLength;
            this.maximumVariableModificationIsoforms = maximumVariableModificationIsoforms;
            this.lp = lp;
            this.addCompIons = addCompIons;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Number of proteins: " + proteinList.Count);
            sb.AppendLine("Number of fixed mods: " + fixedModifications.Count);
            sb.AppendLine("Number of variable mods: " + variableModifications.Count);
            sb.AppendLine("lp: " + string.Join(",", lp));
            sb.AppendLine("protease: " + protease);
            sb.AppendLine("initiatorMethionineBehavior: " + initiatorMethionineBehavior);
            sb.AppendLine("maximumMissedCleavages: " + maximumMissedCleavages);
            sb.AppendLine("minPeptideLength: " + minPeptideLength);
            sb.AppendLine("maxPeptideLength: " + maxPeptideLength);
            sb.AppendLine("maximumVariableModificationIsoforms: " + maximumVariableModificationIsoforms);
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
            var level3_observed = new HashSet<string>();
            var level4_observed = new HashSet<string>();
            int proteinsSeen = 0;
            int old_progress = 0;
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                var myInnerDictionary = new Dictionary<float, List<int>>(100000);
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications, addCompIons).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.NumKnownPossibleLocMods == 0)
                        {
                            lock (level3_observed)
                            {
                                var hc = peptide.BaseLeucineSequence;
                                var observed = level3_observed.Contains(hc);
                                if (observed)
                                    continue;
                                level3_observed.Add(hc);
                            }
                        }

                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            if (peptide.NumKnownPossibleLocMods > 0)
                            {
                                lock (level4_observed)
                                {
                                    var hc = yyy.Sequence;
                                    var observed = level4_observed.Contains(hc);
                                    if (observed)
                                        continue;
                                    level4_observed.Add(hc);
                                }
                            }

                            var ps = new CompactPeptide(yyy, addCompIons);

                            int index;
                            lock (myDictionary)
                            {
                                index = myDictionary.Count;
                                myDictionary.Add(ps);
                            }

                            foreach (var huhu in ps.ProductMassesMightHaveDuplicatesAndNaNs(lp))
                            {
                                if (!double.IsNaN(huhu))
                                {
                                    var rounded = (float)Math.Round(huhu, decimalDigitsForFragmentMassRounding);
                                    List<int> value;
                                    if (myInnerDictionary.TryGetValue(rounded, out value))
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
                        List<int> value;
                        foreach (var hhhh in huhu.Value)
                        {
                            if (myFragmentDictionary.TryGetValue(huhu.Key, out value))
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