using MetaMorpheus;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace IndexSearchAndAnalyze
{
    public class IndexEngine : MyEngine
    {
        public IndexEngine(IndexParams indexParams)
        {
            this.myParams = indexParams;
        }

        protected override MyResults RunSpecific()
        {
            var indexParams = (IndexParams)myParams;
            var myDictionary = new List<CompactPeptide>();
            var myFragmentDictionary = new Dictionary<float, List<int>>(100000);
            int numProteins = 0;
            int totalProteins = indexParams.proteinList.Count;
            HashSet<string> level3_observed = new HashSet<string>();
            HashSet<string> level4_observed = new HashSet<string>();

            var lp = new List<ProductType>() { ProductType.b, ProductType.y };
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<float, List<int>> myInnerDictionary = new Dictionary<float, List<int>>(100000);
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = indexParams.proteinList[i];
                    var digestedList = protein.Digest(indexParams.protease, 2, InitiatorMethionineBehavior.Variable).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length == 1 || peptide.Length > 252)
                            continue;

                        if (peptide.OneBasedPossibleLocalizedModifications.Count == 0)
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

                        peptide.SetFixedModifications(indexParams.fixedModifications);

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(indexParams.variableModifications, 4098, 3, indexParams.localizeableModifications).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            if (peptide.OneBasedPossibleLocalizedModifications.Count > 0)
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

                            var ps = new CompactPeptide(yyy, indexParams.variableModifications, indexParams.localizeableModifications);

                            int index;
                            lock (myDictionary)
                            {
                                index = myDictionary.Count;
                                myDictionary.Add(ps);
                            }

                            foreach (var huhu in yyy.FastUnsortedProductMasses(lp))
                            {
                                float rounded = (float)Math.Round(huhu, 3);
                                List<int> value;
                                if (myInnerDictionary.TryGetValue(rounded, out value))
                                    value.Add(index);
                                else
                                    myInnerDictionary.Add(rounded, new List<int>() { index });
                            }
                            ps.MonoisotopicMass = (float)yyy.MonoisotopicMass;
                        }
                    }
                    numProteins++;
                    if (numProteins % 100 == 0)
                        Console.WriteLine("Proteins: " + numProteins + " / " + totalProteins);
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
                                myFragmentDictionary.Add(huhu.Key, new List<int>() { hhhh });
                        }
                    }
                }
            });

            Console.WriteLine("finished generating peptide index");

            return new IndexResults(myDictionary, myFragmentDictionary, indexParams);
        }
    }
}