using OldInternalLogic;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    public class IndexEngine : MyEngine
    {
        public List<MorpheusModification> fixedModifications { get; private set; }
        public List<Protein> proteinList { get; private set; }
        public List<MorpheusModification> localizeableModifications { get; private set; }
        public Protease protease { get; private set; }
        public List<MorpheusModification> variableModifications { get; private set; }

        public IndexEngine(List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease) : base(2)
        {
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.protease = protease;
        }


        protected override MyResults RunSpecific()
        {
            var myDictionary = new List<CompactPeptide>();
            var myFragmentDictionary = new Dictionary<float, List<int>>(100000);
            int numProteins = 0;
            int totalProteins = proteinList.Count;
            HashSet<string> level3_observed = new HashSet<string>();
            HashSet<string> level4_observed = new HashSet<string>();

            var lp = new List<ProductType>() { ProductType.b, ProductType.y };
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<float, List<int>> myInnerDictionary = new Dictionary<float, List<int>>(100000);
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable).ToList();
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

                        peptide.SetFixedModifications(fixedModifications);

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, 4098, 3, localizeableModifications).ToList();
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

                            var ps = new CompactPeptide(yyy, variableModifications, localizeableModifications);

                            int index;
                            lock (myDictionary)
                            {
                                index = myDictionary.Count;
                                myDictionary.Add(ps);
                            }

                            foreach (var huhu in yyy.FastSortedProductMasses(lp))
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
                    //if (numProteins % 100 == 0)
                    //    output("Proteins: " + numProteins + " / " + totalProteins);
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

            //output("finished generating peptide index");

            return new IndexResults(myDictionary, myFragmentDictionary, this);
        }

        public override void ValidateParams()
        {
            throw new NotImplementedException();
        }
    }
}