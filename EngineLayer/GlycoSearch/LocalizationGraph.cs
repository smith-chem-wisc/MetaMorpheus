using System;
using System.Collections.Generic;
using System.Linq;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics;

namespace EngineLayer.GlycoSearch
{
    public class LocalizationGraph
    {
        public AdjNode[][] array;
        public LocalizationGraph(int posCount, int glycanCount)
        {
            array = new AdjNode[posCount][];
            for (int i = 0; i < posCount; i++)
            {
                array[i] = new AdjNode[glycanCount];
            }
        }

        public void LocalizeOGlycan(int[] modPos, GlycanBox glycanBox, GlycanBox[] boxes, HashSet<int> allPeaks, List<Product> products, int peptideLength)
        {
            var boxSatisfyBox = BoxSatisfyBox(boxes);

            for (int i = 0; i < modPos.Length; i++)
            {
                int maxLength = i + 1;
                int minlength = glycanBox.ModIds.Length - (modPos.Length - 1 - i);

                for (int j = 0; j < boxes.Length; j++)
                {
                    if (boxes[j].NumberOfMods <= maxLength && boxes[j].NumberOfMods >= minlength)
                    {
                        AdjNode adjNode = new AdjNode(i, j, modPos[i], boxes[j]);
                        var cost = LocalizationGraph.CalculateCost(allPeaks, products, peptideLength, modPos, i, glycanBox, boxes[j], 1000);
                        if (i == 0)
                        {
                            //Get cost                             
                            adjNode.maxCost = cost;

                        }
                        else
                        {
                            double maxCost = 0;
                            for (int prej = 0; prej <= j; prej++)
                            {
                                //TO DO: The condition here could be wrong, please change to function BoxSatisfyBox.
                                if (boxSatisfyBox[j][prej] && array[i - 1][prej] != null)
                                {
                                    var tempCost = cost + array[i - 1][prej].maxCost;
                                    if (tempCost > maxCost)
                                    {
                                        adjNode.Sources.Clear();
                                        adjNode.Costs.Clear();

                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(tempCost);
                                        maxCost = tempCost;
                                    }
                                    else if (tempCost == maxCost)
                                    {
                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(tempCost);
                                    }

                                }
                            }
                            //if (adjNode.Costs.Any())
                            {
                                adjNode.maxCost = adjNode.Costs.Max();
                            }
                        }

                        array[i][j] = adjNode;
                    }
                }

            }
        }

        public static double CalculateCost(HashSet<int> allPeaksForLocalization, List<Product> products, int peptideLength, int[] modPos, int modInd, GlycanBox OGlycanBox, GlycanBox box, int FragmentBinsPerDalton)
        {
            if (modInd == modPos.Length - 1)
            {
                return 0;
            }

            var fragmentHash = GlycoPeptides.GetLocalFragmentHash(products, peptideLength, modPos, modInd, OGlycanBox, box, FragmentBinsPerDalton);

            int currentLocalizationScore = allPeaksForLocalization.Intersect(fragmentHash).Count();

            return (double)currentLocalizationScore;
        }

        //The modification problem is turned into a Directed Acyclic Graph. The Graph was build with matrix, and dynamic programming is used.
        public void LocalizeMod(int[] modPos, ModBox totalBox, ModBox[] boxes, HashSet<int> allPeaks, List<Product> products, PeptideWithSetModifications peptide)
        {
            var boxSatisfyBox = BoxSatisfyBox(boxes);

            for (int i = 0; i < modPos.Length; i++)
            {

                for (int j = 0; j < boxes.Length; j++)
                {
                    if (BoxSatisfyModPos(totalBox, boxes[j], modPos[i], peptide))
                    {
                        AdjNode adjNode = new AdjNode(i, j, modPos[i], boxes[j]);
                        var cost = LocalizationGraph.CalculateCostMod(allPeaks, products, peptide.Length, modPos, i, totalBox, boxes[j], 1000);
                        if (i == 0)
                        {
                            //Get cost                             
                            adjNode.maxCost = cost;

                        }
                        else
                        {
                            double maxCost = 0;
                            for (int prej = 0; prej <= j; prej++)
                            {
                                if (boxSatisfyBox[j][prej] && array[i - 1][prej] != null)
                                {
                                    var tempCost = cost + array[i - 1][prej].maxCost;
                                    if (tempCost > maxCost)
                                    {
                                        adjNode.Sources.Clear();
                                        adjNode.Costs.Clear();

                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(tempCost);
                                        maxCost = tempCost;
                                    }
                                    else if (tempCost == maxCost)
                                    {
                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(tempCost);
                                    }

                                }
                            }
                            //if (adjNode.Costs.Any())
                            {
                                adjNode.maxCost = adjNode.Costs.Max();
                            }
                        }

                        array[i][j] = adjNode;
                    }
                }

            }
        }

        private static bool TryGetLeft(int[] array1, int[] array2)
        {
            //Get compliment box
            var gx = array1.GroupBy(p => p).ToDictionary(p => p.Key, p => p.ToList());
            foreach (var iy in array2)
            {
                if (!gx.ContainsKey(iy))
                {
                    return false;
                }
                else if (gx[iy].Count == 0)
                {
                    return false;
                }
                gx[iy].RemoveAt(gx[iy].Count - 1);
            }
            return true;
        }

        public static Dictionary<int, bool[]> BoxSatisfyBox(ModBox[] boxes)
        {
            Dictionary<int, bool[]> boxIdBoxes = new Dictionary<int, bool[]>();
            for (int i = 0; i < boxes.Length; i++)
            {
                bool[] idBoxes = new bool[boxes.Length];
                for (int j = 0; j <= i; j++)
                {
                    if (boxes[i].NumberOfMods <= boxes[j].NumberOfMods + 1 && (boxes[j].NumberOfMods ==0 || TryGetLeft(boxes[i].ModIds, boxes[j].ModIds)))
                    {
                        idBoxes[j] = true;
                    }
                }
                boxIdBoxes.Add(i, idBoxes);
            }

            return boxIdBoxes;
        }

        public static double CalculateCostMod(HashSet<int> allPeaksForLocalization, List<Product> products, int peptideLength, int[] modPos, int modInd, ModBox totalBox, ModBox localBox, int FragmentBinsPerDalton)
        {
            if (modInd == modPos.Length - 1)
            {
                return 0;
            }

            var localFragmentHash = ModBox.GetLocalFragmentHash(products, peptideLength, modPos, modInd, totalBox, localBox, FragmentBinsPerDalton);
 
            int currentLocalizationScore = allPeaksForLocalization.Intersect(localFragmentHash).Count();

            return (double)currentLocalizationScore;
        }

        //For current ModPos at Ind, is the child boxsatify the condition.
        public static bool BoxSatisfyModPos(ModBox totalBox, ModBox box, int Ind, PeptideWithSetModifications peptide)
        {
            //Satisfy left
            foreach (var mn in box.MotifNeeded)
            {
                List<int> possibleModSites = new List<int>();

                ModificationMotif.TryGetMotif(mn.Key, out ModificationMotif motif);
                Modification modWithMotif = new Modification(_target: motif, _locationRestriction: "Anywhere.");

                for (int r = 0; r < Ind-1; r++)
                {
                    if (peptide.AllModsOneIsNterminus.Keys.Contains(r + 2))
                    {
                        continue;
                    }

                    if (ModificationLocalization.ModFits(modWithMotif, peptide.BaseSequence, r + 1, peptide.Length, r + 1))
                    {
                        possibleModSites.Add(r + 2);
                    }
                }

                if (possibleModSites.Count < mn.Value.Count)
                {
                    return false;
                }
            }

            //Get compliment box
            var gx = totalBox.ModIds.GroupBy(p => p).ToDictionary(p => p.Key, p => p.ToList());
            foreach (var iy in box.ModIds)
            {
                gx[iy].RemoveAt(gx[iy].Count - 1);
            }
            var left = gx.SelectMany(p => p.Value).ToArray();

            var complimentBox = new ModBox(left.ToArray());

            //Satify right
            foreach (var mn in complimentBox.MotifNeeded)
            {
                List<int> possibleModSites = new List<int>();

                ModificationMotif.TryGetMotif(mn.Key, out ModificationMotif motif);
                Modification modWithMotif = new Modification(_target: motif, _locationRestriction: "Anywhere.");

                for (int r = Ind-1; r < peptide.Length; r++)
                {
                    if (peptide.AllModsOneIsNterminus.Keys.Contains(r + 2))
                    {
                        continue;
                    }

                    if (ModificationLocalization.ModFits(modWithMotif, peptide.BaseSequence, r + 1, peptide.Length, r + 1))
                    {
                        possibleModSites.Add(r + 2);
                    }
                }

                if (possibleModSites.Count < mn.Value.Count)
                {
                    return false;
                }
            }

            return true;
        }

        public static int[] GetLeft(int[] array1, int[] array2)
        {
            //Get compliment box
            var gx = array1.GroupBy(p => p).ToDictionary(p => p.Key, p => p.ToList());
            foreach (var iy in array2)
            {
                gx[iy].RemoveAt(gx[iy].Count - 1);
            }
            var left = gx.SelectMany(p => p.Value).ToArray();
            return left;
        }

        public static Tuple<int, int>[] GetLocalizedPeptide(AdjNode[][] array, int[] modPos, ModBox[] boxes, int[] indexes)
        {           
            int length = modPos.Length - 1;

            List<Tuple<int, int>> tuples = new List<Tuple<int, int>>();
            //Add first.
            if (boxes[indexes[0]].ModIds.Count()!=0)
            {
                tuples.Add(new Tuple<int, int>(modPos[0], boxes[indexes[0]].ModIds.First()));
            }

            for (int i = 1; i < indexes.Length; i++)
            {
                if (indexes[i] != indexes[i-1])
                {
                    var left = GetLeft(array[i][indexes[i]].ModBox.ModIds, array[i-1][indexes[i-1]].ModBox.ModIds).First();

                    tuples.Add(new Tuple<int, int>( modPos[i], left));
                }
            }

            return tuples.ToArray();
        }

        //Get path of Directed Acyclic Graph by recursion. 
        public static List<int[]> GetAllPaths(AdjNode[][] array, ModBox[] boxes)
        {
            List<int[]> allPaths = new List<int[]>();

            int xlength = array.Length;
            int ylength = array.First().Length;

            int[] temp = new int[xlength];

            temp[xlength - 1] = ylength -1;
            
            PathHelper(allPaths, array, xlength -1, ylength -1, temp);

            return allPaths;
        }

        private static void PathHelper(List<int[]> allPaths, AdjNode[][] array, int xind, int yind, int[] temp)
        {
            if (xind == 0)
            {
                allPaths.Add((int[])temp.Clone());
                return;
            }

            foreach (var pre in array[xind][yind].Sources)
            {
                xind--;
                yind = pre;
                temp[xind] = yind;
                PathHelper(allPaths, array, xind, yind, temp);

                xind++;
            }
        }
    }

    public class AdjNode
    {
        public AdjNode(int i, int j, int modPos, ModBox box)
        {
            PosX = i;
            PosY = j;
            ModPos = modPos;
            ModBox = box;
        }

        public int PosX { get; }
        public int PosY { get; }

        public int ModPos { get; }
        public ModBox ModBox { get;  }

        //sources are represented by index.
        public List<int> Sources = new List<int>();

        public List<double> Costs = new List<double>();

        public double maxCost { get; set; }

    }

}
