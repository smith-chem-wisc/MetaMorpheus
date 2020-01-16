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

        //The modification problem is turned into a Directed Acyclic Graph. The Graph was build with matrix, and dynamic programming is used.
        //The function goes through the AdjNode[][] array from left to right, assign weight to each AdjNode, keep track of the heaviest previous AdjNode.
        public void LocalizeOGlycan(int[] modPos, GlycanBox glycanBox, GlycanBox[] childBoxes, HashSet<int> allPeaks, List<Product> products, int peptideLength)
        {
            var boxSatisfyBox = BoxSatisfyBox(childBoxes);

            for (int i = 0; i < modPos.Length; i++)
            {
                //maxLength: the most mods we can have up to current mod pos; minlengtt: the least mods we can have up to current mod pos.
                int maxLength = i + 1;
                int minlength = glycanBox.ModIds.Length - (modPos.Length - 1 - i);

                for (int j = 0; j < childBoxes.Length; j++)
                {
                    if (childBoxes[j].NumberOfMods <= maxLength && childBoxes[j].NumberOfMods >= minlength)
                    {
                        AdjNode adjNode = new AdjNode(i, j, modPos[i], childBoxes[j]);
                        var cost = LocalizationGraph.CalculateCost(allPeaks, products, peptideLength, modPos, i, glycanBox, childBoxes[j], 1000);

                        //The first line of the graph didnot have Sources.
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
                                //Check if a previous AdjNode exist and the current AdjNode could link to previous AdjNode. 
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

        #region LocalizeMod not limited to OGlycan.
        //Tt is possible to Merge this function to LocalizdOGlycan; but there is possible no need to do that.

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

        //For current ModPos at Ind, is the childbox satify the condition.
        //The function is for ModBox contains Mod that have different motif. 
        public static bool BoxSatisfyModPos(ModBox totalBox, ModBox childBox, int Ind, PeptideWithSetModifications peptide)
        {
            //Satisfy left
            foreach (var mn in childBox.MotifNeeded)
            {
                List<int> possibleModSites = new List<int>();

                ModificationMotif.TryGetMotif(mn.Key, out ModificationMotif motif);
                Modification modWithMotif = new Modification(_target: motif, _locationRestriction: "Anywhere.");

                for (int r = 0; r < Ind - 1; r++)
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
            foreach (var iy in childBox.ModIds)
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

                for (int r = Ind - 1; r < peptide.Length; r++)
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


        #endregion

        //Check if array1 contains array2 with repeats numbers.
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

        //The function defines how a childBox could be linked from all childBoxes.
        public static Dictionary<int, bool[]> BoxSatisfyBox(ModBox[] childBoxes)
        {
            Dictionary<int, bool[]> boxIdBoxes = new Dictionary<int, bool[]>();
            for (int i = 0; i < childBoxes.Length; i++)
            {
                bool[] idBoxes = new bool[childBoxes.Length];
                for (int j = 0; j <= i; j++)
                {
                    if (childBoxes[i].NumberOfMods <= childBoxes[j].NumberOfMods + 1 && (childBoxes[j].NumberOfMods ==0 || TryGetLeft(childBoxes[i].ModIds, childBoxes[j].ModIds)))
                    {
                        idBoxes[j] = true;
                    }
                }
                boxIdBoxes.Add(i, idBoxes);
            }

            return boxIdBoxes;
        }

        //Get path of Directed Acyclic Graph by recursion. 
        //Start from the last AdjNode[row-1 ][col-1], go back to it Sources, which contains the previous AdjNode with the highest cost.
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

        //The original path we get is just an array of AdjNode positions. This function here is to transfer the path into localized path. 
        //The output note: Tuple<(mod site)int, (glycanId)int>[glycanBox.Count] 
        //Basicly,  any change from left to right of the path indicates a modification. For example, the path = [1, 1, 2, 2] which means there is a modification at path[0] and path[3]
        public static Tuple<int, int>[] GetLocalizedPath(AdjNode[][] array, int[] modPos, ModBox[] childBoxes, int[] path)
        {
            int length = modPos.Length - 1;

            List<Tuple<int, int>> tuples = new List<Tuple<int, int>>();
            //Add first mod. If the childBoxes[path[0]].ModIds.Count == 0, means this is an empty childBox. 
            //Otherwise childBoxes[path[0]].ModIds.Count == 1 and childBoxes[path[0]].ModIds only contains one ModId.
            if (childBoxes[path[0]].ModIds.Count() != 0)
            {                
                tuples.Add(new Tuple<int, int>(modPos[0], childBoxes[path[0]].ModIds.First()));
            }

            for (int i = 1; i < path.Length; i++)
            {
                //If there is a change of the path, get the difference between the two Adjnodes of the array.
                if (path[i] != path[i - 1])
                {
                    var left = GetLeft(array[i][path[i]].ModBox.ModIds, array[i - 1][path[i - 1]].ModBox.ModIds).First();

                    tuples.Add(new Tuple<int, int>(modPos[i], left));
                }
            }

            return tuples.ToArray();
        }

        //Get the difference between array 1 and array 2 with repeat numbers.
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
