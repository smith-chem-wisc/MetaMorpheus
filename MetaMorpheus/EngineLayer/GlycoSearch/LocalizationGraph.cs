using System;
using System.Collections.Generic;
using System.Linq;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using MzLibUtil;

namespace EngineLayer.GlycoSearch
{
    public class LocalizationGraph
    {
        public AdjNode[][] array { get; set; }
        public int[] ModPos { get; }

        public int ModBoxId { get; }
        public ModBox ModBox { get; }
        public ModBox[] ChildModBoxes { get; set; }

        public double NoLocalCost{get; set;} //Note that we have node for each glycosite, the matched ions before the first node and after the last node is scored here.
        public double TotalScore { get; set; } //Total score is the score of matched ions that are used for localization. For O-glycan, it is the score of all matched c/zDot ions. 

        public LocalizationGraph(int[] modPos, ModBox modBox, ModBox[] childModBoxes, int id)
        {
            ModPos = modPos;
            ModBox = modBox;
            ModBoxId = id;
            ChildModBoxes = childModBoxes;

            //array is localization graph matrix. array is composed of 2d array of node. From left to right, node is build under a glycosite. From up to down, node is build for each child box.
            array = new AdjNode[modPos.Length][];
            for (int i = 0; i < modPos.Length; i++)
            {
                array[i] = new AdjNode[ChildModBoxes.Length];
            }
        }

        //The modification problem is turned into a Directed Acyclic Graph. The Graph was build with matrix, and dynamic programming is used.
        //The function goes through the AdjNode[][] array from left to right, assign weight to each AdjNode, keep track of the heaviest previous AdjNode.
        public static void LocalizeOGlycan(LocalizationGraph localizationGraph, Ms2ScanWithSpecificMass theScan, Tolerance productTolerance, List<Product> products)
        {
            var boxSatisfyBox = BoxSatisfyBox(localizationGraph.ChildModBoxes);

            for (int i = 0; i < localizationGraph.ModPos.Length; i++)
            {
                //maxLength: the most mods we can have up to current mod pos; minlengtt: the least mods we can have up to current mod pos.
                int maxLength = i + 1;
                int minlength = localizationGraph.ModBox.ModIds.Length - (localizationGraph.ModPos.Length - 1 - i);

                for (int j = 0; j < localizationGraph.ChildModBoxes.Length; j++)
                {
                    if (localizationGraph.ChildModBoxes[j].NumberOfMods <= maxLength && localizationGraph.ChildModBoxes[j].NumberOfMods >= minlength)
                    {
                        AdjNode adjNode = new AdjNode(i, j, localizationGraph.ModPos[i], localizationGraph.ChildModBoxes[j]);

                        double cost = 0;
                        if (i != localizationGraph.ModPos.Length - 1)
                        {              
                            var fragments = GlycoPeptides.GetLocalFragment(products, localizationGraph.ModPos, i, localizationGraph.ModBox, localizationGraph.ChildModBoxes[j]);
                            cost = CalculateCost(theScan, productTolerance, fragments);
                        }

                        adjNode.CurrentCost = cost;
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
                                if (boxSatisfyBox[j][prej] && localizationGraph.array[i - 1][prej] != null)
                                {
                                    adjNode.AllSources.Add(prej);

                                    var tempCost = cost + localizationGraph.array[i - 1][prej].maxCost;
                                    if (tempCost > maxCost)
                                    {
                                        adjNode.CummulativeSources.Clear();
                
                                        adjNode.CummulativeSources.Add(prej);
             
                                        maxCost = tempCost;
                                    }
                                    else if (tempCost == maxCost)
                                    {
                                        adjNode.CummulativeSources.Add(prej);              
                                    }

                                }
                            }

                             adjNode.maxCost = maxCost;

                        }

                        localizationGraph.array[i][j] = adjNode;
                    }
                }

            }

            var unlocalFragments = GlycoPeptides.GetUnlocalFragment(products, localizationGraph.ModPos, localizationGraph.ModBox);
            var noLocalScore = CalculateCost(theScan, productTolerance, unlocalFragments);
            localizationGraph.NoLocalCost = noLocalScore;
            localizationGraph.TotalScore = localizationGraph.array[localizationGraph.ModPos.Length - 1][localizationGraph.ChildModBoxes.Length - 1].maxCost + noLocalScore;
        }

        //Based on our implementation of Graph localization. We need to calculate cost between two nearby nodes (glycosites) 
        public static double CalculateCost(Ms2ScanWithSpecificMass theScan, Tolerance productTolerance, List<double> fragments)
        {
            double score = 0;

            foreach (var f in fragments)
            {
                var closestExperimentalMass = theScan.GetClosestExperimentalIsotopicEnvelope(f);

                // is the mass error acceptable?
                if (productTolerance.Within(closestExperimentalMass.MonoisotopicMass, f) && closestExperimentalMass.Charge <= theScan.PrecursorCharge)
                {
                    score += 1 + closestExperimentalMass.Peaks.Sum(p => p.intensity) / theScan.TotalIonCurrent;
                }
            }
            return score;
        }

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

        //The Directed Acyclic Graph is build from left to right. In the process, we need to know which node can linked to nodes from its left. 
        //Since node contains Childbox. We name this function as BoxSatisfyBox.
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

        //Get all path with hightest score of Directed Acyclic Graph by recursion. 
        //Start from the last AdjNode[row-1 ][col-1], go back to it Sources, which contains the previous AdjNode with the highest cost.
        public static List<int[]> GetAllHighestScorePaths(AdjNode[][] array, ModBox[] boxes)
        {
            List<int[]> allPaths = new List<int[]>();

            int xlength = array.Length;
            int ylength = array.First().Length;

            int[] temp = new int[xlength];

            temp[xlength - 1] = ylength -1;

            GetAllHighestScorePathHelper(allPaths, array, xlength -1, ylength -1, temp);

            return allPaths;
        }

        private static void GetAllHighestScorePathHelper(List<int[]> allPaths, AdjNode[][] array, int xind, int yind, int[] temp)
        {
            if (xind == 0)
            {
                allPaths.Add((int[])temp.Clone());
                return;
            }

            foreach (var pre in array[xind][yind].CummulativeSources)
            {
                xind--;
                yind = pre;
                temp[xind] = yind;
                GetAllHighestScorePathHelper(allPaths, array, xind, yind, temp);

                xind++;
            }
        }

        //Get one path of Directed Acyclic Graph by recursion.
        public static int[] GetFirstPath(AdjNode[][] array, ModBox[] boxes)
        {

            int xlength = array.Length;
            int ylength = array.First().Length;

            int[] temp = new int[xlength];

            temp[xlength - 1] = ylength - 1;

            FirstPathHelper(array, xlength - 1, ylength - 1, temp);

            return temp;
        }

        private static void FirstPathHelper(AdjNode[][] array, int xind, int yind, int[] temp)
        {
            if (xind == 0)
            {
                return;
            }

            var pre = array[xind][yind].CummulativeSources.First();
            xind--;
            yind = pre;
            temp[xind] = yind;
            FirstPathHelper(array, xind, yind, temp);
        }

        //The original path we get is just an array of AdjNode positions. For example, path = [1, 1, 2, 2] means the best nodes are at array[0][1], array[1][1], array[2][2], array[3][2]
        //This function here is to transfer the path into localized Route. Route contains each glycosite with glycanId.
        //Basicly, any change from left to right of the path indicates a modification. For example, the path = [1, 1, 2, 2] which means there is a modification at ModPos[0] and ModPos[2]
        public static Route GetLocalizedPath(LocalizationGraph localizationGraph, int[] path)
        {
            Route route = new Route();

            if (path.Length == 1)
            {
                bool onlyOneLocalized = false;
                if (localizationGraph.TotalScore > 0)
                {
                    onlyOneLocalized = true;
                }
                route.AddPos(localizationGraph.ModPos[0], localizationGraph.ChildModBoxes[path[0]].ModIds.First(), onlyOneLocalized);
                return route;
            }

            //Add first mod. If the childBoxes[path[0]].ModIds.Count == 0, means this is an empty childBox. 
            //Otherwise childBoxes[path[0]].ModIds.Count == 1 and childBoxes[path[0]].ModIds only contains one ModId.
            if (localizationGraph.ChildModBoxes[path[0]].ModIds.Count() != 0)
            {                
                route.AddPos(localizationGraph.ModPos[0], localizationGraph.ChildModBoxes[path[0]].ModIds.First(), localizationGraph.array[0][path[0]].CurrentCost > 0);
            }

            for (int i = 1; i < path.Length; i++)
            {
                //If there is a change of the path, get the difference between the two Adjnodes of the array.
                if (path[i] != path[i - 1])
                {
                    var left = GetLeft(localizationGraph.array[i][path[i]].ModBox.ModIds, localizationGraph.array[i - 1][path[i - 1]].ModBox.ModIds).First();

                    var localPeakExist = localizationGraph.array[i - 1][path[i - 1]].CurrentCost > 0 && (localizationGraph.array[i][path[i]].CurrentCost > 0 || i == path.Length -1);
                    route.AddPos(localizationGraph.ModPos[i], left, localPeakExist);
                }
            }

            return route;
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

        //To understand this funciton, ref to "phosphoRS" papar. It is complicated unless you understand how 'phosphoRS' works.  
        //In order to calculate localization probability for each glycosite, one need to get all possible modifications combinations; which is all Routes from a Graph.
        //The function is to get all routes and calculate the 1/P value for each route which is used to calculate localization probability later. 
        public static List<Route> GetAllPaths_CalP(LocalizationGraph localizationGraph, double p, int n)
        {
            List<Route> allPaths = new List<Route>();

            int xlength = localizationGraph.array.Length;
            int ylength = localizationGraph.array.First().Length;

            //temp is a path. check function GetLocalizedPath.
            int[] temp = new int[xlength];
            double[] temp_cost = new double[xlength];

            //A path in graph localization is always the end of the matrix.
            temp[xlength - 1] = ylength - 1;

            PathHelper_CalP(allPaths, localizationGraph, xlength - 1, ylength - 1, temp, temp_cost, p, n);

            return allPaths;
        }
        private static void PathHelper_CalP(List<Route> allPaths, LocalizationGraph localizationGraph, int xind, int yind, int[] temp, double[] temp_costs, double p, int n)
        {
            if (xind == 0)
            {
                var k = temp_costs.Sum() + localizationGraph.NoLocalCost;

                //To understand the math, ref to "phosphoRS" papar.      
                var cp = 1/(1-MathNet.Numerics.Distributions.Binomial.CDF(p, n, k) + MathNet.Numerics.Distributions.Binomial.PMF(p, n, (int)k));               
    
                var route = GetLocalizedPath(localizationGraph, temp);
                route.Score = k;
                route.ReversePScore = cp;
                allPaths.Add(route);
                return;
            }

            foreach (var pre in localizationGraph.array[xind][yind].AllSources)
            {
                xind--;
                yind = pre;
                temp[xind] = yind;
                temp_costs[xind] = localizationGraph.array[xind][yind].CurrentCost;
                PathHelper_CalP(allPaths, localizationGraph, xind, yind, temp, temp_costs, p, n);

                xind++;
            }
        }

        //Dictionary<int, List<Tuple<int, double>>> is <modPos, List<glycanId, site probability>>
        public static Dictionary<int, List<Tuple<int, double>>> CalSiteSpecificLocalizationProbability(List<Route> routes, int[] modPos)
        {
            Dictionary<int, List<Tuple<int, double>>> probabilityMatrix = new Dictionary<int, List<Tuple<int, double>>>();

            Tuple<int, int, double>[][] matrix = new Tuple<int, int, double>[modPos.Length][];

            for (int i = 0; i < modPos.Length; i++)
            {
                matrix[i] = new Tuple<int, int, double>[routes.Count];
                for (int j = 0; j < routes.Count; j++)
                {
                    foreach (var m in routes[j].Mods)
                    {
                        if (m.Item1 == modPos[i])
                        {
                            matrix[i][j] = new Tuple<int, int, double>(m.Item1, m.Item2, routes[j].ReversePScore);
                        }
                    }
                }
            }

            var sum = routes.Sum(p => p.ReversePScore);

            for (int i = 0; i < modPos.Length; i++)
            {
                var gs = matrix[i].Where(p => p!=null).GroupBy(p => p.Item2);

                List<Tuple<int, double>> value = new List<Tuple<int, double>>();

                foreach (var g in gs)
                {
                    var prob = g.Sum(p => p.Item3) / sum;

                    value.Add(new Tuple<int, double>(g.Key, prob));

                }

                probabilityMatrix.Add(modPos[i], value);
            }

            return probabilityMatrix;
        }

    }

}
