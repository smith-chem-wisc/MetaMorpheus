using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Fragmentation;
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

        public double NoLocalCost{get; set;}   // Note that we have node for each glycosite, the matched ions before the first node and after the last node is scored here.
        public double TotalScore { get; set; } // Total score is the score of matched ions that are used for localization. For O-glycan, it is the score of all matched c/zDot ions. 

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
        /// <summary>
        /// The function goes through the AdjNode[][] array from left to right, assign weight to each AdjNode, keep track of the heaviest previous AdjNode.
        /// </summary>
        /// <param name="localizationGraph"> The space to store the data </param>
        /// <param name="theScan"> The MS2 scan</param>
        /// <param name="productTolerance"></param>
        /// <param name="products"></param>
        public static void LocalizeOGlycan(LocalizationGraph localizationGraph, Ms2ScanWithSpecificMass theScan, Tolerance productTolerance, List<Product> products)
        {
            var boxSatisfyBox = BoxSatisfyBox(localizationGraph.ChildModBoxes);

            for (int i = 0; i < localizationGraph.ModPos.Length; i++)
            {
                //maxLength: the most mods we can have up to current mod pos; minlengtt: the least mods we can have up to current mod pos.
                int maxLength = i + 1; //For the first node, the maxlength is 1. Means we max have one glycan in this positioin.
                int minlength = localizationGraph.ModBox.ModIds.Length - (localizationGraph.ModPos.Length - 1 - i); //In order to get min number, the min = number of glycan in the box - number of node from the last.
                                                                                                                    // Total 3 glycan in the box, end position is 7, then for position 5, the min = 3 - (7-5) = 1.
                for (int j = 0; j < localizationGraph.ChildModBoxes.Length; j++)
                {
                    if (localizationGraph.ChildModBoxes[j].NumberOfMods <= maxLength && localizationGraph.ChildModBoxes[j].NumberOfMods >= minlength)
                    {
                        AdjNode adjNode = new AdjNode(i, j, localizationGraph.ModPos[i], localizationGraph.ChildModBoxes[j]); //chekc the num of glycan in this node is make sense.

                        double cost = 0;
                        if (i != localizationGraph.ModPos.Length - 1) // check the node is not the last one.
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

                                    var tempCost = cost + localizationGraph.array[i - 1][prej].maxCost; //Try to get the max cost from previous AdjNode.
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

        /// <summary>
        /// Calculate the cost/Score of the Scan.
        /// </summary>
        /// <param name="theScan"></param>
        /// <param name="productTolerance"></param>
        /// <param name="fragments"></param>
        /// <returns> The Score </returns>
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

        /// <summary>
        /// Check does the node1 contain everything in another node2? 
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        /// <returns>Ture, False </returns>
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


        /// <summary>
        /// Build a chart for the node connection rule. Used the chart to check if the next node could be linked to the previous node.
        /// </summary>
        /// <param name="childBoxes"></param>
        /// <returns> Chart (one column is previous, one column is current, the value is boolean)</returns>
        public static Dictionary<int, bool[]> BoxSatisfyBox(ModBox[] childBoxes)
        {
            Dictionary<int, bool[]> boxIdBoxes = new Dictionary<int, bool[]>();
            for (int i = 0; i < childBoxes.Length; i++)
            {
                bool[] idBoxes = new bool[childBoxes.Length];
                for (int j = 0; j <= i; j++)
                {
                    if (childBoxes[i].NumberOfMods <= childBoxes[j].NumberOfMods + 1 && (childBoxes[j].NumberOfMods ==0 || TryGetLeft(childBoxes[i].ModIds, childBoxes[j].ModIds)))
                    { //Check the next node could be the same or one more mod than the previous node. Besdies, the next node should contain all mods that the previous node has.
                        idBoxes[j] = true;
                    }
                }
                boxIdBoxes.Add(i, idBoxes);
            }

            return boxIdBoxes;
        }


        /// <summary>
        /// Try to ll the highest score path in the graph. Start from the last AdjNode[row-1 ][col-1], go back to it Sources, which contains the previous AdjNode with the highest cost.
        /// </summary>
        /// <param name="array"></param>
        /// <param name="boxes"></param>
        /// <returns> The path (one or more) with the higgest Score</returns>
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

        /// <summary>
        /// Get The toppest position path of in the localGraph by recursion Method.
        /// </summary>
        /// <param name="array"></param>
        /// <param name="boxes"></param>
        /// <returns></returns>
        public static int[] GetFirstPath(AdjNode[][] array, ModBox[] boxes)
        {

            int xlength = array.Length;
            int ylength = array.First().Length;

            int[] temp = new int[xlength];

            temp[xlength - 1] = ylength - 1; // That is the last node in the graph, position is last one, and the childBpx is also the last one means the whole glycan.

            FirstPathHelper(array, xlength - 1, ylength - 1, temp);

            return temp;
        }

        private static void FirstPathHelper(AdjNode[][] array, int xind, int yind, int[] temp)
        {
            if (xind == 0) //xind = 0 means, there is just one glycosite. So the node must be the last one in the childBox = whole glycan.
            {
                return; // temp[0] = last one in the childBox = length-1.
            }

            var pre = array[xind][yind].CummulativeSources.First(); // The first one in the CummulativeSources is the toppest previous node.
            xind--;
            yind = pre;
            temp[xind] = yind;
            FirstPathHelper(array, xind, yind, temp);
        }

        /// <summary>
        /// Convert the path inforation into Route object.
        /// </summary>
        /// <param name="localizationGraph"></param>
        /// <param name="path"> ex.[1,1,2,2,5] means the node in the localGraph, first node is ModBox1...last Node is modBox5</param>
        /// <returns> Route object, present in glycosite-glycan pait format </returns>
        public static Route GetLocalizedPath(LocalizationGraph localizationGraph, int[] path)
        {
            Route route = new Route();

            if (path.Length == 1) //If there is only one number in the path, we will assined "the first glycan in the childBox" to the glycosite.
            {
                bool onlyOneLocalized = false;
                if (localizationGraph.TotalScore > 0)
                {
                    onlyOneLocalized = true;
                }
                route.AddPos(localizationGraph.ModPos[0], localizationGraph.ChildModBoxes[path[0]].ModIds.First(), onlyOneLocalized);
                return route;
            }

            //Add first mod in the first glycosite.
            //If the childBoxes[path[0]].ModIds.Count == 0, means this is an empty childBox. 
            //Otherwise childBoxes[path[0]].ModIds.Count == 1 and childBoxes[path[0]].ModIds only contains one ModId.
            if (localizationGraph.ChildModBoxes[path[0]].ModIds.Count() != 0)
            {                
                route.AddPos(localizationGraph.ModPos[0], localizationGraph.ChildModBoxes[path[0]].ModIds.First(), localizationGraph.array[0][path[0]].CurrentCost > 0);
            }

            for (int i = 1; i < path.Length; i++)
            {
                // If there is a change of the path, get the difference between the two Adjnodes of the array.
                // If the node is the same childBox as the previous node. That means there is no modification at this glycosite. We can move on to the next glycosite.
                if (path[i] != path[i - 1])
                {
                    var left = GetLeft(localizationGraph.array[i][path[i]].ModBox.ModIds, localizationGraph.array[i - 1][path[i - 1]].ModBox.ModIds).First();

                    var localPeakExist = localizationGraph.array[i - 1][path[i - 1]].CurrentCost > 0 && (localizationGraph.array[i][path[i]].CurrentCost > 0 || i == path.Length -1);
                    route.AddPos(localizationGraph.ModPos[i], left, localPeakExist);
                }
            }

            return route;
        }


        /// <summary>
        /// Get the difference in glycan between two node.
        /// </summary>
        /// <param name="array1"> The composition in this node. Ex. (0,0,1,2) means the cumulative glycoBox is composed of glycan0 + glycan0 + glycan 1 + glycan 2 </param>
        /// <param name="array2"></param>
        /// <returns> The difference of the glycan composition between the two node.</returns>
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
        /// <summary>
        /// Generate the localization probability chart for each glycosite.
        /// </summary>
        /// <param name="routes"></param>
        /// <param name="modPos"></param>
        /// <returns> A dictionary represent the chart for glycosite Probility. Ex. key = 2 (ModPos), [(0,0.1),(1,0.3),(2,0.6)] means glycan 0 is 10 %, glycan 1 is 30%, glycan 2 is 60% </returns>
        public static Dictionary<int, List<Tuple<int, double>>> CalSiteSpecificLocalizationProbability(List<Route> routes, int[] modPos)
        {
            Dictionary<int, List<Tuple<int, double>>> probabilityMatrix = new Dictionary<int, List<Tuple<int, double>>>();

            Tuple<int, int, double>[][] matrix = new Tuple<int, int, double>[modPos.Length][];

            for (int i = 0; i < modPos.Length; i++) // There are all localization set in the route, we just try to sort the certain glycosite-glycan pairs into the corresponding glycosite.
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
