using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using MzLibUtil;
using Nett;
using EngineLayer.GlycoSearch;

namespace Test
{
    [TestFixture]
    public class XLTestGraphLocalization
    {
        //The AdjNode, ModBox, Route and Localization Graph are tested here.

        [Test]
        public static void GraphTest_GetLeft()
        {
            int[] array1 = new int[6] { 0, 0, 0, 1, 1, 2 };
            int[] array2 = new int[3] { 0, 0, 1 };
            var left = LocalizationGraph.GetLeft(array1, array2);

            var knowLeft = new int[3] { 0, 1, 2 };
            Assert.That(Enumerable.SequenceEqual(left, knowLeft));
        }

        [Test]
        public static void GraphTest_GetPath()
        {
            int[] modPos = new int[3] { 2, 4, 6 };
            int[] ids = new int[3] { 0, 1, 1 };
            string[] motifs = new string[3] { "X", "X", "X" };
            var modBox = new ModBox(ids, motifs);
            var childBox = ModBox.BuildChildModBoxes(ids, motifs).ToArray();
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modBox, childBox, 0);

            for (int i = 0; i < modPos.Length; i++)
            {
                for (int j = 0; j < childBox.Length; j++)
                {
                    localizationGraph.array[i][j] = new AdjNode(i, j, modPos[i], childBox[j]);
                    localizationGraph.array[i][j].CummulativeSources = new HashSet<int> { j };
                    localizationGraph.array[i][j].CummulativeCost = 1;
                }
            }
            localizationGraph.array[2][5].CummulativeSources = new HashSet<int> { 4, 5 };

            localizationGraph.array[1][4].CummulativeSources = new HashSet<int> { 2, 4 };

            var allPaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, childBox);

            Assert.That(allPaths.Count == 3);
            Assert.That(Enumerable.SequenceEqual(allPaths[0], new int[3] { 2, 4, 5 }));
            Assert.That(Enumerable.SequenceEqual(allPaths[1], new int[3] { 4, 4, 5 }));
            Assert.That(Enumerable.SequenceEqual(allPaths[2], new int[3] { 5, 5, 5 }));

            var firstPath = LocalizationGraph.GetFirstPath(localizationGraph.array, childBox);
            Assert.That(Enumerable.SequenceEqual(firstPath, new int[3] { 2, 4, 5 }));
        }

        [Test]
        public static void GraphTest_GetLocalizedPath()
        {
            int[] modPos = new int[1] { 4 };
            int[] ids = new int[1] { 1 };
            string[] motifs = new string[1] { "X" };
            var modBox = new ModBox(ids, motifs);
            var childBox = ModBox.BuildChildModBoxes(ids, motifs).ToArray();
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modBox, childBox, 0);

            for (int i = 0; i < modPos.Length; i++)
            {
                for (int j = 0; j < childBox.Length; j++)
                {
                    localizationGraph.array[i][j] = new AdjNode(i, j, modPos[i], childBox[j]);
                    localizationGraph.array[i][j].CummulativeSources = new HashSet<int> { j };
                }
            }
            localizationGraph.TotalScore = 1;

            var allPaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, childBox);

            var route = LocalizationGraph.GetLocalizedPath(localizationGraph, allPaths.First());

            Assert.That(route.Mods.First().Item3);
        }

        [Test]
        public static void GraphTest_Graph()
        {          
            int[] modPos = new int[4] { 1,2,3,4 };
            int[] modIds = new int[2] { 8, 9 }; //modIds can be any number, just to index mod in real modification array.
            string[] motifs = new string[2] { "X", "Y" }; 

            ModBox modBox = new ModBox(modIds, motifs);
            var childBox = ModBox.BuildChildModBoxes(modIds, motifs).ToArray();
            Assert.That(childBox.Count() == 4);

            LocalizationGraph graph = new LocalizationGraph(modPos, modBox, childBox, 0);

            for (int i = 0; i < graph.ModPos.Length; i++)
            {

                for (int j = 0; j < graph.ChildModBoxes.Length; j++)
                {

                        AdjNode adjNode = new AdjNode(i, j, graph.ModPos[i], graph.ChildModBoxes[j]);

                        graph.array[i][j] = adjNode;
                    
                }
            }

            //Test AdjNode.
            Assert.That(graph.array[1][2].ModBox.ModIds.First() == 9);
            Assert.That(graph.array[1][3].PosX == 1 && graph.array[1][3].PosY == 3 && graph.array[1][3].ModPos == 2);

            //Test route and GetAnyOnePath.
            var route = LocalizationGraph.GetAnyOnePath(graph);
            Assert.That(route.Mods.First().Item1 == 1 && route.Mods.First().Item2 == 8);
            Assert.That(route.Score == 0);

            //Test Get First path
            //mod on site 2 and site 3.
            graph.array[1][1].CummulativeSources = new HashSet<int> { 0 };
            graph.array[2][3].CummulativeSources = new HashSet<int> { 1 };
            graph.array[3][3].CummulativeSources = new HashSet<int> { 3 };
            var path = LocalizationGraph.GetFirstPath(graph.array, childBox);
            Assert.That(path[0] == 0 && path[1] == 1 && path[2] == 3 && path[3] == 3);
        }

        [Test]
        public static void GraphTest_BoxSatisfyModPos()
        {
            int[] ids = new int[] {2, 2, 3 };          
            string[] motifs = new string[] { "X", "X", "Y" };
            ModBox modBox = new ModBox(ids, motifs);
            ModBox[] childModBox = ModBox.BuildChildModBoxes(ids, motifs).ToArray();

            int[] modPos = new int[4] { 2, 4, 6, 8 };
            string[] boxMotifs = new string[4] { "X", "Y", "X", "X" };
            //LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modMotifs, modBox, childModBox);

            var test_pos0 = LocalizationGraph.BoxSatisfyModPos(boxMotifs, 0, modBox, childModBox[4]);
            Assert.That(test_pos0 == false);

            var test_pos1 = LocalizationGraph.BoxSatisfyModPos(boxMotifs, 1, modBox, childModBox[4]);
            Assert.That(test_pos1 == true);

            var test_pos2 = LocalizationGraph.BoxSatisfyModPos(boxMotifs, 2, modBox, childModBox[4]);
            Assert.That(test_pos2 == true);

            var test_pos3 = LocalizationGraph.BoxSatisfyModPos(boxMotifs, 3, modBox, childModBox[4]);
            Assert.That(test_pos3 == false);
        }



    }
}