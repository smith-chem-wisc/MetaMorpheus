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
        public static void GraphTest_Graph()
        {          
            int[] modPos = new int[4] { 1,2,3,4 };
            int[] modIds = new int[2] { 8, 9 }; //modIds can be any number, just to index mod in real modification array.
            ModBox modBox = new ModBox(modIds);
            var childBox = ModBox.BuildChildModBoxes(modIds).ToArray();
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
            graph.array[1][1].CummulativeSources = new List<int> { 0 };
            graph.array[2][3].CummulativeSources = new List<int> { 1 };
            graph.array[3][3].CummulativeSources = new List<int> { 3 };
            var path = LocalizationGraph.GetFirstPath(graph.array, childBox);
            Assert.That(path[0] == 0 && path[1] == 1 && path[2] == 3 && path[3] == 3);
        }
    }
}