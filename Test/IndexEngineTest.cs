using EngineLayer;
using EngineLayer.Indexing;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using System;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class IndexEngineTest
    {
        #region Public Methods

        [Test]
        public static void TestIndexEngine()
        {
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            
            Protease p = new Protease("Custom Protease2", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            GlobalVariables.ProteaseDictionary.Add(p.Name,p);
            CommonParameters CommonParameters = new CommonParameters(ScoreCutoff: 1, DigestionParams: new DigestionParams(protease: p.Name, MinPeptideLength: 1));
           
            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.BnoB1ions, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(5, results.PeptideIndex.Count);

            var digestedList = proteinList[0].Digest(CommonParameters.DigestionParams, new List<ModificationWithMass>(), variableModifications).ToList();

            Assert.AreEqual(5, digestedList.Count);
            foreach (var fdfd in digestedList)
            {
                Assert.Contains(fdfd.CompactPeptide(TerminusType.None), results.PeptideIndex);
            }
        }

        [Test]
        public static void TestIndexEngineWithWeirdSeq()
        {
            var proteinList = new List<Protein> { new Protein("MQXQ", null) };
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
            {
                modsDictionary.Add(mod, 0);
            }
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }


            Protease protease = new Protease("Custom Protease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            GlobalVariables.ProteaseDictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                DigestionParams: new DigestionParams(
                    protease: protease.Name, 
                    MinPeptideLength: 1, 
                    InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain), 
                ScoreCutoff: 1);
            
            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.BnoB1ions, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(1, results.PeptideIndex.Count);

            Assert.IsNaN(results.PeptideIndex[0].MonoisotopicMassIncludingFixedMods);
            Assert.AreEqual(30000000 + 1, results.FragmentIndex.Length);
        }

        #endregion Public Methods
    }
}