using EngineLayer;
using NUnit.Framework;
using TaskLayer;
using EngineLayer.Indexing;
using System.Collections.Generic;
using Proteomics;
using MassSpectrometry;
using UsefulProteomicsDatabases;
using Proteomics.ProteolyticDigestion;
using System.IO;

namespace Test
{
    [TestFixture]
    public static class SearchTaskTest
    {
        [Test]
        public static void MassDiffAceptorTest()
        {
            SearchTask searchTask = new SearchTask();
            var result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, searchTask.SearchParameters.MassDiffAcceptorType, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("1mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.TwoMM, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("2mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.ThreeMM, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("3mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.ModOpen, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("-187andUp"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Open, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("OpenSearch"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "custom ppmAroundZero 4");
            Assert.That(result.FileNameAddition.Equals("4ppmAroundZero"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Exact, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("5ppmAroundZero"));
        }

        [Test]
        public static void ParseSearchModeTest()
        {
            SearchTask searchTask = new SearchTask();
            var result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom dot 5 ppm 0,1.0029,2.0052");
            Assert.That(result.NumNotches == 3);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom dot 5 da 0,1.0029,2.0052");
            Assert.That(result.NumNotches == 3);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom interval [0;5],[0;5]");
            Assert.That(result.NumNotches == 1);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom OpenSearch 5");
            Assert.That(result.FileNameAddition.Equals("OpenSearch"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom daltonsAroundZero 5");
            Assert.That(result.FileNameAddition.Equals("5daltonsAroundZero"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom ppmAroundZero 5");
            Assert.That(result.FileNameAddition.Equals("5ppmAroundZero"));

            try
            {
                result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom Test 5");
            }
            catch(MetaMorpheusException)
            {
                return;
            }
            Assert.Fail();
        }

        [Test]
        public static void SemiSpecificTest()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    SearchType = SearchType.NonSpecific
                }
            };
        }          
    }
}
 