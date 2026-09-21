using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// mzIdentML must say whether a search was semi-specific: <c>Enzyme/@semiSpecific</c>.
    /// </summary>
    /// <remarks>
    /// <para>The writer computed the value from the protease's own specificity, which is Full for trypsin even in a
    /// semi-specific search (the semi-specificity comes from <c>SearchModeType</c>), and never set
    /// <c>semiSpecificSpecified</c>, without which the XML serializer leaves the attribute out entirely. So no mzIdentML file
    /// MetaMorpheus wrote ever carried it, and a reader assumed a fully specific search. Now that "semi-trypsin" is gone,
    /// every semi-specific search is trypsin + SearchModeType Semi, so this matters for all of them.</para>
    /// <para>mzIdentML defines semiSpecific as exactly one terminus following the enzyme rules, so a non-specific search
    /// (SearchModeType None) is not semi-specific.</para>
    /// </remarks>
    [TestFixture]
    public static class MzIdentMlSemiSpecificTests
    {
        private static string EnzymeSemiSpecificAttribute(string path)
        {
            using var reader = XmlReader.Create(path);
            while (reader.Read())
            {
                if (reader.NodeType == XmlNodeType.Element && reader.LocalName == "Enzyme")
                    return reader.GetAttribute("semiSpecific");
            }
            Assert.Fail("the file has no Enzyme element");
            return null;
        }

        private static string WriteAndReadSemiSpecific(DigestionAgent protease, bool? semiSpecific, string testName)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(MzIdentMlSemiSpecificTests)}_{testName}.mzID");
            try
            {
                MzIdentMLWriter.WriteMzIdentMl(new List<SpectralMatch>(), new List<ProteinGroup>(), new List<Modification>(), new List<Modification>(),
                    null, new List<DigestionAgent> { protease }, new PpmTolerance(20), new PpmTolerance(20), 2, path, true, semiSpecific);
                return EnzymeSemiSpecificAttribute(path);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        /// <summary>
        /// The search's own specificity, when the caller gives it, decides. A semi-specific search writes
        /// <c>semiSpecific="true"</c>; any other search leaves the attribute out, exactly as every file did before.
        /// </summary>
        /// <remarks>
        /// Review of #2812 (nbollis): writing <c>semiSpecific="false"</c> for every fully specific and non-specific search
        /// changed those files for readers that tell a missing attribute from an explicit false, for no gain. The schema
        /// makes the attribute optional, so only a semi-specific search, the case that was wrong, now writes it.
        /// </remarks>
        [Test]
        [TestCase(true, "true")]
        [TestCase(false, null)]
        public static void WriteMzIdentMl_WithTheSearchSpecificity_WritesSemiSpecificOnlyWhenTrue(bool semiSpecific, string expected)
        {
            Assert.That(WriteAndReadSemiSpecific(ProteaseDictionary.Dictionary["trypsin"], semiSpecific, $"Given{semiSpecific}"), Is.EqualTo(expected));
        }

        /// <summary>
        /// Without it, the protease's own specificity decides, and again only a semi-specific protease writes the attribute.
        /// </summary>
        [Test]
        public static void WriteMzIdentMl_WithoutTheSearchSpecificity_WritesTheProteaseSpecificity()
        {
            Assert.That(WriteAndReadSemiSpecific(ProteaseDictionary.Dictionary["trypsin"], null, "FullProtease"), Is.Null,
                "a fully specific protease leaves the attribute out, as before #2812");

            var semiProtease = new Protease("MzIdentMlSemiSpecificTests-semi", CleavageSpecificity.Semi, null, null, DigestionMotif.ParseDigestionMotifsFromString("K|,R|"));
            Assert.That(WriteAndReadSemiSpecific(semiProtease, null, "SemiProtease"), Is.EqualTo("true"));
        }

        /// <summary>What a search passes for its mzIdentML files, for every search mode and a protease that is itself semi.</summary>
        [Test]
        [TestCase("trypsin", CleavageSpecificity.Full, FragmentationTerminus.Both, false)]
        [TestCase("trypsin", CleavageSpecificity.Semi, FragmentationTerminus.Both, true)]
        [TestCase("trypsin", CleavageSpecificity.Semi, FragmentationTerminus.N, true)]
        [TestCase("trypsin", CleavageSpecificity.None, FragmentationTerminus.N, false)]
        public static void IsSemiSpecificForMzIdentMl_IsTrueForASemiSpecificSearch(string protease, CleavageSpecificity searchModeType, FragmentationTerminus terminus, bool expected)
        {
            var digestionParams = new DigestionParams(protease, searchModeType: searchModeType, fragmentationTerminus: terminus);
            Assert.That(PostSearchAnalysisTask.IsSemiSpecificForMzIdentMl(digestionParams), Is.EqualTo(expected));
        }
    }
}
