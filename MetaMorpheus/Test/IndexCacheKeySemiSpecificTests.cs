using EngineLayer;
using EngineLayer.Indexing;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// A cached index built for a semi-specific search before mzLib #1303 must never be reused.
    /// </summary>
    /// <remarks>
    /// <para>MetaMorpheus reuses a saved index when <see cref="IndexingEngine.ToString"/> of the new search equals the text
    /// saved with the index. That text names the search mode and terminus but no mzLib version. Before mzLib #1303,
    /// SearchModeType Semi with FragmentationTerminus Both digested into C-terminal seeds; since then it digests into the
    /// semi-specific peptides. So an index of seeds cached by an older MetaMorpheus had exactly the key a new semi-specific
    /// Modern or Glyco search computes, and would be reused, silently searching the seeds again.</para>
    /// <para>The key for Semi + Both therefore carries one more line, which no older index has. Every other combination's key
    /// is unchanged, so their cached indexes stay valid.</para>
    /// </remarks>
    [TestFixture]
    public static class IndexCacheKeySemiSpecificTests
    {
        private const string SemiSpecificPeptidesLine = "semiSpecificDigestion: peptides, not seeds (mzLib #1303)";

        private static string CacheKey(CleavageSpecificity searchModeType, FragmentationTerminus terminus)
        {
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus));
            var engine = new IndexingEngine(new List<IBioPolymer> { new Protein("PEPTIDEKPEPTIDER", "P1") }, new List<Modification>(), new List<Modification>(),
                null, null, null, 1, DecoyType.None, commonParameters, null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            return engine.ToString();
        }

        [Test]
        public static void CacheKey_ForSemiWithBothTermini_NamesThatItIndexesPeptides()
        {
            Assert.That(CacheKey(CleavageSpecificity.Semi, FragmentationTerminus.Both), Does.Contain(SemiSpecificPeptidesLine));
        }

        [Test]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.Both)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.N)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.C)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.N)]
        public static void CacheKey_ForEveryOtherCombination_IsUnchanged(CleavageSpecificity searchModeType, FragmentationTerminus terminus)
        {
            Assert.That(CacheKey(searchModeType, terminus), Does.Not.Contain("semiSpecificDigestion"),
                "digestion for these did not change in mzLib #1303, so their cached indexes must stay reusable");
        }
    }
}
