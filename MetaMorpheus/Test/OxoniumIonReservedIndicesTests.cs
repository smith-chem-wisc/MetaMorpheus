using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public class OxoniumIonReservedIndicesTests
    {
        /// <summary>
        /// Ordering-freeze guard. Production code (GlycoSearchEngine, GlycoPeptides.DiagonsticFilter)
        /// reaches into Glycan.AllOxoniumIons by the named indices in OxoniumIonReservedIndices. If
        /// anyone reorders or edits AllOxoniumIons, the scaled-int mass at a reserved index changes
        /// and this test fails loudly. Values are the int masses (monoisotopic m/z * 1e5) that must
        /// live at each reserved slot.
        /// </summary>
        [Test]
        public void Reserved_Indices_Match_AllOxoniumIons_Positions()
        {
            Assert.Multiple(() =>
            {
                Assert.That(Glycan.AllOxoniumIons[OxoniumIonReservedIndices.R138], Is.EqualTo(13805550), "R138 (138.055)");
                Assert.That(Glycan.AllOxoniumIons[OxoniumIonReservedIndices.R144], Is.EqualTo(14406607), "R144 (144.066)");
                Assert.That(Glycan.AllOxoniumIons[OxoniumIonReservedIndices.HexNAc204], Is.EqualTo(20408720), "HexNAc204 (204.087)");
                Assert.That(Glycan.AllOxoniumIons[OxoniumIonReservedIndices.NeuAc274], Is.EqualTo(27409268), "NeuAc274 (274.093)");
                Assert.That(Glycan.AllOxoniumIons[OxoniumIonReservedIndices.NeuAc292], Is.EqualTo(29210324), "NeuAc292 (292.103)");
                Assert.That(Glycan.AllOxoniumIons[OxoniumIonReservedIndices.HexHexNAc366], Is.EqualTo(36614002), "HexHexNAc366 (366.140)");
            });
        }
    }
}
