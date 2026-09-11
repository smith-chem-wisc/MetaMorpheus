using System;
using System.Linq;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class MassDiffAcceptorTypeTests
    {
        [Test]
        public void IsMostAbundant_ReturnsTrue_ForMostAbundantTypes()
        {
            var mostAbundantTypes = new[]
            {
                MassDiffAcceptorType.MostAbundant_Exact,
                MassDiffAcceptorType.MostAbundant_PlusMinusOne,
                MassDiffAcceptorType.MostAbundant_PlusMinusTwo,
            };

            foreach (var type in mostAbundantTypes)
            {
                Assert.That(type.IsMostAbundant(), Is.True, $"{type} should be most abundant");
            }
        }

        [Test]
        public void IsMostAbundant_ReturnsFalse_ForNonMostAbundantTypes()
        {
            var nonMostAbundantTypes = new[]
            {
                MassDiffAcceptorType.Exact,
                MassDiffAcceptorType.OneMM,
                MassDiffAcceptorType.TwoMM,
                MassDiffAcceptorType.ThreeMM,
                MassDiffAcceptorType.ModOpen,
                MassDiffAcceptorType.Open,
                MassDiffAcceptorType.Custom,
                MassDiffAcceptorType.PlusOrMinusThreeMM,
            };

            foreach (var type in nonMostAbundantTypes)
            {
                Assert.That(type.IsMostAbundant(), Is.False, $"{type} should not be most abundant");
            }
        }

        [Test]
        public void IsMostAbundant_CoversAllEnumValues()
        {
            var allTypes = Enum.GetValues<MassDiffAcceptorType>();
            foreach (var type in allTypes)
            {
                var expected = type.ToString().StartsWith("MostAbundant");
                Assert.That(type.IsMostAbundant(), Is.EqualTo(expected), $"Unexpected result for {type}");
            }
        }
    }
}
