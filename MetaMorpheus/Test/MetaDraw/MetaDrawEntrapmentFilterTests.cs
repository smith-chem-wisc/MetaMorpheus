using System;
using System.IO;
using System.Linq;
using GuiFunctions;
using NUnit.Framework;
using Readers;

namespace Test.MetaDraw
{
    /// <summary>
    /// The psmtsv writer marks entrapment ET / ED. MetaDraw's filter compared the label to "T", "D"
    /// and "C" exactly, so every entrapment target PSM would have vanished from the display, and an
    /// entrapment decoy stayed hidden even with decoys shown.
    /// </summary>
    [TestFixture]
    public class MetaDrawEntrapmentFilterTests
    {
        [TearDown]
        public void Reset() => MetaDrawSettings.ResetSettings();

        [TestCase("T", false, false, true)]
        [TestCase("ET", false, false, true)]
        [TestCase("T|ET", false, false, true)]
        [TestCase("D", false, false, false)]
        [TestCase("ED", false, false, false)]
        [TestCase("D", true, false, true)]
        [TestCase("ED", true, false, true)]
        [TestCase("ED|D", true, false, true)]
        [TestCase("C", false, true, true)]
        [TestCase("C", false, false, false)]
        // Mixed labels were hidden before entrapment existed, and still are.
        [TestCase("T|D", true, false, false)]
        [TestCase("T|C", false, true, false)]
        [TestCase("ET|D", true, false, false)]
        public void EntrapmentIsFilteredAsTheTargetOrDecoyItIs(string label, bool showDecoys, bool showContaminants, bool accepted)
        {
            MetaDrawSettings.ResetSettings();
            MetaDrawSettings.ShowDecoys = showDecoys;
            MetaDrawSettings.ShowContaminants = showContaminants;

            SpectrumMatchFromTsv psm = PsmLabelled(label);

            Assert.That(MetaDrawSettings.FilterAcceptsPsm(psm), Is.EqualTo(accepted));
        }

        /// <summary>A real PSM that passes every other filter, with its label rewritten.</summary>
        private static SpectrumMatchFromTsv PsmLabelled(string label)
        {
            string source = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\noSAreaderTest.psmtsv");
            string[] lines = File.ReadAllLines(source);
            string[] header = lines[0].Split('\t');
            int labelColumn = Array.IndexOf(header, SpectrumMatchFromTsvHeader.DecoyContaminantTarget);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"entrapmentFilter_{Guid.NewGuid():N}.psmtsv");
            File.WriteAllLines(path, lines.Select((line, i) =>
            {
                if (i == 0) return line;
                string[] fields = line.Split('\t');
                fields[labelColumn] = label;
                return string.Join('\t', fields);
            }));

            try
            {
                var psm = SpectrumMatchTsvReader.ReadTsv(path, out _)
                    .First(p => p.QValue <= MetaDrawSettings.QValueFilter
                                && (p.QValueNotch == null || p.QValueNotch <= MetaDrawSettings.QValueFilter));
                Assert.That(psm.DecoyContamTarget, Is.EqualTo(label));
                return psm;
            }
            finally
            {
                File.Delete(path);
            }
        }
    }
}
