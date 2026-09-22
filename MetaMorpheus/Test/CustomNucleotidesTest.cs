using EngineLayer;
using NUnit.Framework;
using System.IO;
using Transcriptomics;

namespace Test
{
    [TestFixture]
    public static class CustomNucleotidesTest
    {
        [Test]
        public static void TestCustomNucleotideReading()
        {
            string nucleotidePath = Path.Combine(GlobalVariables.DataDir, "CustomNucleotides", "CustomNucleotides.txt");
            string[] originalLines = File.ReadAllLines(nucleotidePath);
            string[] updatedLines = new string[originalLines.Length + 1];
            originalLines.CopyTo(updatedLines, 0);
            updatedLines[^1] = "TestNucleotide\tq\tTnt\tC5H5N2O2";
            File.WriteAllLines(nucleotidePath, updatedLines);

            try
            {
                GlobalVariables.LoadCustomNucleotides();

                Assert.That(Nucleotide.TryGetResidue('q', out Nucleotide nucleotide), Is.True);
                Assert.That(nucleotide.Name, Is.EqualTo("TestNucleotide"));
                Assert.That(nucleotide.Symbol, Is.EqualTo("Tnt"));
                Assert.That(nucleotide.BaseChemicalFormula.Formula, Is.EqualTo("C5H5N2O2"));
            }
            finally
            {
                File.WriteAllLines(nucleotidePath, originalLines);
            }
        }
    }
}
