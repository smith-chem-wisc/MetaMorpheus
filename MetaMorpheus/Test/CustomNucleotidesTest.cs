using EngineLayer;
using NUnit.Framework;
using System;
using System.IO;
using System.Linq;
using Transcriptomics;

namespace Test
{
    [TestFixture]
    public static class CustomNucleotidesTest
    {
        [Test]
        public static void TestCustomNucleotideReading()
        {
            RunWithTemporaryNucleotideFile(nucleotidePath =>
            {
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
            });
        }

        [Test]
        public static void LoadCustomNucleotidesCreatesFileWhenMissing()
        {
            RunWithTemporaryNucleotideFile(nucleotidePath =>
            {
                if (File.Exists(nucleotidePath))
                    File.Delete(nucleotidePath);

                GlobalVariables.LoadCustomNucleotides();

                Assert.That(File.Exists(nucleotidePath), Is.True);
                string[] lines = File.ReadAllLines(nucleotidePath);
                Assert.That(lines.Length, Is.EqualTo(1));
                Assert.That(lines[0], Is.EqualTo("Name\tOneLetterAbbr.\tSymbol\tBaseChemicalFormula"));
            });
        }

        [Test]
        public static void LoadCustomNucleotidesAllowsIdenticalExistingResidue()
        {
            RunWithTemporaryNucleotideFile(nucleotidePath =>
            {
                Assert.That(Nucleotide.TryGetResidue('A', out Nucleotide adenine), Is.True,
                    "Precondition: 'A' nucleotide must exist.");

                string[] lines =
                {
                    "Name\tOneLetterAbbr.\tSymbol\tBaseChemicalFormula",
                    $"{adenine.Name}\t{adenine.Letter}\t{adenine.Symbol}\t{adenine.BaseChemicalFormula.Formula}"
                };
                File.WriteAllLines(nucleotidePath, lines);

                Assert.DoesNotThrow(GlobalVariables.LoadCustomNucleotides);
            });
        }

        [Test]
        public static void LoadCustomNucleotidesThrowsForConflictingDuplicateAssignment()
        {
            RunWithTemporaryNucleotideFile(nucleotidePath =>
            {
                Assert.That(Nucleotide.TryGetResidue('A', out Nucleotide adenine), Is.True,
                    "Precondition: 'A' nucleotide must exist.");

                char unusedLetter = GetUnusedNucleotideLetter();
                string[] lines =
                {
                    "Name\tOneLetterAbbr.\tSymbol\tBaseChemicalFormula",
                    $"{adenine.Name}\t{unusedLetter}\tX{unusedLetter}x\t{adenine.BaseChemicalFormula.Formula}"
                };
                File.WriteAllLines(nucleotidePath, lines);

                MetaMorpheusException exception = Assert.Throws<MetaMorpheusException>(GlobalVariables.LoadCustomNucleotides);
                Assert.That(exception.Message, Does.Contain("The nucleotide name, letter, or symbol is already assigned to a different nucleotide."));
                Assert.That(exception.Message, Does.Contain("Line 2"));
            });
        }

        [Test]
        public static void LoadCustomNucleotidesThrowsForInvalidFormula()
        {
            RunWithTemporaryNucleotideFile(nucleotidePath =>
            {
                string[] lines =
                {
                    "Name\tOneLetterAbbr.\tSymbol\tBaseChemicalFormula",
                    "BadFormula\tq\tBdf\tnot_a_formula"
                };
                File.WriteAllLines(nucleotidePath, lines);

                MetaMorpheusException exception = Assert.Throws<MetaMorpheusException>(GlobalVariables.LoadCustomNucleotides);
                Assert.That(exception.Message, Does.StartWith("Error while reading 'CustomNucleotides.txt'. Line 2 was not in the correct format:"));
            });
        }

        [Test]
        public static void LoadCustomNucleotidesSkipsBlankAndShortLines()
        {
            RunWithTemporaryNucleotideFile(nucleotidePath =>
            {
                char unusedLetter = GetUnusedNucleotideLetter();
                string customName = "EdgeCaseNucleotide_" + Guid.NewGuid().ToString("N")[..8];
                string customSymbol = "S" + Guid.NewGuid().ToString("N")[..2];

                string[] lines =
                {
                    "Name\tOneLetterAbbr.\tSymbol\tBaseChemicalFormula",
                    "",
                    "TooShort\tz",
                    $"{customName}\t{unusedLetter}\t{customSymbol}\tC5H5N2O2"
                };
                File.WriteAllLines(nucleotidePath, lines);

                Assert.DoesNotThrow(GlobalVariables.LoadCustomNucleotides);
                Assert.That(Nucleotide.TryGetResidue(unusedLetter, out Nucleotide loaded), Is.True);
                Assert.That(loaded.Name, Is.EqualTo(customName));
                Assert.That(loaded.Symbol, Is.EqualTo(customSymbol));
            });
        }

        private static void RunWithTemporaryNucleotideFile(Action<string> action)
        {
            string nucleotidePath = Path.Combine(GlobalVariables.DataDir, "CustomNucleotides", "CustomNucleotides.txt");
            bool hadFile = File.Exists(nucleotidePath);
            string[] originalLines = hadFile ? File.ReadAllLines(nucleotidePath) : Array.Empty<string>();

            try
            {
                action(nucleotidePath);
            }
            finally
            {
                if (hadFile)
                    File.WriteAllLines(nucleotidePath, originalLines);
                else if (File.Exists(nucleotidePath))
                    File.Delete(nucleotidePath);
            }
        }

        private static char GetUnusedNucleotideLetter() => Enumerable.Range('a', 'z' - 'a' + 1)
            .Select(value => (char)value)
            .First(letter => !Nucleotide.TryGetResidue(letter, out _));
    }
}
