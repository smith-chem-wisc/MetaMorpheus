using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using Transcriptomics.Digestion;

namespace Test.GuiTests
{
    [TestFixture]
    public static class LoadDigestionAgentsTest
    {
        private static string GetHeaderLine(string[] tsvLines)
        {
            // The TSV files may have comment lines starting with '#' before the header
            return tsvLines.First(line => !line.StartsWith("#") && line.TrimStart().StartsWith("Name\t"));
        }

        private static string GetHeaderLine(Stream tsvStream)
        {
            var reader = new StreamReader(tsvStream);

            string fileContent = reader.ReadToEnd();
            string[] lines = fileContent.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            return GetHeaderLine(lines);
        }

        /// <summary>
        /// Builds a TSV data row that matches the column layout in the given header.
        /// Populates Name, Motif/Sequences Inducing Cleavage, and Specificity/Cleavage Specificity columns.
        /// </summary>
        private static string BuildDataRow(string header, string name, string motif, string specificity)
        {
            var columns = header.Split('\t');
            var values = new string[columns.Length];
            for (int i = 0; i < columns.Length; i++)
            {
                string col = columns[i].Trim().ToLowerInvariant();
                if (col == "name")
                    values[i] = name;
                else if (col == "motif" || col == "sequences inducing cleavage")
                    values[i] = motif;
                else if (col == "specificity" || col == "cleavage specificity")
                    values[i] = specificity;
                else
                    values[i] = "";
            }
            return string.Join("\t", values);
        }

        /// <summary>
        /// Verifies that custom proteases defined in proteases_custom.tsv are merged
        /// into the main ProteaseDictionary when LoadDigestionAgents runs.
        /// </summary>
        [Test]
        public static void TestCustomProteasesAreMerged()
        {
            // Arrange: write a custom protease to the custom file
            string customPath = GlobalVariables.CustomProteasePath;

            // Read the header from the standard proteases.tsv
            var assembly = typeof(ProteaseDictionary).Assembly;
            string EmbeddedProteaseResourceName = "Proteomics.ProteolyticDigestion.proteases.tsv";

            var stream = assembly.GetManifestResourceStream(EmbeddedProteaseResourceName);
            string header = GetHeaderLine(stream);

            // Back up any existing custom file
            string backupPath = customPath + ".bak";
            bool hadExistingCustomFile = File.Exists(customPath);
            if (hadExistingCustomFile)
                File.Copy(customPath, backupPath, overwrite: true);

            try
            {
                // Write a custom protease entry
                var lines = new List<string>
                {
                    header,
                    BuildDataRow(header, "test-custom-protease", "K|", "full")
                };
                File.WriteAllLines(customPath, lines);

                // Act: reload
                GlobalVariables.SetUpGlobalVariables();

                // Assert: the custom protease should be present
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("test-custom-protease"), Is.True,
                    "Custom protease should be merged into the dictionary");
            }
            finally
            {
                // Restore original custom file
                if (hadExistingCustomFile)
                    File.Copy(backupPath, customPath, overwrite: true);
                else if (File.Exists(customPath))
                    File.Delete(customPath);

                if (File.Exists(backupPath))
                    File.Delete(backupPath);

                // Reset state
                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// Verifies that when a custom protease has the same name as a default protease,
        /// the default entry wins (is not overwritten by the custom one).
        /// </summary>
        [Test]
        public static void TestDefaultProteaseWinsOnNameCollision()
        {
            // Arrange
            string customPath = GlobalVariables.CustomProteasePath;

            // Read the header from the standard proteases.tsv
            var assembly = typeof(ProteaseDictionary).Assembly;
            string EmbeddedProteaseResourceName = "Proteomics.ProteolyticDigestion.proteases.tsv";

            var stream = assembly.GetManifestResourceStream(EmbeddedProteaseResourceName);
            string header = GetHeaderLine(stream);

            string backupPath = customPath + ".bak";
            bool hadExistingCustomFile = File.Exists(customPath);
            if (hadExistingCustomFile)
                File.Copy(customPath, backupPath, overwrite: true);

            try
            {
                // Write a custom entry with the same name as a default protease ("trypsin")
                // but with a different cleavage motif (Asp-N style "|D")
                var lines = new List<string>
                {
                    header,
                    BuildDataRow(header, "trypsin", "|D", "full")
                };
                File.WriteAllLines(customPath, lines);

                // Act
                GlobalVariables.SetUpGlobalVariables();

                // Assert: trypsin should still exist and use the default motif (K|,R|), not |D
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("trypsin"), Is.True);
                var trypsin = ProteaseDictionary.Dictionary["trypsin"];

                // The default trypsin cleaves at K and R; verify it wasn't replaced by the custom |D
                Assert.That(trypsin.DigestionMotifs.Count, Is.GreaterThanOrEqualTo(2),
                    "Default trypsin should have at least 2 motifs (K| and R|). If only 1, the custom entry overwrote the default.");
            }
            finally
            {
                if (hadExistingCustomFile)
                    File.Copy(backupPath, customPath, overwrite: true);
                else if (File.Exists(customPath))
                    File.Delete(customPath);

                if (File.Exists(backupPath))
                    File.Delete(backupPath);

                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// Verifies that custom RNases defined in rnases_custom.tsv are merged
        /// into the main RnaseDictionary when LoadDigestionAgents runs.
        /// </summary>
        [Test]
        public static void TestCustomRnasesAreMerged()
        {
            // Arrange
            string customPath = GlobalVariables.CustomRnasePath;

            // Read the header from the standard proteases.tsv
            var assembly = typeof(RnaseDictionary).Assembly;
            string EmbeddedProteaseResourceName = "Transcriptomics.Digestion.rnases.tsv";

            var stream = assembly.GetManifestResourceStream(EmbeddedProteaseResourceName);
            string header = GetHeaderLine(stream);

            string backupPath = customPath + ".bak";
            bool hadExistingCustomFile = File.Exists(customPath);
            if (hadExistingCustomFile)
                File.Copy(customPath, backupPath, overwrite: true);

            try
            {
                var lines = new List<string>
                {
                    header,
                    BuildDataRow(header, "test-custom-rnase", "A|", "full")
                };
                File.WriteAllLines(customPath, lines);

                // Act
                GlobalVariables.SetUpGlobalVariables();

                // Assert
                Assert.That(RnaseDictionary.Dictionary.ContainsKey("test-custom-rnase"), Is.True,
                    "Custom RNase should be merged into the dictionary");
            }
            finally
            {
                if (hadExistingCustomFile)
                    File.Copy(backupPath, customPath, overwrite: true);
                else if (File.Exists(customPath))
                    File.Delete(customPath);

                if (File.Exists(backupPath))
                    File.Delete(backupPath);

                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// Verifies that when a custom RNase has the same name as a default RNase,
        /// the default entry wins (is not overwritten by the custom one).
        /// </summary>
        [Test]
        public static void TestDefaultRnaseWinsOnNameCollision()
        {
            // Arrange
            string customPath = GlobalVariables.CustomRnasePath;

            // Read the header from the standard proteases.tsv
            var assembly = typeof(RnaseDictionary).Assembly;
            string EmbeddedProteaseResourceName = "Transcriptomics.Digestion.rnases.tsv";

            var stream = assembly.GetManifestResourceStream(EmbeddedProteaseResourceName);
            string header = GetHeaderLine(stream);

            string backupPath = customPath + ".bak";
            bool hadExistingCustomFile = File.Exists(customPath);
            if (hadExistingCustomFile)
                File.Copy(customPath, backupPath, overwrite: true);

            try
            {
                // Write a custom entry colliding with "RNase T1" (which normally cleaves at G|)
                // but with a different motif (C|)
                var lines = new List<string>
                {
                    header,
                    BuildDataRow(header, "RNase T1", "C|", "full")
                };
                File.WriteAllLines(customPath, lines);

                // Act
                GlobalVariables.SetUpGlobalVariables();

                // Assert: RNase T1 should still exist and use the default motif (G|), not C|
                Assert.That(RnaseDictionary.Dictionary.ContainsKey("RNase T1"), Is.True);
                var rnaseT1 = RnaseDictionary.Dictionary["RNase T1"];

                // Default RNase T1 cleaves at G; if the custom C| replaced it, there would be
                // a single motif matching C instead of G
                Assert.That(rnaseT1.DigestionMotifs.Any(m => m.InducingCleavage == "G"), Is.True,
                    "Default RNase T1 should retain its G cleavage motif, not be overwritten by the custom entry");
            }
            finally
            {
                if (hadExistingCustomFile)
                    File.Copy(backupPath, customPath, overwrite: true);
                else if (File.Exists(customPath))
                    File.Delete(customPath);

                if (File.Exists(backupPath))
                    File.Delete(backupPath);

                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// When proteases_custom.tsv does not exist, LoadDigestionAgents creates it
        /// with a header line. Verify this behavior.
        /// </summary>
        [Test]
        public static void TestCustomProteaseFileCreatedWhenMissing()
        {
            string customPath = GlobalVariables.CustomProteasePath;

            string backupPath = customPath + ".bak";
            bool hadExistingCustomFile = File.Exists(customPath);
            if (hadExistingCustomFile)
                File.Copy(customPath, backupPath, overwrite: true);

            try
            {
                // Remove the custom file
                if (File.Exists(customPath))
                    File.Delete(customPath);

                // Act
                GlobalVariables.SetUpGlobalVariables();

                // Assert: the custom file should now exist with at least a header
                Assert.That(File.Exists(customPath), Is.True,
                    "proteases_custom.tsv should be created when it doesn't exist");
                var lines = File.ReadAllLines(customPath);
                Assert.That(lines.Length, Is.GreaterThanOrEqualTo(1),
                    "Custom file should contain at least a header line");
                Assert.That(lines.Any(l => l.TrimStart().StartsWith("Name\t")), Is.True,
                    "Custom file should contain the TSV header");
            }
            finally
            {
                if (hadExistingCustomFile)
                    File.Copy(backupPath, customPath, overwrite: true);

                if (File.Exists(backupPath))
                    File.Delete(backupPath);

                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// When rnases_custom.tsv does not exist, LoadDigestionAgents creates it
        /// with a header line. Verify this behavior.
        /// </summary>
        [Test]
        public static void TestCustomRnaseFileCreatedWhenMissing()
        {
            string customPath = GlobalVariables.CustomRnasePath;

            string backupPath = customPath + ".bak";
            bool hadExistingCustomFile = File.Exists(customPath);
            if (hadExistingCustomFile)
                File.Copy(customPath, backupPath, overwrite: true);

            try
            {
                if (File.Exists(customPath))
                    File.Delete(customPath);

                // Act
                GlobalVariables.SetUpGlobalVariables();

                // Assert
                Assert.That(File.Exists(customPath), Is.True,
                    "rnases_custom.tsv should be created when it doesn't exist");
                var lines = File.ReadAllLines(customPath);
                Assert.That(lines.Length, Is.GreaterThanOrEqualTo(1),
                    "Custom file should contain at least a header line");
                Assert.That(lines.Any(l => l.TrimStart().StartsWith("Name\t")), Is.True,
                    "Custom file should contain the TSV header");
            }
            finally
            {
                if (hadExistingCustomFile)
                    File.Copy(backupPath, customPath, overwrite: true);

                if (File.Exists(backupPath))
                    File.Delete(backupPath);

                GlobalVariables.SetUpGlobalVariables();
            }
        }

      
    }
}
