using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.ExceptionServices;
using Transcriptomics.Digestion;

namespace Test.GuiTests
{
    [TestFixture]
    [NonParallelizable]
    public static class LoadDigestionAgentsTest
    {
        private static void InvokeLoadDigestionAgents()
        {
            var method = typeof(GlobalVariables).GetMethod("LoadDigestionAgents", BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(method, Is.Not.Null, "Could not find GlobalVariables.LoadDigestionAgents via reflection.");

            try
            {
                method.Invoke(null, null);
            }
            catch (TargetInvocationException ex) when (ex.InnerException != null)
            {
                ExceptionDispatchInfo.Capture(ex.InnerException).Throw();
            }
        }

        private static void BackupPath(string path, out string backupPath, out bool hadExistingFile, out bool hadExistingDirectory)
        {
            backupPath = path + ".bak";
            hadExistingFile = File.Exists(path);
            hadExistingDirectory = Directory.Exists(path);

            if (File.Exists(backupPath))
                File.Delete(backupPath);
            else if (Directory.Exists(backupPath))
                Directory.Delete(backupPath, true);

            if (hadExistingFile)
                File.Copy(path, backupPath, overwrite: true);
            else if (hadExistingDirectory)
                Directory.Move(path, backupPath);
        }

        private static void RestorePath(string path, string backupPath, bool hadExistingFile, bool hadExistingDirectory)
        {
            if (File.Exists(path))
                File.Delete(path);
            else if (Directory.Exists(path))
                Directory.Delete(path, true);

            if (hadExistingFile)
                File.Copy(backupPath, path, overwrite: true);
            else if (hadExistingDirectory)
                Directory.Move(backupPath, path);

            if (File.Exists(backupPath))
                File.Delete(backupPath);
            else if (Directory.Exists(backupPath))
                Directory.Delete(backupPath, true);
        }

        private static void DeletePathIfExists(string path)
        {
            if (File.Exists(path))
                File.Delete(path);
            else if (Directory.Exists(path))
                Directory.Delete(path, true);
        }

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

            BackupPath(customPath, out string backupPath, out bool hadExistingCustomFile, out bool hadExistingCustomDirectory);

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
                RestorePath(customPath, backupPath, hadExistingCustomFile, hadExistingCustomDirectory);

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

            BackupPath(customPath, out string backupPath, out bool hadExistingCustomFile, out bool hadExistingCustomDirectory);

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
                RestorePath(customPath, backupPath, hadExistingCustomFile, hadExistingCustomDirectory);

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

            BackupPath(customPath, out string backupPath, out bool hadExistingCustomFile, out bool hadExistingCustomDirectory);

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
                RestorePath(customPath, backupPath, hadExistingCustomFile, hadExistingCustomDirectory);

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

            BackupPath(customPath, out string backupPath, out bool hadExistingCustomFile, out bool hadExistingCustomDirectory);

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
                RestorePath(customPath, backupPath, hadExistingCustomFile, hadExistingCustomDirectory);

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

            BackupPath(customPath, out string backupPath, out bool hadExistingCustomFile, out bool hadExistingCustomDirectory);

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
                RestorePath(customPath, backupPath, hadExistingCustomFile, hadExistingCustomDirectory);

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

            BackupPath(customPath, out string backupPath, out bool hadExistingCustomFile, out bool hadExistingCustomDirectory);

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
                RestorePath(customPath, backupPath, hadExistingCustomFile, hadExistingCustomDirectory);

                GlobalVariables.SetUpGlobalVariables();
            }
        }

        [Test]
        public static void TestCustomProteaseLoadCatchIsHitWhenCustomFileIsLocked()
        {
            string customPath = GlobalVariables.CustomProteasePath;
            BackupPath(customPath, out string backupPath, out bool hadExistingFile, out bool hadExistingDirectory);

            try
            {
                using var lockedFile = new FileStream(customPath, FileMode.Create, FileAccess.ReadWrite, FileShare.None);

                var ex = Assert.Throws<MetaMorpheusException>(() => InvokeLoadDigestionAgents());

                Assert.That(ex!.Message, Does.Contain("Error loading custom proteases with error message:"));
                Assert.That(ex.InnerException, Is.Not.Null);
            }
            finally
            {
                RestorePath(customPath, backupPath, hadExistingFile, hadExistingDirectory);
                GlobalVariables.SetUpGlobalVariables();
            }
        }

        [Test]
        public static void TestCustomProteaseCreationCatchIsHitWhenCustomPathIsDirectory()
        {
            string customPath = GlobalVariables.CustomProteasePath;
            BackupPath(customPath, out string backupPath, out bool hadExistingFile, out bool hadExistingDirectory);

            try
            {
                DeletePathIfExists(customPath);
                Directory.CreateDirectory(customPath);

                var ex = Assert.Throws<MetaMorpheusException>(() => InvokeLoadDigestionAgents());

                Assert.That(ex!.Message, Does.Contain("Error creating default custom protease file with error message:"));
                Assert.That(ex.InnerException, Is.Not.Null);
            }
            finally
            {
                RestorePath(customPath, backupPath, hadExistingFile, hadExistingDirectory);
                GlobalVariables.SetUpGlobalVariables();
            }
        }

        [Test]
        public static void TestCustomRnaseLoadCatchIsHitWhenCustomFileIsLocked()
        {
            string customPath = GlobalVariables.CustomRnasePath;
            BackupPath(customPath, out string backupPath, out bool hadExistingFile, out bool hadExistingDirectory);

            try
            {
                using var lockedFile = new FileStream(customPath, FileMode.Create, FileAccess.ReadWrite, FileShare.None);

                var ex = Assert.Throws<MetaMorpheusException>(() => InvokeLoadDigestionAgents());

                Assert.That(ex!.Message, Does.Contain("Error loading custom rnases with error message:"));
                Assert.That(ex.InnerException, Is.Not.Null);
            }
            finally
            {
                RestorePath(customPath, backupPath, hadExistingFile, hadExistingDirectory);
                GlobalVariables.SetUpGlobalVariables();
            }
        }

        [Test]
        public static void TestCustomRnaseCreationCatchIsHitWhenCustomPathIsDirectory()
        {
            string customPath = GlobalVariables.CustomRnasePath;
            BackupPath(customPath, out string backupPath, out bool hadExistingFile, out bool hadExistingDirectory);

            try
            {
                DeletePathIfExists(customPath);
                Directory.CreateDirectory(customPath);

                var ex = Assert.Throws<MetaMorpheusException>(() => InvokeLoadDigestionAgents());

                Assert.That(ex!.Message, Does.Contain("Error creating default custom rnase file with error message:"));
                Assert.That(ex.InnerException, Is.Not.Null);
            }
            finally
            {
                RestorePath(customPath, backupPath, hadExistingFile, hadExistingDirectory);
                GlobalVariables.SetUpGlobalVariables();
            }
        }

      
    }
}
