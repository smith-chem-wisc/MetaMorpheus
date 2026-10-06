using EngineLayer;
using EngineLayer.DatabaseLoading;
using Nett;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// FlashLFQ's Bayesian protein fold-change step, run from a search task. Four copies of one small spectra file
    /// make a two-condition design (control and treated, two bioreps each), so every protein the search quantifies
    /// has peptides in both conditions.
    /// </summary>
    [TestFixture]
    public static class BayesianProteinQuantTests
    {
        private const string BayesianFileName = "BayesianFoldChangeAnalysis.tsv";

        private static string SpectraSource => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
        private static string Database => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");

        /// <summary>
        /// Copies the spectra four times into <paramref name="folder"/> and, unless <paramref name="designRows"/> is
        /// null, writes an experimental design beside them. Returns the spectra paths.
        /// </summary>
        private static List<string> MakeFourFileExperiment(string folder, string[] designRows)
        {
            Directory.CreateDirectory(folder);
            var files = new List<string>();
            foreach (string name in new[] { "ctrl1", "ctrl2", "trt1", "trt2" })
            {
                string path = Path.Combine(folder, name + ".mzml");
                File.Copy(SpectraSource, path, overwrite: true);
                files.Add(path);
            }

            if (designRows != null)
            {
                File.WriteAllLines(Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName),
                    new[] { "FileName\tCondition\tBiorep\tFraction\tTechrep" }.Concat(designRows));
            }

            return files;
        }

        private static readonly string[] TwoConditions =
        {
            "ctrl1.mzml\tcontrol\t1\t1\t1",
            "ctrl2.mzml\tcontrol\t2\t1\t1",
            "trt1.mzml\ttreated\t1\t1\t1",
            "trt2.mzml\ttreated\t2\t1\t1",
        };

        private static SearchTask BayesianSearchTask(string controlCondition) => new SearchTask
        {
            SearchParameters = new SearchParameters
            {
                Normalize = true,
                DoBayesianProteinQuant = true,
                BayesianControlCondition = controlCondition,
            },
            CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 1),
        };

        /// <summary>Runs the task and returns every warning it raised.</summary>
        private static List<string> RunCollectingWarnings(SearchTask task, string outputFolder, List<string> spectra)
        {
            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                Directory.CreateDirectory(outputFolder);
                task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(Database, false) }, spectra, "bayesian");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }
            return warnings;
        }

        [Test]
        public static void ByDefaultTheBayesianStepIsOffWithSeed42AndCutoff01()
        {
            var parameters = new SearchParameters();
            Assert.That(parameters.DoBayesianProteinQuant, Is.False);
            Assert.That(parameters.BayesianControlCondition, Is.Null);
            Assert.That(parameters.BayesianFoldChangeCutoff, Is.EqualTo(0.1));
            Assert.That(parameters.BayesianRandomSeed, Is.EqualTo(42));

            // A task file written before these settings existed loads with the same defaults.
            string oldToml = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml");
            Assert.That(File.ReadAllText(oldToml), Does.Not.Contain("Bayesian"));
            var loaded = Toml.ReadFile<SearchTask>(oldToml, MetaMorpheusTask.tomlConfig).SearchParameters;
            Assert.That(loaded.DoBayesianProteinQuant, Is.False);
            Assert.That(loaded.BayesianFoldChangeCutoff, Is.EqualTo(0.1));
            Assert.That(loaded.BayesianRandomSeed, Is.EqualTo(42));
        }

        [Test]
        public static void TheSettingsSurviveATomlRoundTrip()
        {
            var task = BayesianSearchTask("control");
            task.SearchParameters.BayesianFoldChangeCutoff = 0.25;
            task.SearchParameters.BayesianRandomSeed = 7;
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "BayesianRoundTrip.toml");
            try
            {
                Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);
                var read = Toml.ReadFile<SearchTask>(path, MetaMorpheusTask.tomlConfig).SearchParameters;

                Assert.That(read.DoBayesianProteinQuant, Is.True);
                Assert.That(read.BayesianControlCondition, Is.EqualTo("control"));
                Assert.That(read.BayesianFoldChangeCutoff, Is.EqualTo(0.25));
                Assert.That(read.BayesianRandomSeed, Is.EqualTo(7));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The step runs, writes one row per protein group that AllQuantifiedProteinGroups.tsv shows and no other,
        /// records its settings and seed in results.txt, and gives the same fold changes on a second run.
        /// </summary>
        [Test]
        public static void ATwoConditionDesignWritesReproducibleFoldChangesForTheProteinTablesGroups()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "BayesianTwoConditions");
            try
            {
                var spectra = MakeFourFileExperiment(Path.Combine(folder, "spectra"), TwoConditions);
                string firstRun = Path.Combine(folder, "run1");
                string secondRun = Path.Combine(folder, "run2");

                var warnings = RunCollectingWarnings(BayesianSearchTask("control"), firstRun, spectra);
                RunCollectingWarnings(BayesianSearchTask("control"), secondRun, spectra);

                Assert.That(warnings, Has.None.Contains("Bayesian"));

                string[] bayesian = File.ReadAllLines(Path.Combine(firstRun, BayesianFileName));
                Assert.That(bayesian[0], Does.StartWith("Protein Group\t"));
                var bayesianGroups = bayesian.Skip(1).Where(line => line.Length > 0).Select(line => line.Split('\t')[0]).ToList();
                Assert.That(bayesianGroups, Is.Not.Empty);
                Assert.That(bayesian.Skip(1).Where(line => line.Length > 0).Select(line => line.Split('\t')[4]), Is.All.EqualTo("treated"));

                var proteinTableGroups = File.ReadAllLines(Path.Combine(firstRun, "AllQuantifiedProteinGroups.tsv"))
                    .Skip(1).Select(line => line.Split('\t')[0]).ToHashSet();
                Assert.That(bayesianGroups, Is.SubsetOf(proteinTableGroups));
                Assert.That(bayesianGroups, Is.Unique);

                Assert.That(File.ReadAllText(Path.Combine(firstRun, "results.txt")),
                    Does.Contain("Bayesian protein quantification: control condition control, fold-change cutoff 0.1, random seed 42"));

                Assert.That(File.ReadAllText(Path.Combine(secondRun, BayesianFileName)),
                    Is.EqualTo(File.ReadAllText(Path.Combine(firstRun, BayesianFileName))), "the fixed seed makes the step reproducible");
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A group the protein table hides is not written to the Bayesian table either. FlashLFQ's own groups keep
        /// it: one protein of the database is loaded as a contaminant, and contaminants are not written.
        /// </summary>
        [Test]
        public static void AGroupTheProteinTableHidesIsNotInTheBayesianTable()
        {
            const string contaminantAccession = "P10662";
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "BayesianHiddenGroup");
            try
            {
                var spectra = MakeFourFileExperiment(Path.Combine(folder, "spectra"), TwoConditions);
                var entries = File.ReadAllText(Database).Split('>', StringSplitOptions.RemoveEmptyEntries);
                string targets = Path.Combine(folder, "targets.fasta");
                string contaminants = Path.Combine(folder, "contaminants.fasta");
                File.WriteAllText(targets, string.Concat(entries.Where(e => !e.Contains("|" + contaminantAccession + "|")).Select(e => ">" + e)));
                File.WriteAllText(contaminants, string.Concat(entries.Where(e => e.Contains("|" + contaminantAccession + "|")).Select(e => ">" + e)));

                var task = BayesianSearchTask("control");
                task.SearchParameters.WriteContaminants = false;
                string output = Path.Combine(folder, "out");
                Directory.CreateDirectory(output);
                task.RunTask(output, new List<DbForTask> { new DbForTask(targets, false), new DbForTask(contaminants, true) }, spectra, "bayesian");

                var proteinTableGroups = File.ReadAllLines(Path.Combine(output, "AllQuantifiedProteinGroups.tsv"))
                    .Skip(1).Select(line => line.Split('\t')[0]).ToList();
                var bayesianGroups = File.ReadAllLines(Path.Combine(output, BayesianFileName))
                    .Skip(1).Where(line => line.Length > 0).Select(line => line.Split('\t')[0]).ToList();

                Assert.That(proteinTableGroups, Has.None.Contains(contaminantAccession), "premise: the protein table hides the contaminant");
                Assert.That(bayesianGroups, Is.Not.Empty);
                Assert.That(bayesianGroups, Has.None.Contains(contaminantAccession));
                Assert.That(bayesianGroups, Is.EquivalentTo(proteinTableGroups));
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// Every design the step cannot run on warns, skips the Bayesian file, and leaves the rest of
        /// quantification alone: FlashLFQ would throw on the first, and a throw loses every label-free table.
        /// </summary>
        [TestCase("absent", "is not one of the experimental design's conditions (control, treated)")]
        [TestCase("oneCondition", "has only one (control)")]
        [TestCase("noDesign", "no ExperimentalDesign.tsv was found")]
        public static void ADesignTheStepCannotUseWarnsAndQuantificationContinues(string scenario, string expectedWarning)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "BayesianRefused_" + scenario);
            try
            {
                string[] design = scenario switch
                {
                    "oneCondition" => new[]
                    {
                        "ctrl1.mzml\tcontrol\t1\t1\t1",
                        "ctrl2.mzml\tcontrol\t2\t1\t1",
                        "trt1.mzml\tcontrol\t3\t1\t1",
                        "trt2.mzml\tcontrol\t4\t1\t1",
                    },
                    "noDesign" => null,
                    _ => TwoConditions,
                };

                var task = BayesianSearchTask(scenario == "absent" ? "absent" : "control");
                task.SearchParameters.Normalize = scenario != "noDesign"; // normalization would skip quant with no design
                var spectra = MakeFourFileExperiment(Path.Combine(folder, "spectra"), design);
                string output = Path.Combine(folder, "out");

                var warnings = RunCollectingWarnings(task, output, spectra);

                Assert.That(warnings, Has.Some.Contains(expectedWarning));
                Assert.That(warnings, Has.Some.Contains("Skipping Bayesian protein quantification; the rest of quantification continues"));
                Assert.That(File.Exists(Path.Combine(output, BayesianFileName)), Is.False);
                Assert.That(File.ReadAllLines(Path.Combine(output, "AllQuantifiedPeptides.tsv")).Length, Is.GreaterThan(1));
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// FlashLFQ quantifies a SILAC search before its files are split into light and heavy samples, and the
        /// split does not rerun the Bayesian step, so a SILAC search warns and writes no Bayesian table.
        /// </summary>
        [Test]
        public static void ASilacSearchWarnsAndWritesNoBayesianTable()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "BayesianSilac");
            try
            {
                Residue heavyLysine = new("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All);
                Residue lightLysine = Residue.GetResidue('K');
                var task = BayesianSearchTask("control");
                task.SearchParameters.SilacLabels = new List<SilacLabel>
                {
                    new SilacLabel(lightLysine.Letter, heavyLysine.Letter, heavyLysine.ThisChemicalFormula.Formula, heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass)
                };
                var spectra = MakeFourFileExperiment(Path.Combine(folder, "spectra"), TwoConditions);
                string output = Path.Combine(folder, "out");

                var warnings = RunCollectingWarnings(task, output, spectra);

                Assert.That(warnings, Has.Some.Contains("Bayesian protein quantification does not support SILAC"));
                Assert.That(File.Exists(Path.Combine(output, BayesianFileName)), Is.False);
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        [Test]
        public static void RunningWithoutNormalizationWarnsThatLoadingDifferencesBecomeFoldChanges()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "BayesianUnnormalized");
            try
            {
                var task = BayesianSearchTask("control");
                task.SearchParameters.Normalize = false;
                var spectra = MakeFourFileExperiment(Path.Combine(folder, "spectra"), TwoConditions);
                string output = Path.Combine(folder, "out");

                var warnings = RunCollectingWarnings(task, output, spectra);

                Assert.That(warnings, Has.Some.Contains("running on unnormalized intensities"));
                Assert.That(File.Exists(Path.Combine(output, BayesianFileName)), Is.True);
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }
    }
}
