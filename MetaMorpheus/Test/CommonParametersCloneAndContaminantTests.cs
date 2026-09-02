using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using MzLibUtil;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    /// <summary>
    /// CloneWithNewValues exists so a front end that edits a subset of the settings does not have to
    /// go through the constructor, which would silently reset every setting it was not passed. That
    /// is the failure it prevents, so it is what these tests pin: an argument that was supplied is
    /// applied, and an argument that was not is carried over untouched.
    /// </summary>
    [TestFixture]
    public static class CommonParametersCloneAndContaminantTests
    {
        /// <summary>
        /// Deliberately non-default on every axis, so "carried over" cannot be confused with
        /// "happened to match the default".
        /// </summary>
        private static CommonParameters DistinctiveParameters() => new CommonParameters(
            dissociationType: DissociationType.ETD,
            doPrecursorDeconvolution: false,
            useProvidedPrecursorInfo: false,
            reportAllAmbiguity: false,
            trimMsMsPeaks: false,
            qValueThreshold: 0.02,
            scoreCutoff: 7,
            precursorMassTolerance: new PpmTolerance(11),
            productMassTolerance: new PpmTolerance(13),
            maxThreadsToUsePerFile: 3,
            digestionParams: new DigestionParams(protease: "Arg-C", maxMissedCleavages: 4),
            listOfModsVariable: new List<(string, string)> { ("Common Variable", "Oxidation on M") },
            listOfModsFixed: new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C") });

        [Test]
        public static void CloneWithNewValuesWithNoArgumentsChangesNothing()
        {
            CommonParameters original = DistinctiveParameters();
            CommonParameters clone = original.CloneWithNewValues();

            Assert.AreNotSame(original, clone);
            Assert.AreEqual(original.DissociationType, clone.DissociationType);
            Assert.AreEqual(original.DoPrecursorDeconvolution, clone.DoPrecursorDeconvolution);
            Assert.AreEqual(original.UseProvidedPrecursorInfo, clone.UseProvidedPrecursorInfo);
            Assert.AreEqual(original.ReportAllAmbiguity, clone.ReportAllAmbiguity);
            Assert.AreEqual(original.TrimMsMsPeaks, clone.TrimMsMsPeaks);
            Assert.AreEqual(original.QValueThreshold, clone.QValueThreshold);
            Assert.AreEqual(original.ScoreCutoff, clone.ScoreCutoff);
            Assert.AreEqual(original.PrecursorMassTolerance.ToString(), clone.PrecursorMassTolerance.ToString());
            Assert.AreEqual(original.ProductMassTolerance.ToString(), clone.ProductMassTolerance.ToString());
            Assert.AreEqual(original.MaxThreadsToUsePerFile, clone.MaxThreadsToUsePerFile);
            Assert.AreEqual(original.DigestionParams, clone.DigestionParams);
            CollectionAssert.AreEqual(original.ListOfModsVariable.ToList(), clone.ListOfModsVariable.ToList());
            CollectionAssert.AreEqual(original.ListOfModsFixed.ToList(), clone.ListOfModsFixed.ToList());
        }

        /// <summary>
        /// One override at a time, each checked against a neighbouring setting that must not have
        /// moved. A test that only asserted the overridden value would still pass if the method
        /// reset everything else, which is the exact bug this method exists to avoid.
        /// </summary>
        [Test]
        public static void CloneWithNewValuesAppliesOnlyTheArgumentSupplied()
        {
            CommonParameters original = DistinctiveParameters();

            CommonParameters dissociation = original.CloneWithNewValues(dissociationType: DissociationType.HCD);
            Assert.AreEqual(DissociationType.HCD, dissociation.DissociationType);
            Assert.AreEqual(original.ScoreCutoff, dissociation.ScoreCutoff);

            CommonParameters deconvolution = original.CloneWithNewValues(doPrecursorDeconvolution: true);
            Assert.IsTrue(deconvolution.DoPrecursorDeconvolution);
            Assert.AreEqual(original.DissociationType, deconvolution.DissociationType);

            CommonParameters providedPrecursor = original.CloneWithNewValues(useProvidedPrecursorInfo: true);
            Assert.IsTrue(providedPrecursor.UseProvidedPrecursorInfo);
            Assert.AreEqual(original.DoPrecursorDeconvolution, providedPrecursor.DoPrecursorDeconvolution);

            CommonParameters ambiguity = original.CloneWithNewValues(reportAllAmbiguity: true);
            Assert.IsTrue(ambiguity.ReportAllAmbiguity);
            Assert.AreEqual(original.TrimMsMsPeaks, ambiguity.TrimMsMsPeaks);

            CommonParameters trim = original.CloneWithNewValues(trimMsMsPeaks: true);
            Assert.IsTrue(trim.TrimMsMsPeaks);
            Assert.AreEqual(original.ReportAllAmbiguity, trim.ReportAllAmbiguity);

            CommonParameters qValue = original.CloneWithNewValues(qValueThreshold: 0.05);
            Assert.AreEqual(0.05, qValue.QValueThreshold);
            Assert.AreEqual(original.ScoreCutoff, qValue.ScoreCutoff);

            CommonParameters score = original.CloneWithNewValues(scoreCutoff: 9);
            Assert.AreEqual(9, score.ScoreCutoff);
            Assert.AreEqual(original.QValueThreshold, score.QValueThreshold);

            CommonParameters precursorTolerance = original.CloneWithNewValues(precursorMassTolerance: new PpmTolerance(25));
            Assert.AreEqual(25, precursorTolerance.PrecursorMassTolerance.Value);
            Assert.AreEqual(original.ProductMassTolerance.Value, precursorTolerance.ProductMassTolerance.Value);

            CommonParameters productTolerance = original.CloneWithNewValues(productMassTolerance: new PpmTolerance(35));
            Assert.AreEqual(35, productTolerance.ProductMassTolerance.Value);
            Assert.AreEqual(original.PrecursorMassTolerance.Value, productTolerance.PrecursorMassTolerance.Value);

            CommonParameters threads = original.CloneWithNewValues(maxThreadsToUsePerFile: 1);
            Assert.AreEqual(1, threads.MaxThreadsToUsePerFile);
            Assert.AreEqual(original.ScoreCutoff, threads.ScoreCutoff);

            IDigestionParams trypsin = new DigestionParams(protease: "trypsin", maxMissedCleavages: 1);
            CommonParameters digestion = original.CloneWithNewValues(digestionParams: trypsin);
            Assert.AreEqual(trypsin, digestion.DigestionParams);
            Assert.AreEqual(original.DissociationType, digestion.DissociationType);

            var newVariable = new List<(string, string)> { ("Common Variable", "Acetylation on K") };
            CommonParameters variableMods = original.CloneWithNewValues(listOfModsVariable: newVariable);
            CollectionAssert.AreEqual(newVariable, variableMods.ListOfModsVariable.ToList());
            CollectionAssert.AreEqual(original.ListOfModsFixed.ToList(), variableMods.ListOfModsFixed.ToList());

            var newFixed = new List<(string, string)> { ("Common Fixed", "Oxidation on M") };
            CommonParameters fixedMods = original.CloneWithNewValues(listOfModsFixed: newFixed);
            CollectionAssert.AreEqual(newFixed, fixedMods.ListOfModsFixed.ToList());
            CollectionAssert.AreEqual(original.ListOfModsVariable.ToList(), fixedMods.ListOfModsVariable.ToList());
        }

        [Test]
        public static void CloneWithNewValuesDoesNotMutateTheOriginal()
        {
            CommonParameters original = DistinctiveParameters();
            DissociationType before = original.DissociationType;

            original.CloneWithNewValues(dissociationType: DissociationType.CID, scoreCutoff: 99);

            Assert.AreEqual(before, original.DissociationType);
            Assert.AreEqual(7, original.ScoreCutoff);
        }

        /// <summary>
        /// CloneWithNewDissociationType is now a one-argument call to CloneWithNewValues. Pinning the
        /// equivalence means the delegation cannot be broken without a test saying so - it has three
        /// call sites in the search tasks.
        /// </summary>
        [Test]
        public static void CloneWithNewDissociationTypeMatchesCloneWithNewValues()
        {
            CommonParameters original = DistinctiveParameters();

            CommonParameters viaNarrowHelper = original.CloneWithNewDissociationType(DissociationType.HCD);
            CommonParameters viaGeneralHelper = original.CloneWithNewValues(dissociationType: DissociationType.HCD);

            Assert.AreEqual(DissociationType.HCD, viaNarrowHelper.DissociationType);
            Assert.AreEqual(viaGeneralHelper.DissociationType, viaNarrowHelper.DissociationType);
            Assert.AreEqual(viaGeneralHelper.ScoreCutoff, viaNarrowHelper.ScoreCutoff);
            Assert.AreEqual(viaGeneralHelper.QValueThreshold, viaNarrowHelper.QValueThreshold);
            Assert.AreEqual(viaGeneralHelper.MaxThreadsToUsePerFile, viaNarrowHelper.MaxThreadsToUsePerFile);
            Assert.AreEqual(viaGeneralHelper.DigestionParams, viaNarrowHelper.DigestionParams);
            Assert.AreEqual(viaGeneralHelper.PrecursorMassTolerance.Value, viaNarrowHelper.PrecursorMassTolerance.Value);
        }

        /// <summary>
        /// The contaminant filename heuristic, which CMD, the WPF GUI and the Avalonia GUI all now
        /// share. It is deliberately ordinal so a machine's locale cannot change which databases a
        /// search treats as contaminants.
        /// </summary>
        [Test]
        [TestCase("contaminants.fasta", true)]
        [TestCase("CONTAMINANTS.fasta", true)]
        [TestCase("Contaminant_list.xml", true)]
        [TestCase("crap.fasta", true)]
        [TestCase("CRAP.fasta", true)]
        [TestCase("cRaP.xml", true)]
        [TestCase(@"C:\databases\contaminants\uniprot.xml", true)]
        [TestCase("uniprot-human.xml", false)]
        [TestCase("smalldb.fasta", false)]
        [TestCase("", false)]
        public static void LooksLikeContaminantMatchesNameCaseInsensitively(string filePath, bool expected)
        {
            Assert.AreEqual(expected, DbForTask.LooksLikeContaminant(filePath));
        }

        /// <summary>
        /// A guard on the reason the heuristic is ordinal rather than culture-sensitive: "contaminant"
        /// contains a dotted i, and ToUpper() under a Turkish locale maps it to a different character.
        /// </summary>
        [Test]
        public static void LooksLikeContaminantIsNotCultureSensitive()
        {
            var originalCulture = System.Threading.Thread.CurrentThread.CurrentCulture;
            try
            {
                System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("tr-TR");
                Assert.IsTrue(DbForTask.LooksLikeContaminant("CONTAMINANTS.fasta"));
                Assert.IsTrue(DbForTask.LooksLikeContaminant("contaminants.fasta"));
            }
            finally
            {
                System.Threading.Thread.CurrentThread.CurrentCulture = originalCulture;
            }
        }

        /// <summary>
        /// The constructor path CMD uses, so the heuristic and the flag it feeds stay wired together.
        /// </summary>
        [Test]
        public static void DbForTaskCarriesTheContaminantFlagItWasGiven()
        {
            var contaminant = new DbForTask("contaminants.fasta", DbForTask.LooksLikeContaminant("contaminants.fasta"));
            var target = new DbForTask("uniprot-human.xml", DbForTask.LooksLikeContaminant("uniprot-human.xml"));

            Assert.IsTrue(contaminant.IsContaminant);
            Assert.IsFalse(target.IsContaminant);
        }
    }
}
