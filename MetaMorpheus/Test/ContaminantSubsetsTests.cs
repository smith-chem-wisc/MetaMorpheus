using EngineLayer;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// The shipped contaminant panel carries the 48 proteins of the Sigma UPS1/UPS2 spike-in standard (the "Sigma UPS" group of
    /// the GPM cRAP list). They are ordinary human proteins (cytochrome c, tau, creatine kinase M, PRDX1 ...), so in non-human
    /// samples they claim native peptides as contamination. Contaminants/Subsets splits the panel into the 216 entries that are
    /// not UPS and the 48 that are. These tests keep the split exact if the full panel is ever edited.
    /// </summary>
    [TestFixture]
    public static class ContaminantSubsetsTests
    {
        // cRAP "Sigma UPS" group, https://www.thegpm.org/crap/
        private static readonly HashSet<string> UpsAccessions = new()
        {
            "P02768", "P01008", "P08758", "P61769", "P55957", "P00915", "P00918", "P04040", "P07339", "P08311", "P01031", "P02741",
            "P00167", "P99999", "P01133", "P05413", "P06396", "P08263", "P09211", "P69905", "P68871", "P01344", "P10145", "P06732",
            "P00709", "P41159", "P61626", "P02144", "Q15843", "P15559", "P16083", "P01127", "P62937", "Q06830", "P01112", "P02753",
            "P00441", "P63165", "P12081", "P10636", "P10599", "P01375", "P02787", "P02788", "P51965", "O00762", "P63279", "P62979",
        };

        private static string ContaminantsDir => Path.Combine(GlobalVariables.DataDir, "Contaminants");

        private static Dictionary<string, string> Load(string path)
        {
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(path, true, DecoyType.None, GlobalVariables.AllModsKnown, true,
                new List<string>(), out Dictionary<string, Modification> _);
            // the loader also expands annotated sequence variants (P00441_A5S ...); compare the canonical entries
            return proteins.Where(p => !p.AppliedSequenceVariations.Any()).ToDictionary(p => p.Accession, p => p.BaseSequence);
        }

        [Test]
        public static void SubsetsPartitionTheFullPanelOnTheUpsAccessions()
        {
            var full = Load(Path.Combine(ContaminantsDir, "MetaMorpheusContaminants.xml"));
            var noUps = Load(Path.Combine(ContaminantsDir, "Subsets", "MetaMorpheusContaminants_NoUPS.xml"));
            var ups = Load(Path.Combine(ContaminantsDir, "Subsets", "MetaMorpheusContaminants_UPS.xml"));

            Assert.That(UpsAccessions, Has.Count.EqualTo(48));
            Assert.That(ups.Keys, Is.EquivalentTo(UpsAccessions));
            Assert.That(noUps.Keys, Is.EquivalentTo(full.Keys.Except(UpsAccessions)));
            Assert.That(noUps.Count + ups.Count, Is.EqualTo(full.Count));

            foreach (var (accession, sequence) in noUps.Concat(ups))
            {
                Assert.That(sequence, Is.EqualTo(full[accession]), accession);
            }
        }

        [Test]
        public static void AddDefaultContaminantsStillFindsOnlyTheFullPanel()
        {
            // The GUI's +ADD DEFAULT CONTAMINANTS adds every file directly in Contaminants, so the subsets must stay in a sub folder
            string[] topLevel = Directory.GetFiles(ContaminantsDir);
            Assert.That(topLevel.Select(Path.GetFileName), Is.EquivalentTo(new[] { "MetaMorpheusContaminants.xml" }));
        }
    }
}
