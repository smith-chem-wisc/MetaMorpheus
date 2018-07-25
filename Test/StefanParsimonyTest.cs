using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class StefanParsimonyTest
    {
        
        [Test]
        public static void ParsimonyVariableTreatAsUnique()
        {
            bool modPeptidesAreUnique = true;

            // set up mods
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            var mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);

            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(new List<ProductType> { ProductType.B, ProductType.Y });

            Protease protease = new Protease("kprotease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            // modified version of protein
            var protein1 = new Protein("PEPTIDEM", "accession1");
            // unmodified version of protein
            var protein2 = new Protein("YYYKPEPTIDEM", "accession2");

            var pep1 = protein1.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<ModificationWithMass> { mod }, new List<ModificationWithMass>()).First();
            var pep2 = protein2.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<ModificationWithMass> { mod }, new List<ModificationWithMass>()).ToList()[1];

            // check to make sure mod is present
            Assert.That(pep1.Sequence.Equals(pep2.Sequence));
            Assert.That(pep1.NumMods == 1);
            Assert.That(pep2.NumMods == 1);

            // build the dictionary for input to parsimony
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            var cp1 = pep1.CompactPeptide(terminusType);
            var cp2 = pep2.CompactPeptide(terminusType);
            Assert.That(cp1.Equals(cp2));

            compactPeptideToProteinPeptideMatching.Add(pep1.CompactPeptide(terminusType), new HashSet<PeptideWithSetModifications> { pep1 });
            Assert.That(compactPeptideToProteinPeptideMatching.ContainsKey(cp2));
            compactPeptideToProteinPeptideMatching[cp2].Add(pep2);

            // apply parsimony
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique,new CommonParameters(), new List<string>());
            pae.Run();

            // check to make sure both peptides are associated with both proteins
            Assert.That(compactPeptideToProteinPeptideMatching.Count == 1);
            Assert.That(compactPeptideToProteinPeptideMatching.First().Value.Count == 2);
            var seq = compactPeptideToProteinPeptideMatching.First().Value.First().Sequence;
            foreach (var sequence in compactPeptideToProteinPeptideMatching.First().Value)
            {
                Assert.That(sequence.Sequence.Equals(seq));
            }
        }

        [Test]
        public static void ParsimonyVariableDontTreatAsUnique()
        {
            bool modPeptidesAreUnique = false;

            // set up mods
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            var mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            Protease protease = new Protease("k Protease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(new List<ProductType> { ProductType.B, ProductType.Y });

            // modified version of protein
            var protein1 = new Protein("PEPTIDEM", "accession1");
            // unmodified version of protein
            var protein2 = new Protein("YYYKPEPTIDEM", "accession2");

            var pep1 = protein1.Digest(new DigestionParams(protease: "k Protease", minPeptideLength: 1), new List<ModificationWithMass> { mod }, new List<ModificationWithMass>()).First();
            var pep2 = protein2.Digest(new DigestionParams(protease: "k Protease", minPeptideLength: 1), new List<ModificationWithMass> { mod }, new List<ModificationWithMass>()).ToList()[1];

            // check to make sure mod is present
            Assert.That(pep1.Sequence.Equals(pep2.Sequence));
            Assert.That(pep1.NumMods == 1);
            Assert.That(pep2.NumMods == 1);

            // build the dictionary for input to parsimony
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            var cp1 = pep1.CompactPeptide(terminusType);
            var cp2 = pep2.CompactPeptide(terminusType);
            Assert.That(cp1.Equals(cp2));

            compactPeptideToProteinPeptideMatching.Add(pep1.CompactPeptide(terminusType), new HashSet<PeptideWithSetModifications> { pep1 });
            Assert.That(compactPeptideToProteinPeptideMatching.ContainsKey(cp2));
            compactPeptideToProteinPeptideMatching[cp2].Add(pep2);

            // apply parsimony
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique, new CommonParameters(), new List<string>());
            pae.Run();


            // check to make sure both peptides are associated with both proteins
            Assert.That(compactPeptideToProteinPeptideMatching.Count == 1);
            Assert.That(compactPeptideToProteinPeptideMatching.First().Value.Count == 2);
            var seq = compactPeptideToProteinPeptideMatching.First().Value.First().Sequence;
            foreach (var sequence in compactPeptideToProteinPeptideMatching.First().Value)
            {
                Assert.That(sequence.Sequence.Equals(seq));
            }
        }

        [Test]
        public static void ParsimonyLocalizeableTreatAsUnique()
        {
            bool modPeptidesAreUnique = true;

            // set up mods
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            var mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);

            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(new List<ProductType> { ProductType.B, ProductType.Y });
            Protease protease = new Protease("k protease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            // modified version of protein
            var protein1 = new Protein("PEPTIDEM", "accession1");
            // unmodified version of protein
            var protein2 = new Protein("YYYKPEPTIDEM", "accession2");

            var pep1 = protein1.Digest(new DigestionParams(protease.Name, minPeptideLength: 1), new List<ModificationWithMass> { mod }, new List<ModificationWithMass>()).First();
            var pep2 = protein2.Digest(new DigestionParams(protease.Name, minPeptideLength: 1), new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList()[1];

            // check to make sure mod is present
            Assert.That(pep1.Sequence != pep2.Sequence);
            Assert.That(pep1.NumMods == 1);
            Assert.That(pep2.NumMods == 0);

            // build the dictionary for input to parsimony
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            compactPeptideToProteinPeptideMatching.Add(pep1.CompactPeptide(terminusType), new HashSet<PeptideWithSetModifications> { pep1 });
            compactPeptideToProteinPeptideMatching.Add(pep2.CompactPeptide(terminusType), new HashSet<PeptideWithSetModifications> { pep2 });

            // apply parsimony
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique, new CommonParameters(), new List<string>());
            pae.Run();

            // check to make sure both peptides are NOT associated with both proteins
            Assert.That(compactPeptideToProteinPeptideMatching.Count == 2);
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                Assert.That(kvp.Value.Count == 1);
            }
        }

        [Test]
        public static void ParsimonyLocalizeableDontTreatAsUnique()
        {
            bool modPeptidesAreUnique = false;

            // set up mods
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            var mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);

            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(new List<ProductType> { ProductType.B, ProductType.Y });
            Protease protease = new Protease("kProtease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            // modified version of protein
            var protein1 = new Protein("PEPTIDEM", "accession1");
            // unmodified version of protein
            var protein2 = new Protein("YYYKPEPTIDEM", "accession2");

            var pep1 = protein1.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<ModificationWithMass> { mod }, new List<ModificationWithMass>()).First();
            var pep2 = protein2.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList()[1];

            // check to make sure mod is present
            Assert.That(pep1.Sequence != pep2.Sequence);
            Assert.That(pep1.NumMods == 1);
            Assert.That(pep2.NumMods == 0);

            // build the dictionary for input to parsimony
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            compactPeptideToProteinPeptideMatching.Add(pep1.CompactPeptide(terminusType), new HashSet<PeptideWithSetModifications> { pep1 });
            compactPeptideToProteinPeptideMatching.Add(pep2.CompactPeptide(terminusType), new HashSet<PeptideWithSetModifications> { pep2 });

            // apply parsimony
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, modPeptidesAreUnique, new CommonParameters(), new List<string>());
            pae.Run();

            // check to make sure both peptides are associated with both proteins
            Assert.That(compactPeptideToProteinPeptideMatching.Count == 2);
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                Assert.That(kvp.Value.Count == 2);
                var seq = kvp.Value.First().Sequence;

                foreach (var peptide in kvp.Value)
                {
                    Assert.That(peptide.Sequence.Equals(seq));
                }
            }
            var test1 = compactPeptideToProteinPeptideMatching.First().Value.ToList();
            Assert.That(test1[0].OneBasedStartResidueInProtein != test1[1].OneBasedStartResidueInProtein);
            Assert.That(test1[0].OneBasedEndResidueInProtein != test1[1].OneBasedEndResidueInProtein);
        }

        [Test]
        public static void ParsimonyWeirdCatch()
        {
            Protein protein1 = new Protein("MATSIK", "protein1", isDecoy: true);
            Protein protein2 = new Protein("MATSLK", "protein2");
            Protein protein3 = new Protein("MTASIK", "protein3");

            IEnumerable<ModificationWithMass> allKnownFixedModifications = new List<ModificationWithMass>();
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var pep1 = protein1.Digest(digestionParams, allKnownFixedModifications, variableModifications).First();
            var pep2 = protein2.Digest(digestionParams, allKnownFixedModifications, variableModifications).First();
            var pep3 = protein3.Digest(digestionParams, allKnownFixedModifications, variableModifications).First();

            CompactPeptide compactPeptide1 = pep1.CompactPeptide(TerminusType.None);
            CompactPeptide compactPeptide2 = pep2.CompactPeptide(TerminusType.None);
            CompactPeptide compactPeptide3 = pep3.CompactPeptide(TerminusType.None);

            Assert.AreEqual(compactPeptide1, compactPeptide2);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {compactPeptide1, new HashSet<PeptideWithSetModifications>{pep1, pep2} },
                {compactPeptide3, new HashSet<PeptideWithSetModifications>{pep3} }
            };

            var cool = (ProteinParsimonyResults)new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, false, new CommonParameters(), new List<string>()).Run();

            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching.Count);

            // Only 1 because the target is removed!!!
            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide1].Count);
            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide2].Count);

            Assert.AreEqual(1, compactPeptideToProteinPeptideMatching[compactPeptide3].Count);
        }

        
        private static Tuple<List<PeptideSpectralMatch>, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>, MassDiffAcceptor, bool, CompactPeptideBase, CompactPeptideBase> GetInfo(bool localizeable)
        {
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: 2, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain, maxModsForPeptides: 1), scoreCutoff: 1);


            // Alanine = Glycine + CH2
            Protein protein1 = new Protein("MA", "protein1");
            Protein protein2 = new Protein("MG", "protein2");
            Protein protein3;
            double monoisotopicMass = Chemistry.ChemicalFormula.ParseFormula("CH2").MonoisotopicMass;
            ModificationMotif.TryGetMotif("G", out ModificationMotif motif1);
            ModificationMotif.TryGetMotif("A", out ModificationMotif motif2);
            TerminusLocalization modificationSites = TerminusLocalization.Any;
            List<ModificationWithMass> allKnownFixedModifications = new List<ModificationWithMass>
            {
                new ModificationWithMass("CH2 on Glycine", null, motif1, modificationSites, monoisotopicMass)
            };
            List<ModificationWithMass> variableModifications;

            ModificationWithMass alanineMod = new ModificationWithMass("CH2 on Alanine", null, motif2, modificationSites, monoisotopicMass);

            if (localizeable)
            {
                variableModifications = new List<ModificationWithMass>();
                IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{alanineMod} }
                };
                protein3 = new Protein("MA", "protein3", oneBasedModifications: oneBasedModifications);
            }
            else
            {
                variableModifications = new List<ModificationWithMass>();
                variableModifications = new List<ModificationWithMass> { alanineMod };
                protein3 = new Protein("MA", "protein3");
            }

            var pepWithSetModifications1 = protein1.Digest(CommonParameters.DigestionParams, allKnownFixedModifications, variableModifications).First();

            var pepWithSetModifications2 = protein2.Digest(CommonParameters.DigestionParams, allKnownFixedModifications, variableModifications).First();

            var pepWithSetModifications3 = protein3.Digest(CommonParameters.DigestionParams, allKnownFixedModifications, variableModifications).Last();

            CompactPeptide compactPeptide1 = new CompactPeptide(pepWithSetModifications1, TerminusType.None);
            CompactPeptide compactPeptideDuplicate = new CompactPeptide(pepWithSetModifications2, TerminusType.None);
            Assert.AreEqual(compactPeptide1, compactPeptideDuplicate);
            CompactPeptide compactPeptide2 = new CompactPeptide(pepWithSetModifications3, TerminusType.None);

            string fullFilePath = null;
            int precursorCharge = 0;
            TestDataFile testDataFile = new TestDataFile();
            MsDataScan mzLibScan = testDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(mzLibScan, 0, precursorCharge, fullFilePath);
            int scanIndex = 0;
            double score = 0;
            int notch = 0;
            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(compactPeptide1, notch, score, scanIndex, scan, CommonParameters.DigestionParams);
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(compactPeptide1, notch, score, scanIndex, scan, CommonParameters.DigestionParams);
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            PeptideSpectralMatch psm3 = new PeptideSpectralMatch(compactPeptide2, notch, score, scanIndex, scan, CommonParameters.DigestionParams);
            psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            var newPsms = new List<PeptideSpectralMatch>
            {
                psm1,
                psm2,
                psm3
            };

            MassDiffAcceptor massDiffAcceptors = new SinglePpmAroundZeroSearchMode(5);
            SequencesToActualProteinPeptidesEngine stappe = new SequencesToActualProteinPeptidesEngine(newPsms, new List<Protein> { protein1, protein2, protein3 },
                allKnownFixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>());

            var haha = (SequencesToActualProteinPeptidesEngineResults)stappe.Run();
            var compactPeptideToProteinPeptideMatching = haha.CompactPeptideToProteinPeptideMatching;

            Assert.AreEqual(2, compactPeptideToProteinPeptideMatching.Count);

            psm1.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            bool noOneHitWonders = false;

            return new Tuple<List<PeptideSpectralMatch>, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>, MassDiffAcceptor, bool, CompactPeptideBase, CompactPeptideBase>
            (
                newPsms, compactPeptideToProteinPeptideMatching, massDiffAcceptors, noOneHitWonders, compactPeptide1, compactPeptide2
            );
        }       
    }
}