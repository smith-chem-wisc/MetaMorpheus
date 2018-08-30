using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class StefanParsimonyTest
    {
        [Test]
        public static void ParsimonyTreatModificationsAsUnique()
        {
            bool modPeptidesAreUnique = true;

            // set up mods
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            var mod = new Modification(_id: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);

            Protease protease = new Protease("kprotease", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            
            // modified version of protein
            var protein1 = new Protein("PEPTIDEM", "accession1");

            // unmodified version of protein
            var protein2 = new Protein("YYYKPEPTIDEM", "accession2");

            List<PeptideWithSetModifications> pwsmsFromProtein1 = protein1.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<Modification> { mod }, new List<Modification>()).ToList();  //this is a fixed mod
            List<PeptideWithSetModifications> pwsmsFromProtein2 = protein2.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<Modification>(), new List<Modification>()).ToList();

            // check to make sure mod is present
            PeptideWithSetModifications modifiedPeptide = pwsmsFromProtein1[0];
            PeptideWithSetModifications unmodifiedPeptide = pwsmsFromProtein2[1];

            Assert.That(!modifiedPeptide.Sequence.Equals(unmodifiedPeptide.Sequence)); // sequences should not be equal (one has a mod)
            Assert.That(modifiedPeptide.BaseSequence.Equals(unmodifiedPeptide.BaseSequence)); // base sequences should be equal
            Assert.That(modifiedPeptide.NumMods == 1); // methionine was oxidized on this protein
            Assert.That(unmodifiedPeptide.NumMods == 0); // there was no modification on this protein

            // build PSMs for parsimony
            List<PeptideSpectralMatch> psmsForParsimony = new List<PeptideSpectralMatch>();

            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(modifiedPeptide, 0, 10, 1, null, new DigestionParams(), new List<MatchedFragmentIon>());
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(unmodifiedPeptide, 0, 10, 2, null, new DigestionParams(), new List<MatchedFragmentIon>());
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

            //one peptide is unused
            //PeptideSpectralMatch psm3 = new PeptideSpectralMatch(pep2[0], 0, 10, 3, null, new DigestionParams());
            //psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            
            psmsForParsimony.Add(psm1);
            psmsForParsimony.Add(psm2);

            // apply parsimony
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(psmsForParsimony, modPeptidesAreUnique, new CommonParameters(), new List<string>());

            // because the two chosen peptides are the same, we should end up with both protein accessions still in the list
            var proteinParsimonyResult = (ProteinParsimonyResults)pae.Run();
            int countOfProteinGroups = proteinParsimonyResult.ProteinGroups.Count;
            
            // because modified peptides were considered as unique, then there should be two protein groups after parsimony, and one protein accession for each peptide
            Assert.That(countOfProteinGroups == 2);
            Assert.That(psm1.BestMatchingPeptideWithSetMods.First().Pwsm.Protein.Accession == "accession1");
            Assert.That(psm2.BestMatchingPeptideWithSetMods.First().Pwsm.Protein.Accession == "accession2");
        }

        [Test]
        public static void ParsimonyVariableDontTreatAsUnique()
        {
            bool modPeptidesAreUnique = false;

            // set up mods
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            var mod = new Modification(_id: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);

            Protease protease = new Protease("kprotease", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            // modified version of protein
            var protein1 = new Protein("PEPTIDEM", "accession1");

            // unmodified version of protein
            var protein2 = new Protein("YYYKPEPTIDEM", "accession2");

            List<PeptideWithSetModifications> pwsmsFromProtein1 = protein1.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<Modification> { mod }, new List<Modification>()).ToList();  //this is a fixed mod
            List<PeptideWithSetModifications> pwsmsFromProtein2 = protein2.Digest(new DigestionParams(protease: protease.Name, minPeptideLength: 1), new List<Modification>(), new List<Modification>()).ToList();

            // check to make sure mod is present
            PeptideWithSetModifications modifiedPeptide = pwsmsFromProtein1[0];
            PeptideWithSetModifications unmodifiedPeptide = pwsmsFromProtein2[1];

            Assert.That(!modifiedPeptide.Sequence.Equals(unmodifiedPeptide.Sequence)); // sequences should not be equal (one has a mod)
            Assert.That(modifiedPeptide.BaseSequence.Equals(unmodifiedPeptide.BaseSequence)); // base sequences should be equal
            Assert.That(modifiedPeptide.NumMods == 1); // methionine was oxidized on this protein
            Assert.That(unmodifiedPeptide.NumMods == 0); // there was no modification on this protein

            // build PSMs for parsimony
            List<PeptideSpectralMatch> psmsForParsimony = new List<PeptideSpectralMatch>();

            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(modifiedPeptide, 0, 10, 1, null, new DigestionParams(), new List<MatchedFragmentIon>());
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(unmodifiedPeptide, 0, 10, 2, null, new DigestionParams(), new List<MatchedFragmentIon>());
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

            //one peptide is unused
            //PeptideSpectralMatch psm3 = new PeptideSpectralMatch(pep2[0], 0, 10, 3, null, new DigestionParams());
            //psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

            psmsForParsimony.Add(psm1);
            psmsForParsimony.Add(psm2);

            // apply parsimony
            ProteinParsimonyEngine pae = new ProteinParsimonyEngine(psmsForParsimony, modPeptidesAreUnique, new CommonParameters(), new List<string>());

            // 
            var proteinParsimonyResult = (ProteinParsimonyResults)pae.Run();
            int countOfProteinGroups = proteinParsimonyResult.ProteinGroups.Count;

            // 
            Assert.That(countOfProteinGroups == 1);
            Assert.That(psm1.BestMatchingPeptideWithSetMods.Select(p => p.Pwsm.Protein.Accession).Distinct().Count() == 2);
            Assert.That(psm2.BestMatchingPeptideWithSetMods.Select(p => p.Pwsm.Protein.Accession).Distinct().Count() == 2);
        }
    }
}