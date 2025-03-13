using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System;
using System.Collections.Generic;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;

namespace Test
{
    [TestFixture]
    internal class PeptideSpectralMatchTest
    {
        [Test]
        public static void GetAminoAcidCoverageTest()
        {
            CommonParameters commonParams = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            MsDataScan scanNumberOne = new MsDataScan(new MzSpectrum(new double[] { 10 }, new double[] { 1 }, false), 1, 2, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", 10, 2, 100, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass ms2ScanOneMzTen = new Ms2ScanWithSpecificMass(scanNumberOne, 10, 2, "File", new CommonParameters());

            string sequence = "";
            Dictionary<string, Modification> allKnownMods = new();
            int numFixedMods = 0;
            DigestionParams digestionParams = new DigestionParams();
            Protein myProtein = new Protein(sequence, "ACCESSION");
            int oneBasedStartResidueInProtein = 0;
            int oneBasedEndResidueInProtein = Math.Max(sequence.Length, 0);
            int missedCleavages = 0;
            CleavageSpecificity cleavageSpecificity = CleavageSpecificity.Full;
            string peptideDescription = null;
            string pairedTargetDecoySequence = null;

            PeptideWithSetModifications pwsmNoBaseSequence = new(sequence, allKnownMods, numFixedMods, digestionParams, myProtein,
                oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity,
                peptideDescription, pairedTargetDecoySequence);
            PeptideSpectralMatch psmNoBaseSequenceNoMFI = new(pwsmNoBaseSequence, 0, 10, 0, ms2ScanOneMzTen, commonParams,
                new List<MatchedFragmentIon>());
            psmNoBaseSequenceNoMFI.ResolveAllAmbiguities();

            //PSM has neither sequence nor matched fragment ions
            var b = psmNoBaseSequenceNoMFI.BaseSequence;
            Assert.That(b, Is.EqualTo(""));
            var m = psmNoBaseSequenceNoMFI.MatchedFragmentIons;
            Assert.That(m.Count, Is.EqualTo(0));
            psmNoBaseSequenceNoMFI.GetAminoAcidCoverage();

            sequence = "PEPTIDE";
            oneBasedEndResidueInProtein = Math.Max(sequence.Length, 0);
            myProtein = new Protein(sequence, "ACCESSION");
            var test = new PeptideWithSetModifications(sequence, allKnownMods);
            PeptideWithSetModifications pwsmBaseSequence = new PeptideWithSetModifications(sequence, allKnownMods, numFixedMods, digestionParams, myProtein,
                oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity,
                peptideDescription, pairedTargetDecoySequence);
            PeptideSpectralMatch psmBaseSequenceNoMFI = new(pwsmBaseSequence, 0, 10, 0, ms2ScanOneMzTen, commonParams,
                new List<MatchedFragmentIon>());

            //PSM has sequence but does not have matched fragment ions
            psmBaseSequenceNoMFI.ResolveAllAmbiguities();
            b = psmBaseSequenceNoMFI.BaseSequence;
            Assert.That(b, Is.EqualTo(sequence));
            m = psmBaseSequenceNoMFI.MatchedFragmentIons;
            Assert.That(m.Count, Is.EqualTo(0));
            psmBaseSequenceNoMFI.GetAminoAcidCoverage();

            //PSM has no sequence but does have matched fragment ions
            Product ntProduct = new Product(ProductType.y, FragmentationTerminus.N, 1, 1, 1, 0);
            MatchedFragmentIon mfi = new(ntProduct, 1, 1, 1);
            List<MatchedFragmentIon> mfiList = new List<MatchedFragmentIon>() { mfi };
            PeptideSpectralMatch psmNoBaseSequenceMFI = new(pwsmNoBaseSequence, 0, 10, 0, ms2ScanOneMzTen, commonParams,
                mfiList);
            psmNoBaseSequenceMFI.ResolveAllAmbiguities();

            b = psmNoBaseSequenceMFI.BaseSequence;
            Assert.That(b, Is.EqualTo(""));
            m = psmNoBaseSequenceMFI.MatchedFragmentIons;
            Assert.That(m.Count, Is.EqualTo(1));
            psmNoBaseSequenceMFI.GetAminoAcidCoverage();

            //PSM has sequence and matched fragment ions
            PeptideSpectralMatch psmBaseSequenceMFI = new(pwsmBaseSequence, 0, 10, 0, ms2ScanOneMzTen, commonParams,
                mfiList);
            psmBaseSequenceMFI.ResolveAllAmbiguities();

            b = psmBaseSequenceMFI.BaseSequence;
            Assert.That(b, Is.EqualTo("PEPTIDE"));
            m = psmBaseSequenceMFI.MatchedFragmentIons;
            Assert.That(m.Count, Is.EqualTo(1));
            psmBaseSequenceMFI.GetAminoAcidCoverage();

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>() { psmNoBaseSequenceNoMFI, psmBaseSequenceNoMFI, psmNoBaseSequenceMFI, psmBaseSequenceMFI };

            foreach (var psm in psms)
            {
                psm.ToString();
            }
        }
    }
}