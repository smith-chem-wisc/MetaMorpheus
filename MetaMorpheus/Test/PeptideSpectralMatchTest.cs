using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

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
            Dictionary< string, Modification > allKnownMods = new();
            int numFixedMods = 0;
            DigestionParams digestionParams = new DigestionParams();
            Protein myProtein = new Protein(sequence,"ACCESSION");
            int oneBasedStartResidueInProtein = 0;
            int oneBasedEndResidueInProtein = Math.Max(sequence.Length,0);
            int missedCleavages = 0;
            CleavageSpecificity cleavageSpecificity = CleavageSpecificity.Full;
            string peptideDescription = null;
            int? pairedTargetDecoyHash = null;

            //PeptideWithSetModifications pwsmNoBaseSequence = new(sequence, allKnownMods, numFixedMods, digestionParams, myProtein,
            //    oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity,
            //    peptideDescription, pairedTargetDecoyHash);
            //PeptideSpectralMatch psmNoBaseSequence = new(pwsmNoBaseSequence, 0, 10, 0, ms2ScanOneMzTen, commonParams,
            //    new List<MatchedFragmentIon>());

            //var b = psmNoBaseSequence.BaseSequence;
            //var m = psmNoBaseSequence.MatchedFragmentIons;
            //psmNoBaseSequence.GetAminoAcidCoverage();

            sequence = "PEPTIDE";
            oneBasedEndResidueInProtein = Math.Max(sequence.Length , 0);
            myProtein = new Protein(sequence, "ACCESSION");
            PeptideWithSetModifications pwsmBaseSequence = new(sequence, allKnownMods, numFixedMods, digestionParams, myProtein,
                oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity,
                peptideDescription, pairedTargetDecoyHash);
            PeptideSpectralMatch psmBaseSequence = new(pwsmBaseSequence, 0, 10, 0, ms2ScanOneMzTen, commonParams,
                new List<MatchedFragmentIon>());

            var b = psmBaseSequence.BaseSequence;
            var m = psmBaseSequence.MatchedFragmentIons;
            psmBaseSequence.GetAminoAcidCoverage();

            Assert.IsTrue(false);
        }
    }
}
