using Chemistry;
using EngineLayer;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Test
{
    [TestFixture]
    public static class EventArgsTest
    {

        [Test]
        public static void SingleEventArgsTest()
        {
            
            Protein parentProteinForMatch = new Protein("MEK", null);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            List<Modification> variableModifications = new List<Modification> { new Modification(_originalId: "21", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = parentProteinForMatch.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(4, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications[1];
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.GetAllScansList().Last(), pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());

            var theoreticalProducts = ps.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, theoreticalProducts, new CommonParameters());
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(ps, 0, 0, 2, scan, digestionParams, matchedIons);

            LocalizationEngine f = new LocalizationEngine(new List<PeptideSpectralMatch> { newPsm }, myMsDataFile, new CommonParameters(), new List<string>());

            var singleEngine= new SingleEngineEventArgs(f);
            Assert.That(singleEngine.MyEngine.Equals(f));

            var singleFile = new SingleFileEventArgs("", new List<string>());
            Assert.That(singleFile.WrittenFile.Equals(""));

            var stringList = new StringListEventArgs(new List<string>());
            var rr = (stringList.StringList.DefaultIfEmpty().First());
            Assert.That(stringList.StringList.DefaultIfEmpty().First() == null);
            
        }
    }
}
