using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;

namespace Test
{
    internal static class PsmTsvWriterTests
    {
        [Test]
        public static void ResolveModificationsTest()
        {
            double mass = 12.0 + new PeptideWithSetModifications(new Protein("LNLDLDND", "prot1"),new DigestionParams(),1,8,CleavageSpecificity.Full,"",0, new Dictionary<int, Modification>(),0,null).MonoisotopicMass.ToMz(1);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), mass, 1, "", new CommonParameters());

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);

            Dictionary<DissociationType, List<double>> NeutralLosses = new Dictionary<DissociationType, List<double>>();
            NeutralLosses.Add(DissociationType.HCD, new List<double> { 0 });

            Modification modFormula_C1 = new Modification(_originalId: "modC", _accession: "", _modificationType: "mt", _featureType: "", _target: motif1, _locationRestriction: "Anywhere.", _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("C1")), null, null, null, null, _neutralLosses: NeutralLosses, null, null);
            Modification modFormula_H1 = new Modification(_originalId: "modH", _accession: "", _modificationType: "mt", _featureType: "", _target: motif1, _locationRestriction: "Anywhere.", _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("H1")), null, null, null, null, _neutralLosses: NeutralLosses, null, null);

            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification>{ modFormula_C1, modFormula_H1 }},
            };
            Protein protein1 = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: oneBasedModifications);

            Dictionary<int, Modification> allModsOneIsNterminus1 = new Dictionary<int, Modification>
            {
                {2, modFormula_C1},
            };

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus1, 0);

            Dictionary<int, Modification> allModsOneIsNterminus2 = new Dictionary<int, Modification>
            {
                {2,modFormula_H1 },
            };

            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus2, 0);

            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1), scoreCutoff: 1);
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", CommonParameters));

            List<MatchedFragmentIon> mfi = new List<MatchedFragmentIon>();

            //we're adding a neutral loss of 5 to the product to make sure we hit the right spot in the unit test to add that loss to the product ion string
            Product p = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 5);
            mfi.Add(new MatchedFragmentIon(ref p, 1, 1, 1));
            PeptideSpectralMatch myPsm = new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, new CommonParameters(), mfi);

            myPsm.AddOrReplace(pwsm2, 10, 0, true, mfi, 10);

            myPsm.ResolveAllAmbiguities();

            //Here we have a situation where there are two mods at the same position with different chemical formuala. They cannot be resolved and so the return value is null.
            Assert.IsNull(myPsm.ModsChemicalFormula);

            string myPsmString = myPsm.ToString();
            string[] myPsmStringSplit = myPsmString.Split('\t');
            string ppmErrorString = myPsmStringSplit[24];

            //The two different mods produce two separate mass errors, which are both then reported
            Assert.AreEqual("0.00|11801.30", ppmErrorString);

            //Make sure we see produt ion neutral losses in the output.
            string matchedIonSeries = myPsmStringSplit[39];
            Assert.AreEqual("[(b1-5.00)+1]", matchedIonSeries);


            //removing one of the peptides to reset for the next test
            myPsm.RemoveThisAmbiguousPeptide(0, pwsm2);

            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus1, 0);
            myPsm.AddOrReplace(pwsm3, 10, 0, true, mfi, 10);

            myPsm.ResolveAllAmbiguities();

            //Now we have removed one of the peptides with a different chemical formual and replaced it with a mod that has the same chemical formula as the remaining original best peptide
            //Here we have a situation where there are two mods at the same position have the same chemical formuala and they can be resolved and so the return value the chemical formual of the mod.
            Assert.AreEqual("C", myPsm.ModsChemicalFormula.Formula.ToString());

            myPsmString = myPsm.ToString();
            myPsmStringSplit = myPsmString.Split('\t');
            ppmErrorString = myPsmStringSplit[24];

            Assert.AreEqual("0.00", ppmErrorString);
        }
    }
}