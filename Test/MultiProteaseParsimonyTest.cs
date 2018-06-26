using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Text;
using EngineLayer;
using MassSpectrometry;
using Proteomics;
using System.Linq;
using IO.MzML;

namespace Test
{
    [TestFixture]
    class MultiProteaseParsimonyTest
    {
     
       [Test]
       public static void MultiProteaseCompactPeptideMatchingTest()
       {
            string[] sequences = {
                "---ABC",
                "--EFGABC",
            };
                
            IEnumerable<string> sequencesInducingCleavage = new List<string> {  "-" };
            IEnumerable<string> sequencesInducingCleavage2 = new List<string> { "G" };


            var protease = new Protease("test1", sequencesInducingCleavage, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            
            var protease2 = new Protease("test2", sequencesInducingCleavage2, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            
            var peptideList = new HashSet<PeptideWithSetModifications>();
           
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i<sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));

            DigestionParams digestionParams = new DigestionParams(protease: protease, MinPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2, MinPeptideLength: 1);
           
            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        
                        
                        
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        
                        case "ABC": peptideList.Add(peptide); break;
                        case "--EFG": peptideList.Add(peptide); break;
                        
                        

                    }
                }
                
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];
            
            PeptideWithSetModifications[] PWSM= new PeptideWithSetModifications[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i<peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }


            

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[3] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1] });
            dictionary.Add(peptides[2], new HashSet<PeptideWithSetModifications> { PWSM[2] });


            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");


            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC":
                            if (peptide.digestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.digestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }



                        case "EFGABC": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;





                        case "--EFG": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2)); break;



                    }
                }
            }

            
            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            List<DigestionParams> digestionParamsList = new List<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
           // digestionParamsList.Add(digestionParams3);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Oxidation of M", "Common Variable", motif2, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, false, null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();
          
            List < ProteinGroup > proteinGroups = proteinAnalysisResults.ProteinGroups;
            Assert.AreEqual(2, proteinGroups.Count);
            Assert.AreEqual(1, proteinGroups.ElementAt(0).AllPeptides.Count);
            Assert.AreEqual(3, proteinGroups.ElementAt(1).AllPeptides.Count);




        }
        
        
    }
}
