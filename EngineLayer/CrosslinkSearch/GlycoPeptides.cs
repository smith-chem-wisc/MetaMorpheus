using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class GlycoPeptides
    {
        public static IEnumerable<Tuple<int, List<Product>>> GlyGetTheoreticalFragments(DissociationType dissociationType,
    List<int> possibleModPositions, PeptideWithSetModifications peptide, Glycan glycan)
        {
            Modification modification = GlycanToModification(glycan);

            foreach (var position in possibleModPositions)
            {               
                Dictionary<int, Modification> testMods = new Dictionary<int, Modification> { { position + 1, modification } };

                var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                    peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

                List<Product> theoreticalProducts = testPeptide.Fragment(dissociationType, FragmentationTerminus.Both).ToList();

                yield return new Tuple<int, List<Product>>(position, theoreticalProducts);
            }

        }

        public static Modification GlycanToModification(Glycan glycan)
        {          
            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();
            List<double> lossMasses = glycan.Ions.Where(p => p.IonMass <= 1000).Select(p => glycan.Mass - p.IonMass).OrderBy(p => p).ToList();
            lossMasses.Add(glycan.Mass - 83.038194);
            lossMasses.Add(glycan.Mass);        
            neutralLosses.Add(DissociationType.HCD, lossMasses);
            neutralLosses.Add(DissociationType.CID, lossMasses);
            //TO DO: add diagnosticIons as a property for Glycan
            Dictionary<DissociationType, List<double>> diagnosticIons = new Dictionary<DissociationType, List<double>>();
            diagnosticIons.Add(DissociationType.HCD, glycan.GetDiagnosticIons().Values.ToList());
            diagnosticIons.Add(DissociationType.CID, glycan.GetDiagnosticIons().Values.ToList());

            Modification modification = new Modification(_originalId:"Glycan", _monoisotopicMass: glycan.Mass, _neutralLosses: neutralLosses, _diagnosticIons : diagnosticIons);
            return modification;
        }
    }
}
