using EngineLayer.Indexing;
using EngineLayer;
using EngineLayer.NonSpecificEnzymeSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using Chemistry;

namespace Test
{
    internal class NonSpecificSearchEngineTests
    {
        [Test]
        public void CheckDissociationTypeSupport_AddCompIons_NotSupported_ThrowsNotImplementedException()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                SearchTarget = true,
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                LocalFdrCategories = new List<FdrCategory>
                {
                    FdrCategory.NonSpecific
                }
            };
            DigestionParams dp = new DigestionParams("singleC", minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus.C, searchModeType: CleavageSpecificity.None);

            //Presently, the NonSpecificEnzymeSearchEngine does not support AddCompIons with Dissociation Type IRMPD
            CommonParameters CommonParameters = new CommonParameters(
                dissociationType: DissociationType.IRMPD,
                digestionParams: dp,
                scoreCutoff: 5,
                precursorMassTolerance: new PpmTolerance(5),
                addCompIons: true);

            PeptideWithSetModifications guiltyPwsm = new PeptideWithSetModifications("DQPKLLGIETPLPKKE", null);
            var fragments = new List<Product>();
            guiltyPwsm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, fragments);

            var myMsDataFile = new TestDataFile(guiltyPwsm.MonoisotopicMass, fragments.Select(x => x.NeutralMass.ToMz(1)).ToArray());
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            fixedModifications.Add(mod2);

            var proteinList = new List<Protein> { new Protein("GGGGGCDQPKLLGIETPLPKKEGGGGG", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.None, CommonParameters,
                null, SearchParameters.MaxFragmentSize, true, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;
            var precursorIndexDict = indexResults.PrecursorIndex;

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            SpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[3][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            SpectralMatch[] allPsmsArray = allPsmsArrays[2];

            //coisolation index
            List<int>[] coisolationIndex = new List<int>[listOfSortedms2Scans.Length];
            for (int i = 0; i < listOfSortedms2Scans.Length; i++)
            {
                coisolationIndex[i] = new List<int> { i };
            }

            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, coisolationIndex, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, null, variableModifications, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());

            // Act & Assert
            //Throws errero because AddCompIons is not supported with Dissociation Type IRMPD
            Assert.Throws<NotImplementedException>(() => engine.Run());
        }

        [Test]
        public void TestSnesScoringWithCompIons()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                SearchTarget = true,
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                LocalFdrCategories = new List<FdrCategory>
                {
                    FdrCategory.NonSpecific
                }
            };
            DigestionParams dp = new DigestionParams("singleC", minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus.C, searchModeType: CleavageSpecificity.None);

            //Presently, the NonSpecificEnzymeSearchEngine does not support AddCompIons with Dissociation Type IRMPD
            CommonParameters CommonParameters = new CommonParameters(
                dissociationType: DissociationType.LowCID,
                digestionParams: dp,
                scoreCutoff: 5,
                precursorMassTolerance: new PpmTolerance(5),
                addCompIons: true);

            PeptideWithSetModifications guiltyPwsm = new PeptideWithSetModifications("DQPKLLGIETPLPKKE", null);
            var fragments = new List<Product>();
            guiltyPwsm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, fragments);

            var myMsDataFile = new TestDataFile(guiltyPwsm.MonoisotopicMass, fragments.Select(x => x.NeutralMass.ToMz(1)).ToArray());
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            fixedModifications.Add(mod2);

            var proteinList = new List<Protein> { new Protein("GGGGGCDQPKLLGIETPLPKKEGGGGG", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.None, CommonParameters,
                null, SearchParameters.MaxFragmentSize, true, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            //This next line is the key bit to hit the code path that is being tested
            fragmentIndexDict[1655] = new List<int> { 1 };

            var precursorIndexDict = indexResults.PrecursorIndex;

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            SpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[3][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            SpectralMatch[] allPsmsArray = allPsmsArrays[2];

            //coisolation index
            List<int>[] coisolationIndex = new List<int>[listOfSortedms2Scans.Length];
            for (int i = 0; i < listOfSortedms2Scans.Length; i++)
            {
                coisolationIndex[i] = new List<int> { i };
            }

            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, coisolationIndex, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, null, variableModifications, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());

            engine.Run();
        }
    }
}
