using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using NUnit.Framework.Legacy;
using Omics.Digestion;
using Omics.Fragmentation.Peptide;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class ParameterTest
    {
        [Test]
        public static void TestFileSpecificParametersClone()
        {
            var fileSpecificToml = Toml.ReadFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "testFileSpecfic.toml"), MetaMorpheusTask.tomlConfig);
            FileSpecificParameters fsp = new FileSpecificParameters(fileSpecificToml);
            FileSpecificParameters fspClone = fsp.Clone();
            Assert.That(fsp.DissociationType, Is.EqualTo(fspClone.DissociationType));
            Assert.That(fsp.MaxMissedCleavages, Is.EqualTo(fspClone.MaxMissedCleavages));
            Assert.That(fsp.MaxModsForPeptide, Is.EqualTo(fspClone.MaxModsForPeptide));
            Assert.That(fsp.MaxPeptideLength, Is.EqualTo(fspClone.MaxPeptideLength));
            Assert.That(fsp.MinPeptideLength, Is.EqualTo(fspClone.MinPeptideLength));
            Assert.That(fsp.PrecursorMassTolerance, Is.EqualTo(fspClone.PrecursorMassTolerance));
            Assert.That(fsp.ProductMassTolerance, Is.EqualTo(fspClone.ProductMassTolerance));
            Assert.That(fsp.DigestionAgent, Is.EqualTo(fspClone.DigestionAgent));
            Assert.That(fsp.SeparationType, Is.EqualTo(fspClone.SeparationType));
            CollectionAssert.AreEqual(fsp.CustomIons, fspClone.CustomIons);
        }

        [Test]
        //This test exists because file specific parameters code has a tendency to overwrite common parameters and make life horrible 
        public static void TestFileSpecifcParameterOverwrite()
        {
            CommonParameters defaultParameters = new CommonParameters();

            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = 
                new List<ProductType> { ProductType.b, ProductType.y, ProductType.c };
            CommonParameters notDefaultParameters = new CommonParameters(
                dissociationType: DissociationType.ETD,
                doPrecursorDeconvolution: false,
                useProvidedPrecursorInfo: false,
                deconvolutionIntensityRatio: 69,
                deconvolutionMaxAssumedChargeState: 69,
                reportAllAmbiguity: false,
                addCompIons: true,
                totalPartitions: 69,
                scoreCutoff: 69,
                numberOfPeaksToKeepPerWindow: 69,
                minimumAllowedIntensityRatioToBasePeak: 69,
                trimMs1Peaks: true,
                trimMsMsPeaks: false,
                productMassTolerance: new AbsoluteTolerance(69),
                precursorMassTolerance: new AbsoluteTolerance(69),
                deconvolutionMassTolerance: new AbsoluteTolerance(69),
                maxThreadsToUsePerFile: 69,
                digestionParams: new DigestionParams(
                    protease: "Asp-N",
                    maxMissedCleavages: 69,
                    minPeptideLength: 1,
                    maxPeptideLength: 69,
                    maxModificationIsoforms: 69,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                    maxModsForPeptides: 69,
                    searchModeType: CleavageSpecificity.Semi,
                    fragmentationTerminus: FragmentationTerminus.C
                    ),
                listOfModsVariable: new List<(string, string)> { ("asdf", "asdf") },
                listOfModsFixed: new List<(string, string)> { ("asdf", "asdf") }
                );

            //check that the defaults are not the same as the not defaults.
            //IF ONE OF THESE FAILS, PLEASE UPDATE THE "notDefaultParameters"
            Assert.That(defaultParameters.DissociationType, Is.Not.EqualTo(notDefaultParameters.DissociationType));
            Assert.That(defaultParameters.DoPrecursorDeconvolution, Is.Not.EqualTo(notDefaultParameters.DoPrecursorDeconvolution));
            Assert.That(defaultParameters.UseProvidedPrecursorInfo, Is.Not.EqualTo(notDefaultParameters.UseProvidedPrecursorInfo));
            Assert.That(defaultParameters.DeconvolutionIntensityRatio, Is.Not.EqualTo(notDefaultParameters.DeconvolutionIntensityRatio));
            Assert.That(defaultParameters.DeconvolutionMaxAssumedChargeState, Is.Not.EqualTo(notDefaultParameters.DeconvolutionMaxAssumedChargeState));
            Assert.That(defaultParameters.ReportAllAmbiguity, Is.Not.EqualTo(notDefaultParameters.ReportAllAmbiguity));
            Assert.That(defaultParameters.AddCompIons, Is.Not.EqualTo(notDefaultParameters.AddCompIons));
            Assert.That(defaultParameters.TotalPartitions, Is.Not.EqualTo(notDefaultParameters.TotalPartitions));
            Assert.That(defaultParameters.ScoreCutoff, Is.Not.EqualTo(notDefaultParameters.ScoreCutoff));
            Assert.That(defaultParameters.NumberOfPeaksToKeepPerWindow, Is.Not.EqualTo(notDefaultParameters.NumberOfPeaksToKeepPerWindow));
            Assert.That(defaultParameters.MinimumAllowedIntensityRatioToBasePeak, Is.Not.EqualTo(notDefaultParameters.MinimumAllowedIntensityRatioToBasePeak));
            Assert.That(defaultParameters.TrimMs1Peaks, Is.Not.EqualTo(notDefaultParameters.TrimMs1Peaks));
            Assert.That(defaultParameters.TrimMsMsPeaks, Is.Not.EqualTo(notDefaultParameters.TrimMsMsPeaks));
            Assert.That(defaultParameters.ProductMassTolerance, Is.Not.EqualTo(notDefaultParameters.ProductMassTolerance));
            Assert.That(defaultParameters.PrecursorMassTolerance, Is.Not.EqualTo(notDefaultParameters.PrecursorMassTolerance));
            Assert.That(defaultParameters.DeconvolutionMassTolerance, Is.Not.EqualTo(notDefaultParameters.DeconvolutionMassTolerance));
            Assert.That(defaultParameters.MaxThreadsToUsePerFile, Is.Not.EqualTo(notDefaultParameters.MaxThreadsToUsePerFile));
            Assert.That(defaultParameters.DigestionParams, Is.Not.EqualTo(notDefaultParameters.DigestionParams));
            Assert.That(defaultParameters.ListOfModsVariable, Is.Not.EqualTo(notDefaultParameters.ListOfModsVariable));
            Assert.That(defaultParameters.ListOfModsFixed, Is.Not.EqualTo(notDefaultParameters.ListOfModsFixed));
            Assert.That(defaultParameters.CustomIons, Is.Not.EqualTo(notDefaultParameters.CustomIons));

            FileSpecificParameters emptyFileSpecificParameters = new FileSpecificParameters();
            CommonParameters updatedParameters = MetaMorpheusTask.SetAllFileSpecificCommonParams(notDefaultParameters, emptyFileSpecificParameters);

            //CHECK THAT NOTHING CHANGED
            Assert.That(updatedParameters.DissociationType, Is.EqualTo(notDefaultParameters.DissociationType));
            Assert.That(updatedParameters.DoPrecursorDeconvolution, Is.EqualTo(notDefaultParameters.DoPrecursorDeconvolution));
            Assert.That(updatedParameters.UseProvidedPrecursorInfo, Is.EqualTo(notDefaultParameters.UseProvidedPrecursorInfo));
            Assert.That(updatedParameters.DeconvolutionIntensityRatio, Is.EqualTo(notDefaultParameters.DeconvolutionIntensityRatio));
            Assert.That(updatedParameters.DeconvolutionMaxAssumedChargeState, Is.EqualTo(notDefaultParameters.DeconvolutionMaxAssumedChargeState));
            Assert.That(updatedParameters.ReportAllAmbiguity, Is.EqualTo(notDefaultParameters.ReportAllAmbiguity));
            Assert.That(updatedParameters.AddCompIons, Is.EqualTo(notDefaultParameters.AddCompIons));
            Assert.That(updatedParameters.TotalPartitions, Is.EqualTo(notDefaultParameters.TotalPartitions));
            Assert.That(updatedParameters.ScoreCutoff, Is.EqualTo(notDefaultParameters.ScoreCutoff));
            Assert.That(updatedParameters.NumberOfPeaksToKeepPerWindow, Is.EqualTo(notDefaultParameters.NumberOfPeaksToKeepPerWindow));
            Assert.That(updatedParameters.MinimumAllowedIntensityRatioToBasePeak, Is.EqualTo(notDefaultParameters.MinimumAllowedIntensityRatioToBasePeak));
            Assert.That(updatedParameters.TrimMs1Peaks, Is.EqualTo(notDefaultParameters.TrimMs1Peaks));
            Assert.That(updatedParameters.TrimMsMsPeaks, Is.EqualTo(notDefaultParameters.TrimMsMsPeaks));
            Assert.That(updatedParameters.ProductMassTolerance, Is.EqualTo(notDefaultParameters.ProductMassTolerance));
            Assert.That(updatedParameters.PrecursorMassTolerance, Is.EqualTo(notDefaultParameters.PrecursorMassTolerance));
            Assert.That(updatedParameters.DeconvolutionMassTolerance, Is.EqualTo(notDefaultParameters.DeconvolutionMassTolerance));
            Assert.That(updatedParameters.MaxThreadsToUsePerFile, Is.EqualTo(notDefaultParameters.MaxThreadsToUsePerFile));
            Assert.That(updatedParameters.DigestionParams, Is.EqualTo(notDefaultParameters.DigestionParams));
            Assert.That(updatedParameters.ListOfModsVariable, Is.EqualTo(notDefaultParameters.ListOfModsVariable));
            Assert.That(updatedParameters.ListOfModsFixed, Is.EqualTo(notDefaultParameters.ListOfModsFixed));
            Assert.That(updatedParameters.CustomIons, Is.EqualTo(notDefaultParameters.CustomIons));

            FileSpecificParameters basicFileSpecificParameters = new FileSpecificParameters
            {
                PrecursorMassTolerance = new PpmTolerance(10),
                ProductMassTolerance = new PpmTolerance(30),
                DigestionAgent = new Protease("Arg-C", CleavageSpecificity.Full, null, null, new List<DigestionMotif> { new DigestionMotif("K", null, 1, "") }),
                MinPeptideLength = 1,
                MaxPeptideLength = 50,
                MaxMissedCleavages = 2,
                MaxModsForPeptide = 1,
                DissociationType = DissociationType.CID,
                CustomIons = new List<ProductType> { ProductType.b, ProductType.y }
            };
            updatedParameters = MetaMorpheusTask.SetAllFileSpecificCommonParams(notDefaultParameters, basicFileSpecificParameters);
            //CHECK THAT SOMETHINGS CHANGED AND OTHERS DIDN'T
            Assert.That(updatedParameters.DissociationType, Is.EqualTo(basicFileSpecificParameters.DissociationType));
            Assert.That(updatedParameters.ProductMassTolerance, Is.EqualTo(basicFileSpecificParameters.ProductMassTolerance));
            Assert.That(updatedParameters.PrecursorMassTolerance, Is.EqualTo(basicFileSpecificParameters.PrecursorMassTolerance));
            Assert.That(updatedParameters.DigestionParams.MaxMods, Is.EqualTo(basicFileSpecificParameters.MaxModsForPeptide));
            Assert.That(updatedParameters.DigestionParams.MaxMissedCleavages, Is.EqualTo(basicFileSpecificParameters.MaxMissedCleavages));
            Assert.That(updatedParameters.DigestionParams.MinLength, Is.EqualTo(basicFileSpecificParameters.MinPeptideLength));
            Assert.That(updatedParameters.DigestionParams.MaxLength, Is.EqualTo(basicFileSpecificParameters.MaxPeptideLength));
            Assert.That(updatedParameters.DigestionParams.DigestionAgent, Is.EqualTo(basicFileSpecificParameters.DigestionAgent));
            Assert.That(updatedParameters.CustomIons, Is.EqualTo(basicFileSpecificParameters.CustomIons));

            Assert.That(updatedParameters.DoPrecursorDeconvolution, Is.EqualTo(notDefaultParameters.DoPrecursorDeconvolution));
            Assert.That(updatedParameters.UseProvidedPrecursorInfo, Is.EqualTo(notDefaultParameters.UseProvidedPrecursorInfo));
            Assert.That(updatedParameters.DeconvolutionIntensityRatio, Is.EqualTo(notDefaultParameters.DeconvolutionIntensityRatio));
            Assert.That(updatedParameters.DeconvolutionMaxAssumedChargeState, Is.EqualTo(notDefaultParameters.DeconvolutionMaxAssumedChargeState));
            Assert.That(updatedParameters.ReportAllAmbiguity, Is.EqualTo(notDefaultParameters.ReportAllAmbiguity));
            Assert.That(updatedParameters.AddCompIons, Is.EqualTo(notDefaultParameters.AddCompIons));
            Assert.That(updatedParameters.TotalPartitions, Is.EqualTo(notDefaultParameters.TotalPartitions));
            Assert.That(updatedParameters.ScoreCutoff, Is.EqualTo(notDefaultParameters.ScoreCutoff));
            Assert.That(updatedParameters.NumberOfPeaksToKeepPerWindow, Is.EqualTo(notDefaultParameters.NumberOfPeaksToKeepPerWindow));
            Assert.That(updatedParameters.MinimumAllowedIntensityRatioToBasePeak, Is.EqualTo(notDefaultParameters.MinimumAllowedIntensityRatioToBasePeak));
            Assert.That(updatedParameters.TrimMs1Peaks, Is.EqualTo(notDefaultParameters.TrimMs1Peaks));
            Assert.That(updatedParameters.TrimMsMsPeaks, Is.EqualTo(notDefaultParameters.TrimMsMsPeaks));
            Assert.That(updatedParameters.DeconvolutionMassTolerance, Is.EqualTo(notDefaultParameters.DeconvolutionMassTolerance));
            Assert.That(updatedParameters.MaxThreadsToUsePerFile, Is.EqualTo(notDefaultParameters.MaxThreadsToUsePerFile));
            Assert.That(((DigestionParams)updatedParameters.DigestionParams).InitiatorMethionineBehavior, Is.EqualTo(((DigestionParams)notDefaultParameters.DigestionParams).InitiatorMethionineBehavior));
            Assert.That(updatedParameters.ListOfModsVariable, Is.EqualTo(notDefaultParameters.ListOfModsVariable));
            Assert.That(updatedParameters.ListOfModsFixed, Is.EqualTo(notDefaultParameters.ListOfModsFixed));

        }
    }
}