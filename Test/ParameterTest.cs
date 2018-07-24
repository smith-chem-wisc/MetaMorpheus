﻿using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class ParameterTest
    {
        [Test]
        //This test exists because file specific parameters code has a tendency to overwrite common parameters and make life horrible 
        public static void TestFileSpecifcParameterOverwrite()
        {
            CommonParameters defaultParameters = new CommonParameters();
            CommonParameters notDefaultParameters = new CommonParameters(
                bIons: false,
                yIons: false,
                zDotIons: true,
                cIons: true,
                doPrecursorDeconvolution: false,
                useProvidedPrecursorInfo: false,
                deconvolutionIntensityRatio: 69,
                deconvolutionMaxAssumedChargeState: 69,
                reportAllAmbiguity: false,
                addCompIons: true,
                totalPartitions: 69,
                scoreCutoff: 69,
                topNpeaks: 69,
                minRatio: 69,
                trimMs1Peaks: true,
                trimMsMsPeaks: false,
                useDeltaScore: true,
                calculateEValue: true,
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
                    semiProteaseDigestion: true,
                    terminusTypeSemiProtease: TerminusType.C
                    ),
                listOfModsVariable: new List<(string, string)> { ("asdf", "asdf") },
                listOfModsFixed: new List<(string, string)> { ("asdf", "asdf") }
                );

            //check that the defaults are not the same as the not defaults.
            //IF ONE OF THESE FAILS, PLEASE UPDATE THE "notDefaultParameters"
            Assert.AreNotEqual(defaultParameters.BIons, notDefaultParameters.BIons);
            Assert.AreNotEqual(defaultParameters.YIons, notDefaultParameters.YIons);
            Assert.AreNotEqual(defaultParameters.ZdotIons, notDefaultParameters.ZdotIons);
            Assert.AreNotEqual(defaultParameters.CIons, notDefaultParameters.CIons);
            Assert.AreNotEqual(defaultParameters.DoPrecursorDeconvolution, notDefaultParameters.DoPrecursorDeconvolution);
            Assert.AreNotEqual(defaultParameters.UseProvidedPrecursorInfo, notDefaultParameters.UseProvidedPrecursorInfo);
            Assert.AreNotEqual(defaultParameters.DeconvolutionIntensityRatio, notDefaultParameters.DeconvolutionIntensityRatio);
            Assert.AreNotEqual(defaultParameters.DeconvolutionMaxAssumedChargeState, notDefaultParameters.DeconvolutionMaxAssumedChargeState);
            Assert.AreNotEqual(defaultParameters.ReportAllAmbiguity, notDefaultParameters.ReportAllAmbiguity);
            Assert.AreNotEqual(defaultParameters.AddCompIons, notDefaultParameters.AddCompIons);
            Assert.AreNotEqual(defaultParameters.TotalPartitions, notDefaultParameters.TotalPartitions);
            Assert.AreNotEqual(defaultParameters.ScoreCutoff, notDefaultParameters.ScoreCutoff);
            Assert.AreNotEqual(defaultParameters.TopNpeaks, notDefaultParameters.TopNpeaks);
            Assert.AreNotEqual(defaultParameters.MinRatio, notDefaultParameters.MinRatio);
            Assert.AreNotEqual(defaultParameters.TrimMs1Peaks, notDefaultParameters.TrimMs1Peaks);
            Assert.AreNotEqual(defaultParameters.TrimMsMsPeaks, notDefaultParameters.TrimMsMsPeaks);
            Assert.AreNotEqual(defaultParameters.UseDeltaScore, notDefaultParameters.UseDeltaScore);
            Assert.AreNotEqual(defaultParameters.CalculateEValue, notDefaultParameters.CalculateEValue);
            Assert.AreNotEqual(defaultParameters.ProductMassTolerance, notDefaultParameters.ProductMassTolerance);
            Assert.AreNotEqual(defaultParameters.PrecursorMassTolerance, notDefaultParameters.PrecursorMassTolerance);
            Assert.AreNotEqual(defaultParameters.DeconvolutionMassTolerance, notDefaultParameters.DeconvolutionMassTolerance);
            Assert.AreNotEqual(defaultParameters.MaxThreadsToUsePerFile, notDefaultParameters.MaxThreadsToUsePerFile);
            Assert.AreNotEqual(defaultParameters.DigestionParams, notDefaultParameters.DigestionParams);
            Assert.AreNotEqual(defaultParameters.ListOfModsVariable, notDefaultParameters.ListOfModsVariable);
            Assert.AreNotEqual(defaultParameters.ListOfModsFixed, notDefaultParameters.ListOfModsFixed);

            FileSpecificParameters emptyFileSpecificParameters = new FileSpecificParameters();
            CommonParameters updatedParameters = MetaMorpheusTask.SetAllFileSpecificCommonParams(notDefaultParameters, emptyFileSpecificParameters);

            //CHECK THAT NOTHING CHANGED
            Assert.AreEqual(updatedParameters.BIons, notDefaultParameters.BIons);
            Assert.AreEqual(updatedParameters.YIons, notDefaultParameters.YIons);
            Assert.AreEqual(updatedParameters.ZdotIons, notDefaultParameters.ZdotIons);
            Assert.AreEqual(updatedParameters.CIons, notDefaultParameters.CIons);
            Assert.AreEqual(updatedParameters.DoPrecursorDeconvolution, notDefaultParameters.DoPrecursorDeconvolution);
            Assert.AreEqual(updatedParameters.UseProvidedPrecursorInfo, notDefaultParameters.UseProvidedPrecursorInfo);
            Assert.AreEqual(updatedParameters.DeconvolutionIntensityRatio, notDefaultParameters.DeconvolutionIntensityRatio);
            Assert.AreEqual(updatedParameters.DeconvolutionMaxAssumedChargeState, notDefaultParameters.DeconvolutionMaxAssumedChargeState);
            Assert.AreEqual(updatedParameters.ReportAllAmbiguity, notDefaultParameters.ReportAllAmbiguity);
            Assert.AreEqual(updatedParameters.AddCompIons, notDefaultParameters.AddCompIons);
            Assert.AreEqual(updatedParameters.TotalPartitions, notDefaultParameters.TotalPartitions);
            Assert.AreEqual(updatedParameters.ScoreCutoff, notDefaultParameters.ScoreCutoff);
            Assert.AreEqual(updatedParameters.TopNpeaks, notDefaultParameters.TopNpeaks);
            Assert.AreEqual(updatedParameters.MinRatio, notDefaultParameters.MinRatio);
            Assert.AreEqual(updatedParameters.TrimMs1Peaks, notDefaultParameters.TrimMs1Peaks);
            Assert.AreEqual(updatedParameters.TrimMsMsPeaks, notDefaultParameters.TrimMsMsPeaks);
            Assert.AreEqual(updatedParameters.UseDeltaScore, notDefaultParameters.UseDeltaScore);
            Assert.AreEqual(updatedParameters.CalculateEValue, notDefaultParameters.CalculateEValue);
            Assert.AreEqual(updatedParameters.ProductMassTolerance, notDefaultParameters.ProductMassTolerance);
            Assert.AreEqual(updatedParameters.PrecursorMassTolerance, notDefaultParameters.PrecursorMassTolerance);
            Assert.AreEqual(updatedParameters.DeconvolutionMassTolerance, notDefaultParameters.DeconvolutionMassTolerance);
            Assert.AreEqual(updatedParameters.MaxThreadsToUsePerFile, notDefaultParameters.MaxThreadsToUsePerFile);
            Assert.AreEqual(updatedParameters.DigestionParams, notDefaultParameters.DigestionParams);
            Assert.AreEqual(updatedParameters.ListOfModsVariable, notDefaultParameters.ListOfModsVariable);
            Assert.AreEqual(updatedParameters.ListOfModsFixed, notDefaultParameters.ListOfModsFixed);

            FileSpecificParameters basicFileSpecificParameters = new FileSpecificParameters
            {
                PrecursorMassTolerance = new PpmTolerance(10),
                ProductMassTolerance = new PpmTolerance(30),
                Protease = new Protease("Arg-C", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null),
                MinPeptideLength = 1,
                MaxPeptideLength = 50,
                MaxMissedCleavages = 2,
                MaxModsForPeptide = 1,
                BIons = true,
                YIons = true,
                CIons = false,
                ZdotIons = false
            };
            updatedParameters = MetaMorpheusTask.SetAllFileSpecificCommonParams(notDefaultParameters, basicFileSpecificParameters);
            //CHECK THAT SOMETHINGS CHANGED AND OTHERS DIDN'T
            Assert.AreEqual(updatedParameters.BIons, basicFileSpecificParameters.BIons);
            Assert.AreEqual(updatedParameters.YIons, basicFileSpecificParameters.YIons);
            Assert.AreEqual(updatedParameters.ZdotIons, basicFileSpecificParameters.ZdotIons);
            Assert.AreEqual(updatedParameters.ProductMassTolerance, basicFileSpecificParameters.ProductMassTolerance);
            Assert.AreEqual(updatedParameters.PrecursorMassTolerance, basicFileSpecificParameters.PrecursorMassTolerance);
            Assert.AreEqual(updatedParameters.CIons, basicFileSpecificParameters.CIons);
            Assert.AreEqual(updatedParameters.DigestionParams.MaxModsForPeptide, basicFileSpecificParameters.MaxModsForPeptide);
            Assert.AreEqual(updatedParameters.DigestionParams.MaxMissedCleavages, basicFileSpecificParameters.MaxMissedCleavages);
            Assert.AreEqual(updatedParameters.DigestionParams.MinPeptideLength, basicFileSpecificParameters.MinPeptideLength);
            Assert.AreEqual(updatedParameters.DigestionParams.MaxPeptideLength, basicFileSpecificParameters.MaxPeptideLength);
            Assert.AreEqual(updatedParameters.DigestionParams.Protease, basicFileSpecificParameters.Protease);

            Assert.AreEqual(updatedParameters.DoPrecursorDeconvolution, notDefaultParameters.DoPrecursorDeconvolution);
            Assert.AreEqual(updatedParameters.UseProvidedPrecursorInfo, notDefaultParameters.UseProvidedPrecursorInfo);
            Assert.AreEqual(updatedParameters.DeconvolutionIntensityRatio, notDefaultParameters.DeconvolutionIntensityRatio);
            Assert.AreEqual(updatedParameters.DeconvolutionMaxAssumedChargeState, notDefaultParameters.DeconvolutionMaxAssumedChargeState);
            Assert.AreEqual(updatedParameters.ReportAllAmbiguity, notDefaultParameters.ReportAllAmbiguity);
            Assert.AreEqual(updatedParameters.AddCompIons, notDefaultParameters.AddCompIons);
            Assert.AreEqual(updatedParameters.TotalPartitions, notDefaultParameters.TotalPartitions);
            Assert.AreEqual(updatedParameters.ScoreCutoff, notDefaultParameters.ScoreCutoff);
            Assert.AreEqual(updatedParameters.TopNpeaks, notDefaultParameters.TopNpeaks);
            Assert.AreEqual(updatedParameters.MinRatio, notDefaultParameters.MinRatio);
            Assert.AreEqual(updatedParameters.TrimMs1Peaks, notDefaultParameters.TrimMs1Peaks);
            Assert.AreEqual(updatedParameters.TrimMsMsPeaks, notDefaultParameters.TrimMsMsPeaks);
            Assert.AreEqual(updatedParameters.UseDeltaScore, notDefaultParameters.UseDeltaScore);
            Assert.AreEqual(updatedParameters.CalculateEValue, notDefaultParameters.CalculateEValue);
            Assert.AreEqual(updatedParameters.DeconvolutionMassTolerance, notDefaultParameters.DeconvolutionMassTolerance);
            Assert.AreEqual(updatedParameters.MaxThreadsToUsePerFile, notDefaultParameters.MaxThreadsToUsePerFile);
            Assert.AreEqual(updatedParameters.DigestionParams.InitiatorMethionineBehavior, notDefaultParameters.DigestionParams.InitiatorMethionineBehavior);
            Assert.AreEqual(updatedParameters.ListOfModsVariable, notDefaultParameters.ListOfModsVariable);
            Assert.AreEqual(updatedParameters.ListOfModsFixed, notDefaultParameters.ListOfModsFixed);

        }
    }
}