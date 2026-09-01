using System;
using System.Linq;
using GuiFunctions;
using NUnit.Framework;
using SpectralAveraging;

namespace Test
{
    [TestFixture]
    public class AveragingGuiComponentsTest
    {
        [Test]
        public static void TestSpecAveragingParametersConstructor()
        {
            SpectralAveragingParameters specOptions = new SpectralAveragingParameters()
            { Percentile = 2, MinSigmaValue = 1.3, MaxSigmaValue = 3 };
            SpectralAveragingParametersViewModel parametersVm = new(specOptions);
            Assert.That(parametersVm.RejectionType == OutlierRejectionType.NoRejection);
            Assert.That(parametersVm.WeightingType == SpectraWeightingType.WeightEvenly);
            Assert.That(parametersVm.SpectraFileAveragingType == SpectraFileAveragingType.AverageAll);
            Assert.That(parametersVm.PerformNormalization);
            Assert.That(parametersVm.NumberOfScansToAverage == 5);
            Assert.That(Math.Abs(parametersVm.BinSize - 0.01) < 0.0001);
            Assert.That(Math.Abs(parametersVm.MinSigmaValue - 1.3) < 0.0001);
            Assert.That(Math.Abs(parametersVm.MaxSigmaValue - 3) < 0.0001);
            Assert.That(Math.Abs(parametersVm.Percentile - 2) < 0.0001);
            Assert.That(parametersVm.RejectionTypes.SequenceEqual(Enum.GetValues<OutlierRejectionType>()));
            Assert.That(parametersVm.WeightingTypes.Length == 2);
            Assert.That(parametersVm.SpectraFileAveragingTypes.Length == 3);
        }

        [Test]
        public static void TestSpecAveragingParametersResetDefaults()
        {
            // set all values
            SpectralAveragingParameters specOptions = new SpectralAveragingParameters()
            { Percentile = 2, MinSigmaValue = 1.3, MaxSigmaValue = 3 };
            SpectralAveragingParametersViewModel alteredParamsVm = new(specOptions)
            {
                RejectionType = OutlierRejectionType.AveragedSigmaClipping,
                WeightingType = SpectraWeightingType.MrsNoiseEstimation,
                SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
                PerformNormalization = false,
                Percentile = 14,
                MinSigmaValue = 20,
                MaxSigmaValue = 15,
                BinSize = 4,
                NumberOfScansToAverage = 20,
            };
            SpectralAveragingParametersViewModel defaultParamsVm = new(new SpectralAveragingParameters());

            // ensure they were set in the modified version
            Assert.That(Math.Abs(alteredParamsVm.Percentile - 14) < 0.0001);
            Assert.That(alteredParamsVm.RejectionType == OutlierRejectionType.AveragedSigmaClipping);
            Assert.That(alteredParamsVm.WeightingType == SpectraWeightingType.MrsNoiseEstimation);
            Assert.That(alteredParamsVm.SpectraFileAveragingType == SpectraFileAveragingType.AverageDdaScans);
            Assert.That(!alteredParamsVm.PerformNormalization);
            Assert.That(alteredParamsVm.NumberOfScansToAverage == 20);
            Assert.That(Math.Abs(alteredParamsVm.BinSize - 4) < 0.0001);
            Assert.That(Math.Abs(alteredParamsVm.MinSigmaValue - 20) < 0.0001);
            Assert.That(Math.Abs(alteredParamsVm.MaxSigmaValue - 15) < 0.0001);

            // reset all values 
            Assert.That(!alteredParamsVm.Equals(defaultParamsVm));
            alteredParamsVm.ResetDefaults();
            Assert.That(alteredParamsVm.Equals(defaultParamsVm));

            // ensure parameters were reset
            foreach (var property in alteredParamsVm.SpectralAveragingParameters.GetType().GetProperties())
            {
                Assert.That(property.GetValue(alteredParamsVm.SpectralAveragingParameters), 
                    Is.EqualTo(property.GetValue(defaultParamsVm.SpectralAveragingParameters)));
            }

            // ensure view model was reset
            foreach (var property in defaultParamsVm.GetType().GetProperties().Where(p => p.PropertyType != typeof(SpectralAveragingParameters)))
            {
                Assert.That(property.GetValue(alteredParamsVm), Is.EqualTo(property.GetValue(defaultParamsVm)));
            }
        }

        [Test]
        public static void TestSpectralAveragingOptionsSpecialGettersAndSetters()
        {
            SpectralAveragingParameters specOptions = new SpectralAveragingParameters()
            { Percentile = 2, MinSigmaValue = 1.3, MaxSigmaValue = 3 };
            SpectralAveragingParametersViewModel parametersVm = new(specOptions);

            parametersVm.SpectraFileAveragingType = SpectraFileAveragingType.AverageEverynScansWithOverlap;
            Assert.That(specOptions.SpectraFileAveragingType == SpectraFileAveragingType.AverageEverynScansWithOverlap);
            parametersVm.NumberOfScansToAverage = 3;
            Assert.That(parametersVm.NumberOfScansToAverage == 3);
            Assert.That(specOptions.NumberOfScansToAverage == 3);
            Assert.That(specOptions.ScanOverlap == 2);

            parametersVm.NumberOfScansToAverage = 5;
            Assert.That(parametersVm.NumberOfScansToAverage == 5);
            Assert.That(specOptions.NumberOfScansToAverage == 5);
            Assert.That(specOptions.ScanOverlap == 4);

            parametersVm.SpectraFileAveragingType = SpectraFileAveragingType.AverageAll;
            Assert.That(parametersVm.NumberOfScansToAverage == 5);
            Assert.That(specOptions.NumberOfScansToAverage == 5);
            Assert.That(specOptions.ScanOverlap == 4);

            parametersVm.SpectraFileAveragingType = SpectraFileAveragingType.AverageEverynScansWithOverlap;
            Assert.That(specOptions.SpectraFileAveragingType == SpectraFileAveragingType.AverageEverynScansWithOverlap);
            Assert.That(parametersVm.NumberOfScansToAverage == 5);
            Assert.That(specOptions.NumberOfScansToAverage == 5);
            Assert.That(specOptions.ScanOverlap == 4);


            parametersVm.PerformNormalization = true;
            Assert.That(parametersVm.PerformNormalization);
            Assert.That(specOptions.NormalizationType == NormalizationType.RelativeToTics);

            parametersVm.PerformNormalization = false;
            Assert.That(!parametersVm.PerformNormalization);
            Assert.That(specOptions.NormalizationType == NormalizationType.NoNormalization);
        }


        [Test]
        public static void TestViewModelSetOtherParametersAndEquality()
        {
            var viewModel = new SpectralAveragingParametersViewModel(new SpectralAveragingParameters());

            var presetDda1 = new SpectralAveragingParameters()
            {
                OutputType = OutputType.MzML,
                NormalizationType = NormalizationType.RelativeToTics,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly,
                BinSize = 0.01,
                SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
                NumberOfScansToAverage = 5,
                ScanOverlap = 4,
                MaxSigmaValue = 3,
                MinSigmaValue = 0.5,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping
            };

            var presetDda1Vm = new SpectralAveragingParametersViewModel(presetDda1);
            var defaultParamsVm = new SpectralAveragingParametersViewModel(new SpectralAveragingParameters());

            Assert.That(!presetDda1Vm.Equals(defaultParamsVm));
            defaultParamsVm.SetOtherParameters("Dda1");
            Assert.That(presetDda1Vm.Equals(defaultParamsVm));

            var presetDda2 = new SpectralAveragingParameters()
            {
                OutputType = OutputType.MzML,
                NormalizationType = NormalizationType.RelativeToTics,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly,
                BinSize = 0.01,
                SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
                NumberOfScansToAverage = 5,
                ScanOverlap = 4,
                MaxSigmaValue = 3,
                MinSigmaValue = 0.5,
                OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping
            };

            var presetDda2Vm = new SpectralAveragingParametersViewModel(presetDda2);
            defaultParamsVm = new SpectralAveragingParametersViewModel(new SpectralAveragingParameters());

            Assert.That(!presetDda2Vm.Equals(defaultParamsVm));
            defaultParamsVm.SetOtherParameters("Dda2");
            Assert.That(presetDda2Vm.Equals(defaultParamsVm));

            var presetDirectInjection = new SpectralAveragingParameters()
            {
                OutputType = OutputType.MzML,
                NormalizationType = NormalizationType.RelativeToTics,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly,
                BinSize = 0.01,
                SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
                NumberOfScansToAverage = 15,
                ScanOverlap = 14,
                OutlierRejectionType = OutlierRejectionType.MinMaxClipping
            };

            var presetDirectInjectionVm = new SpectralAveragingParametersViewModel(presetDirectInjection);
            defaultParamsVm = new SpectralAveragingParametersViewModel(new SpectralAveragingParameters());

            Assert.That(!presetDirectInjectionVm.Equals(defaultParamsVm));
            defaultParamsVm.SetOtherParameters("DirectInjection");
            Assert.That(presetDirectInjectionVm.Equals(defaultParamsVm));


            Assert.Throws<ArgumentException> (() => defaultParamsVm.SetOtherParameters("-1"));
        }


        [Test]
        public static void TestViewModelEqualityAllProperties()
        {
            var viewModel = new SpectralAveragingParametersViewModel(new SpectralAveragingParameters());

            var numericalTypes = new[] { typeof(int), typeof(double), typeof(float), typeof(OutputType),
                typeof(OutlierRejectionType), typeof(SpectraWeightingType), typeof(SpectraFileAveragingType),
                typeof(NormalizationType), typeof(SpectralAveragingType), };

            // check properties in view model
            foreach (var property in viewModel.GetType().GetProperties())
            {
                var propertyType = property.PropertyType;
                var testViewModel = new SpectralAveragingParametersViewModel(new SpectralAveragingParameters());

                Assert.That(viewModel.Equals(testViewModel));
                if (numericalTypes.Any(p => p == propertyType))
                {
                    property.SetValue(testViewModel, -1, null);
                    Assert.That(!viewModel.Equals(testViewModel));
                }
            }

            // check properties in the parameters the view model wraps
            foreach (var property in viewModel.SpectralAveragingParameters.GetType().GetProperties())
            {
                var propertyType = property.PropertyType;
                var testViewModel = new SpectralAveragingParametersViewModel(new SpectralAveragingParameters());

                Assert.That(viewModel.Equals(testViewModel));
                if (numericalTypes.Any(p => p == propertyType))
                {
                    property.SetValue(testViewModel.SpectralAveragingParameters, -1, null);
                    Assert.That(!viewModel.Equals(testViewModel));
                }
            }
        }

        
    }
}
