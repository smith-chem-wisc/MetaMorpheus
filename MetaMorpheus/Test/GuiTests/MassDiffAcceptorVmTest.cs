using System;
using System.Globalization;
using System.Linq;
using System.Reflection;
using GuiFunctions;
using NUnit.Framework;
using TaskLayer;

namespace Test.GuiTests
{
    [TestFixture]
    public class MassDifferenceAcceptorViewModelTests
    {
        [Test]
        public void Constructor_InitializesProperties_Correctly()
        {
            string defaultCustomText = "Custom dot 5 ppm 0";
            foreach (MassDiffAcceptorType type in Enum.GetValues<MassDiffAcceptorType>())
            {
                var vm = new MassDifferenceAcceptorSelectionViewModel(type, "");

                Assert.That(vm.MassDiffAcceptorTypes.Select(m => m.Type), Is.EquivalentTo(Enum.GetValues<MassDiffAcceptorType>()));
                Assert.That(vm.SelectedType.Type, Is.EqualTo(type));

                if (type == MassDiffAcceptorType.Custom)
                {
                    Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.Notch)); // Default mode
                    Assert.That(vm.CustomName, Is.EqualTo("Custom")); // Default name
                    Assert.That(vm.CustomMdac, Is.EqualTo(defaultCustomText));
                }
                else
                {
                    Assert.That(vm.CustomMdac, Is.EqualTo(string.Empty));
                }
            }
        }

        [Test]
        public void SelectedType_CachesAndRestoresCustomMdac()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Exact, "");
            var customModel = vm.MassDiffAcceptorTypes.First(m => m.Type == MassDiffAcceptorType.Custom);
            var exactModel = vm.MassDiffAcceptorTypes.First(m => m.Type == MassDiffAcceptorType.Exact);

            // Switch to Custom, set CustomMdac
            vm.SelectedType = customModel;
            vm.CustomMdac = "TestNotch dot 0.1 da 0,1.0029,2.0052";
            // Switch away from Custom, should cache and clear
            vm.SelectedType = exactModel;
            Assert.That(vm.CustomMdac, Is.EqualTo(string.Empty));
            // Switch back to Custom, should restore
            vm.SelectedType = customModel;
            Assert.That(vm.CustomMdac, Is.EqualTo("TestNotch dot 0.1 da 0,1.0029,2.0052"));
        }

        [Test]
        public void SelectedType_CrashesOnUnknownType()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Exact, "");
            // Create a fake/invalid enum value not present in MassDiffAcceptorType
            var invalidType = (MassDiffAcceptorType)9999;

            // Use reflection to call the private CreateModel method with the invalid type
            var method = typeof(MassDifferenceAcceptorSelectionViewModel)
                .GetMethod("CreateModel", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            try
            {
                method.Invoke(vm, new object[] { invalidType });
                Assert.Fail("Expected NotImplementedException was not thrown.");
            }
            catch (TargetInvocationException ex)
            {
                Assert.That(ex.InnerException, Is.TypeOf<NotImplementedException>());
                Assert.That(ex.InnerException.Message, Is.EqualTo("No model implemented for type: 9999"));
            }
        }

        [Test]
        public void CustomMdac_Property_SetAndGet()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            bool propertyChangedFired = false;
            vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.CustomMdac)) propertyChangedFired = true; };

            vm.CustomMdac = "test";
            Assert.That(vm.CustomMdac, Is.EqualTo("test"));
            Assert.That(propertyChangedFired, Is.True);
        }

        [Test]
        public void MassDifferenceAcceptorTypeModel_Equality_Works()
        {
            var model1 = new MassDifferenceAcceptorTypeViewModel
            {
                Type = MassDiffAcceptorType.OneMM,
                Label = "A",
                ToolTip = "B"
            };
            var model2 = new MassDifferenceAcceptorTypeViewModel
            {
                Type = MassDiffAcceptorType.OneMM,
                Label = "C",
                ToolTip = "D"
            };
            var model3 = new MassDifferenceAcceptorTypeViewModel
            {
                Type = MassDiffAcceptorType.TwoMM,
                Label = "E",
                ToolTip = "F"
            };

            Assert.That(model1.ToolTip, Is.EqualTo("B"));
            Assert.That(model1.Label, Is.EqualTo("A"));
            Assert.That(model1.Equals(model2), Is.True);
            Assert.That(model1.Equals(model3), Is.False);
            Assert.That(model1.Equals(MassDiffAcceptorType.OneMM), Is.True);
            Assert.That(model1.Equals(MassDiffAcceptorType.TwoMM), Is.False);
            Assert.That(model1.Equals((object)model2), Is.True);
            Assert.That(model1.Equals((object)null), Is.False);
            Assert.That(model1.Equals((object)"not a model"), Is.False);
            Assert.That(model1.GetHashCode(), Is.EqualTo((int)MassDiffAcceptorType.OneMM));
        }

        [Test]
        public void MassDifferenceAcceptorTypeModel_Equality_MissedMonosIsCorrect()
        {
            var host = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.OneMM, "");
            var model1 = host.MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.OneMM);

            Assert.That(model1.PositiveMissedMonos, Is.EqualTo(1));
            Assert.That(model1.NegativeMissedMonos, Is.EqualTo(0));

            var model2 = host.MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.TwoMM);
            Assert.That(model2.PositiveMissedMonos, Is.EqualTo(2));
            Assert.That(model2.NegativeMissedMonos, Is.EqualTo(0));

            var model3 = host.MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.PlusOrMinusThreeMM);
            Assert.That(model3.PositiveMissedMonos, Is.EqualTo(3));
            Assert.That(model3.NegativeMissedMonos, Is.EqualTo(3));

            var model4 = host.MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.ModOpen);
            Assert.That(model4.PositiveMissedMonos, Is.EqualTo(0));
            Assert.That(model4.NegativeMissedMonos, Is.EqualTo(0));

            var model5 = host.MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.Open);
            Assert.That(model5.PositiveMissedMonos, Is.EqualTo(0));
            Assert.That(model5.NegativeMissedMonos, Is.EqualTo(0));

            var model6 = host.MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.Custom);
            Assert.That(model6.PositiveMissedMonos, Is.EqualTo(0));
            Assert.That(model6.NegativeMissedMonos, Is.EqualTo(0));
        }

        [Test]
        public void MassDifferenceAcceptorModel_Instance_ReturnsDefault()
        {
            var instance = MassDifferenceAcceptorSelectionModel.Instance;
            Assert.That(instance, Is.TypeOf<MassDifferenceAcceptorSelectionModel>());
            Assert.That(instance.SelectedType.Type, Is.EqualTo(MassDiffAcceptorType.TwoMM));
            Assert.That(instance.CustomMdac, Is.EqualTo(""));
            Assert.That(instance.CustomMdacModes.Count, Is.EqualTo(3)); 
        }

        [Test]
        public void NotchMode_Properties_Produce_Valid_CustomMdac()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            vm.CustomMode = CustomMdacMode.Notch;
            vm.CustomName = "TestNotch";
            vm.ToleranceValue = "0.1";
            vm.SelectedToleranceType = "da";
            vm.DotMassShifts = "0,1.0029,2.0052";

            string expected = "TestNotch dot 0.1 da 0,1.0029,2.0052";
            Assert.That(vm.CustomMdac, Is.EqualTo(expected));
        }

        [Test]
        public void IntervalMode_Properties_Produce_Valid_CustomMdac()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            vm.CustomMode = CustomMdacMode.Interval;
            vm.CustomName = "TestInterval";
            vm.IntervalRanges = "[0,200];[300,400]";

            string expected = "TestInterval interval [0,200];[300,400]";
            Assert.That(vm.CustomMdac, Is.EqualTo(expected));
        }

        [Test]
        public void AroundZeroMode_Properties_Produce_Valid_CustomMdac_Ppm()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            vm.CustomMode = CustomMdacMode.AroundZero;
            vm.CustomName = "TestAroundZero";
            vm.SelectedToleranceType = "ppm";
            vm.ToleranceValue = "5";

            string expected = "TestAroundZero ppmAroundZero 5";
            Assert.That(vm.CustomMdac, Is.EqualTo(expected));
        }

        [Test]
        public void AroundZeroMode_Properties_Produce_Valid_CustomMdac_Da()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            vm.CustomMode = CustomMdacMode.AroundZero;
            vm.CustomName = "TestAroundZero";
            vm.SelectedToleranceType = "da";
            vm.ToleranceValue = "2.1";

            string expected = "TestAroundZero daAroundZero 2.1";
            Assert.That(vm.CustomMdac, Is.EqualTo(expected));
        }

        [Test]
        public void ManualCustomMdac_Overrides_GuiFields()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            string manual = "ManualName dot 0.2 ppm 0,1.003";
            vm.CustomMdac = manual;
            Assert.That(vm.CustomMdac, Is.EqualTo(manual));
            // Changing a GUI field should switch back to auto-generation
            vm.CustomMode = CustomMdacMode.Notch;
            vm.CustomName = "AutoName";
            vm.ToleranceValue = "0.3";
            vm.SelectedToleranceType = "da";
            vm.DotMassShifts = "0,1.0029";
            string expected = "AutoName dot 0.3 da 0,1.0029";
            Assert.That(vm.CustomMdac, Is.EqualTo(expected));
        }

        [Test]
        public void LegacyParsing_Notch_SetsPropertiesCorrectly()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "LegacyName dot 0.5 ppm 0,1.0033548");
            Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.Notch));
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Custom));
            Assert.That(vm.CustomName, Is.EqualTo("LegacyName"));
            Assert.That(vm.ToleranceValue, Is.EqualTo("0.5"));
            Assert.That(vm.SelectedToleranceType, Is.EqualTo("ppm"));
            Assert.That(vm.DotMassShifts, Is.EqualTo("0"));

            var mm = vm.PredefinedNotches.FirstOrDefault(p => p.Name == "Missed Mono");
            Assert.That(mm.IsSelected, Is.True);
            Assert.That(mm.MaxPositiveFrequency, Is.EqualTo(1));
            Assert.That(mm.MaxNegativeFrequency, Is.EqualTo(0));
        }

        [Test]
        public void LegacyParsing_Notch_SetsSelectableNotchesCorrectly()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "Custom dot 7 ppm -2.00671,-1.00335,0,1.00335,2.00671,3.01006");
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Custom));
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Custom));
            Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.Notch));
            Assert.That(vm.CustomName, Is.EqualTo("Custom"));
            Assert.That(vm.ToleranceValue, Is.EqualTo("7"));
            Assert.That(vm.SelectedToleranceType, Is.EqualTo("ppm"));
            Assert.That(vm.DotMassShifts, Is.EqualTo("0"));

            var mm = vm.PredefinedNotches.FirstOrDefault(p => p.Name == "Missed Mono");
            Assert.That(mm.IsSelected, Is.True);
            Assert.That(mm.MaxPositiveFrequency, Is.EqualTo(3));
            Assert.That(mm.MaxNegativeFrequency, Is.EqualTo(2));
        }

        [Test]
        public void LegacyParsing_Interval_SetsPropertiesCorrectly()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "LegacyName interval [1,140];[18,20]");
            Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.Interval));
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Custom));
            Assert.That(vm.CustomName, Is.EqualTo("LegacyName"));
            Assert.That(vm.DotMassShifts, Is.EqualTo("0"));

            Assert.That(vm.IntervalRanges, Is.EqualTo("[1,140];[18,20]"));
        }

        [Test]
        public void LegacyParsing_AroundZero_PPM_SetsPropertiesCorrectly()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "LegacyName ppmAroundZero 7");
            Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.AroundZero));
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Custom));
            Assert.That(vm.CustomName, Is.EqualTo("LegacyName"));
            Assert.That(vm.ToleranceValue, Is.EqualTo("7"));
            Assert.That(vm.SelectedToleranceType, Is.EqualTo("ppm"));
            Assert.That(vm.DotMassShifts, Is.EqualTo("0"));
        }

        [Test]
        public void LegacyParsing_AroundZero_Da_SetsPropertiesCorrectly()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "LegacyName daAroundZero 17");
            Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.AroundZero));
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Custom));
            Assert.That(vm.CustomName, Is.EqualTo("LegacyName"));
            Assert.That(vm.ToleranceValue, Is.EqualTo("17"));
            Assert.That(vm.SelectedToleranceType, Is.EqualTo("da"));
            Assert.That(vm.DotMassShifts, Is.EqualTo("0"));

            vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "LegacyName daltonsAroundZero 17");
            Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.AroundZero));
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Custom));
            Assert.That(vm.CustomName, Is.EqualTo("LegacyName"));
            Assert.That(vm.ToleranceValue, Is.EqualTo("17"));
            Assert.That(vm.SelectedToleranceType, Is.EqualTo("da"));
            Assert.That(vm.DotMassShifts, Is.EqualTo("0"));
        }

        [Test]
        public void LegacyParsing_Open_SetsPropertiesCorrectly()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "LegacyName OpenSearch");
            Assert.That(vm.SelectedType, Is.EqualTo(MassDiffAcceptorType.Open));
            Assert.That(vm.CustomName, Is.EqualTo("LegacyName"));
        }

        [Test]
        public void PredefinedNotches_Defaults_AreCorrect()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            Assert.That(vm.PredefinedNotches.Count, Is.EqualTo(1));
            var notch = vm.PredefinedNotches[0];
            Assert.That(notch.Name, Is.EqualTo("Missed Mono"));
            Assert.That(notch.MonoisotopicMass, Is.EqualTo(Chemistry.Constants.C13MinusC12));
            Assert.That(notch.MaxPositiveFrequency, Is.EqualTo(0));
            Assert.That(notch.MaxNegativeFrequency, Is.EqualTo(0));
            Assert.That(notch.IsSelected, Is.False);
        }

        [Test]
        public void SelectableNotch_PropertyChanged_TriggersCustomMdacUpdate()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            var notch = vm.PredefinedNotches[0];
            bool customMdacChanged = false;
            vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.CustomMdac)) customMdacChanged = true; };

            notch.IsSelected = true;
            Assert.That(customMdacChanged, Is.True);
            customMdacChanged = false;
            notch.MaxPositiveFrequency = 2;
            Assert.That(customMdacChanged, Is.True);
            customMdacChanged = false;
            notch.MaxNegativeFrequency = 1;
            Assert.That(customMdacChanged, Is.True);
        }

        [Test]
        public void GetAllNotchCombinations_ProducesAllSums()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            // Add a second notch for more complex combinations
            vm.PredefinedNotches.Add(new SelectableNotchViewModel("Other", 2.0));
            var missedMono = vm.PredefinedNotches[0];
            var other = vm.PredefinedNotches[1];

            missedMono.IsSelected = true;
            missedMono.MaxPositiveFrequency = 1;
            missedMono.MaxNegativeFrequency = 2;
            other.IsSelected = true;
            other.MaxPositiveFrequency = 1;
            other.MaxNegativeFrequency = 1;

            // Should produce all combinations of -2, -1, 0, 1 (missedMono) and -1*2.0, 0, 1*2.0 (other)
            var combos = vm.GetType()
                .GetMethod("GetAllNotchCombinations", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(vm, null) as System.Collections.IEnumerable;

            var results = combos.Cast<double>().OrderBy(x => x).ToList();

            // For missedMono: -2, -1, 0, 1; for other: -2, 0, 2
            // All possible sums:
            var expected = (
                from mm in new[] { -2, -1, 0, 1 }.Select(i => i * missedMono.MonoisotopicMass)
                from o in new[] { -2.0, 0.0, 2.0 }
                select mm + o
            ).Distinct().OrderBy(x => x).ToList();

            Assert.That(results, Is.EquivalentTo(expected));
        }

        [Test]
        public void BuildCustomMdac_CombinesNotchCombinationsAndUserInput()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            vm.CustomMode = CustomMdacMode.Notch;
            vm.CustomName = "ComboTest";
            vm.ToleranceValue = "0.1";
            vm.SelectedToleranceType = "da";
            // Add a second notch
            vm.PredefinedNotches.Add(new SelectableNotchViewModel("Other", 2.0));
            var missedMono = vm.PredefinedNotches[0];
            var other = vm.PredefinedNotches[1];

            missedMono.IsSelected = true;
            missedMono.MaxPositiveFrequency = 1;
            missedMono.MaxNegativeFrequency = 2;
            other.IsSelected = true;
            other.MaxPositiveFrequency = 1;
            other.MaxNegativeFrequency = 1;

            // User also enters a custom value
            vm.DotMassShifts = "3.5";

            // Build expected set
            var mmVals = new[] { -2, -1, 0, 1 }.Select(i => i * missedMono.MonoisotopicMass);
            var oVals = new[] { -2.0, 0.0, 2.0 };
            var combos = (
                from mm in mmVals
                from o in oVals
                select mm + o
            ).Concat(new[] { 3.5 }).Distinct().OrderBy(x => x).ToList();

            var expected = $"ComboTest dot 0.1 da {string.Join(",", combos.Select(x => x.ToString("G6", CultureInfo.InvariantCulture)))}";
            Assert.That(vm.CustomMdac, Is.EqualTo(expected));
        }

        [Test]
        public void GetAllNotchCombinations_ReturnsZero_WhenNoneSelected()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            foreach (var n in vm.PredefinedNotches)
            {
                n.IsSelected = false;
            }
            var combos = vm.GetType()
                .GetMethod("GetAllNotchCombinations", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(vm, null) as System.Collections.IEnumerable;
            var results = combos.Cast<double>().ToList();
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0], Is.EqualTo(0));
        }

        [Test]
        public void CartesianProduct_WorksForMultipleRanges()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "");
            var method = vm.GetType().GetMethod("CartesianProduct", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);

            int[][] ranges = new[] {
                new[] { -1, 0, 1 },
                new[] { 0, 2 }
            };
            var result = method.Invoke(null, new object[] { ranges }) as System.Collections.IEnumerable;
            var combos = result.Cast<int[]>().ToList();

            // Should be 3*2 = 6 combinations
            Assert.That(combos.Count, Is.EqualTo(6));
            Assert.That(combos, Does.Contain(new[] { -1, 0 }));
            Assert.That(combos, Does.Contain(new[] { 1, 2 }));
        }
    }
}
