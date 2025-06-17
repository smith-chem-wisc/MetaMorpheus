using System;
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
            var model1 = new MassDifferenceAcceptorTypeModel
            {
                Type = MassDiffAcceptorType.OneMM,
                Label = "A",
                ToolTip = "B"
            };
            var model2 = new MassDifferenceAcceptorTypeModel
            {
                Type = MassDiffAcceptorType.OneMM,
                Label = "C",
                ToolTip = "D"
            };
            var model3 = new MassDifferenceAcceptorTypeModel
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
        public void LegacyParsing_Sets_Properties_Correctly()
        {
            var vm = new MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType.Custom, "LegacyName dot 0.5 ppm 0,1.0029");
            Assert.That(vm.CustomMode, Is.EqualTo(CustomMdacMode.Notch));
            Assert.That(vm.CustomName, Is.EqualTo("LegacyName"));
            Assert.That(vm.ToleranceValue, Is.EqualTo("0.5"));
            Assert.That(vm.SelectedToleranceType, Is.EqualTo("ppm"));
            Assert.That(vm.DotMassShifts, Is.EqualTo("0,1.0029"));
        }
    }
}
