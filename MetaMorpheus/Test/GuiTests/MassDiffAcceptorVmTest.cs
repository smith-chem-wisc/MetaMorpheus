using System;
using System.Linq;
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
            foreach (MassDiffAcceptorType type in Enum.GetValues<MassDiffAcceptorType>())
            {
                string customText = "custom";
                var vm = new MassDifferenceAcceptorSelectionViewModel(type, customText);

                Assert.That(vm.MassDiffAcceptorTypes.Select(m => m.Type), Is.EquivalentTo(Enum.GetValues<MassDiffAcceptorType>()));
                Assert.That(vm.SelectedType.Type, Is.EqualTo(type));
                Assert.That(vm.CustomMdac, Is.EqualTo(customText));
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
            vm.CustomMdac = "abc123";
            // Switch away from Custom, should cache and clear
            vm.SelectedType = exactModel;
            Assert.That(vm.CustomMdac, Is.EqualTo(string.Empty));
            // Switch back to Custom, should restore
            vm.SelectedType = customModel;
            Assert.That(vm.CustomMdac, Is.EqualTo("abc123"));
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
        }
    }
}
