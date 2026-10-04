using System;
using System.Collections.ObjectModel;
using System.ComponentModel;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using NUnit.Framework;

namespace Test.MetaDraw
{
    [TestFixture]
    public class PlotModelStatParametersTests
    {
        #region PlotModelStatParameters Tests

        [Test]
        public void TestDefaultPropertyValues()
        {
            // Arrange & Act
            var parameters = new PlotModelStatParameters();

            // Assert - verify default values
            Assert.That(parameters.UseLogScaleYAxis, Is.EqualTo(false));
            Assert.That(parameters.GroupingProperty, Is.EqualTo("None"));
            Assert.That(parameters.AllowAmbiguousGroups, Is.EqualTo(false));
            Assert.That(parameters.MinRelativeCutoff, Is.EqualTo(0.0));
            Assert.That(parameters.MaxRelativeCutoff, Is.EqualTo(100.0));
            Assert.That(parameters.NormalizeHistogramToFile, Is.EqualTo(false));
            Assert.That(parameters.DisplayFilteredOnly, Is.EqualTo(true));
        }

        [Test]
        public void TestClone_CreatesIndependentCopy()
        {
            // Arrange
            var original = new PlotModelStatParameters
            {
                UseLogScaleYAxis = true,
                GroupingProperty = "Precursor Charge",
                AllowAmbiguousGroups = true,
                MinRelativeCutoff = 5.5,
                MaxRelativeCutoff = 95.5,
                NormalizeHistogramToFile = true,
                DisplayFilteredOnly = false
            };

            // Act
            var clone = original.Clone();

            // Assert - values are copied
            Assert.That(clone.UseLogScaleYAxis, Is.EqualTo(original.UseLogScaleYAxis));
            Assert.That(clone.GroupingProperty, Is.EqualTo(original.GroupingProperty));
            Assert.That(clone.AllowAmbiguousGroups, Is.EqualTo(original.AllowAmbiguousGroups));
            Assert.That(clone.MinRelativeCutoff, Is.EqualTo(original.MinRelativeCutoff));
            Assert.That(clone.MaxRelativeCutoff, Is.EqualTo(original.MaxRelativeCutoff));
            Assert.That(clone.NormalizeHistogramToFile, Is.EqualTo(original.NormalizeHistogramToFile));
            Assert.That(clone.DisplayFilteredOnly, Is.EqualTo(original.DisplayFilteredOnly));

            // Assert - clone is independent (modifying clone doesn't affect original)
            clone.UseLogScaleYAxis = !original.UseLogScaleYAxis;
            clone.GroupingProperty = "File Name";
            clone.AllowAmbiguousGroups = !original.AllowAmbiguousGroups;
            clone.MinRelativeCutoff = 10.0;
            clone.MaxRelativeCutoff = 90.0;
            clone.NormalizeHistogramToFile = !original.NormalizeHistogramToFile;
            clone.DisplayFilteredOnly = !original.DisplayFilteredOnly;

            Assert.That(original.UseLogScaleYAxis, Is.Not.EqualTo(clone.UseLogScaleYAxis));
            Assert.That(original.GroupingProperty, Is.Not.EqualTo(clone.GroupingProperty));
            Assert.That(original.AllowAmbiguousGroups, Is.Not.EqualTo(clone.AllowAmbiguousGroups));
            Assert.That(original.MinRelativeCutoff, Is.Not.EqualTo(clone.MinRelativeCutoff));
            Assert.That(original.MaxRelativeCutoff, Is.Not.EqualTo(clone.MaxRelativeCutoff));
            Assert.That(original.NormalizeHistogramToFile, Is.Not.EqualTo(clone.NormalizeHistogramToFile));
            Assert.That(original.DisplayFilteredOnly, Is.Not.EqualTo(clone.DisplayFilteredOnly));
        }

        [Test]
        public void TestEquals_NullComparison()
        {
            // Arrange
            var parameters = new PlotModelStatParameters();

            // Assert
            Assert.That(parameters.Equals(null), Is.False);
        }

        [Test]
        public void TestEquals_SameValues()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters
            {
                UseLogScaleYAxis = true,
                GroupingProperty = "Precursor Charge",
                AllowAmbiguousGroups = true,
                MinRelativeCutoff = 5.5,
                MaxRelativeCutoff = 95.5,
                NormalizeHistogramToFile = true,
                DisplayFilteredOnly = false
            };

            var parameters2 = new PlotModelStatParameters
            {
                UseLogScaleYAxis = true,
                GroupingProperty = "Precursor Charge",
                AllowAmbiguousGroups = true,
                MinRelativeCutoff = 5.5,
                MaxRelativeCutoff = 95.5,
                NormalizeHistogramToFile = true,
                DisplayFilteredOnly = false
            };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.True);
        }

        [Test]
        public void TestEquals_DifferentUseLogScaleYAxis()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { UseLogScaleYAxis = true };
            var parameters2 = new PlotModelStatParameters { UseLogScaleYAxis = false };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        [Test]
        public void TestEquals_DifferentGroupingProperty()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { GroupingProperty = "None" };
            var parameters2 = new PlotModelStatParameters { GroupingProperty = "Precursor Charge" };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        [Test]
        public void TestEquals_DifferentAllowAmbiguousGroups()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { AllowAmbiguousGroups = true };
            var parameters2 = new PlotModelStatParameters { AllowAmbiguousGroups = false };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        [Test]
        public void TestEquals_DifferentMinRelativeCutoff()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { MinRelativeCutoff = 5.0 };
            var parameters2 = new PlotModelStatParameters { MinRelativeCutoff = 10.0 };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        [Test]
        public void TestEquals_DifferentMaxRelativeCutoff()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { MaxRelativeCutoff = 95.0 };
            var parameters2 = new PlotModelStatParameters { MaxRelativeCutoff = 90.0 };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        [Test]
        public void TestEquals_DifferentNormalizeHistogramToFile()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { NormalizeHistogramToFile = true };
            var parameters2 = new PlotModelStatParameters { NormalizeHistogramToFile = false };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        [Test]
        public void TestEquals_DifferentDisplayFilteredOnly()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { DisplayFilteredOnly = true };
            var parameters2 = new PlotModelStatParameters { DisplayFilteredOnly = false };

            // Assert
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        [Test]
        public void TestEquals_FloatingPointTolerance_MinRelativeCutoff()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { MinRelativeCutoff = 5.0 };
            var parameters2 = new PlotModelStatParameters { MinRelativeCutoff = 5.00005 };

            // Assert - values within tolerance should be equal
            Assert.That(parameters1.Equals(parameters2), Is.True);
        }

        [Test]
        public void TestEquals_FloatingPointTolerance_MaxRelativeCutoff()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { MaxRelativeCutoff = 95.0 };
            var parameters2 = new PlotModelStatParameters { MaxRelativeCutoff = 95.00005 };

            // Assert - values within tolerance should be equal
            Assert.That(parameters1.Equals(parameters2), Is.True);
        }

        [Test]
        public void TestEquals_FloatingPointTolerance_OutsideTolerance()
        {
            // Arrange
            var parameters1 = new PlotModelStatParameters { MinRelativeCutoff = 5.0 };
            var parameters2 = new PlotModelStatParameters { MinRelativeCutoff = 5.001 };

            // Assert - values outside tolerance should not be equal
            Assert.That(parameters1.Equals(parameters2), Is.False);
        }

        #endregion

        #region PlotModelStatParametersViewModel Tests

        [SetUp]
        public void SetUp()
        {
            // Reset settings to ensure clean state
            MetaDrawSettings.ResetSettings();
        }

        [TearDown]
        public void TearDown()
        {
            // Reset settings after each test
            MetaDrawSettings.ResetSettings();
        }

        [Test]
        public void TestSingletonInstance_ReturnsSameInstance()
        {
            // Act
            var instance1 = PlotModelStatParametersViewModel.Instance;
            var instance2 = PlotModelStatParametersViewModel.Instance;

            // Assert
            Assert.That(instance1, Is.SameAs(instance2));
        }

        [Test]
        public void TestGroupingProperties_InitializedCorrectly()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Assert
            Assert.That(viewModel.GroupingProperties, Is.Not.Null);
            Assert.That(viewModel.GroupingProperties.Count, Is.GreaterThan(0));
            Assert.That(viewModel.GroupingProperties, Contains.Item("None"));
            Assert.That(viewModel.GroupingProperties, Contains.Item("Precursor Charge"));
            Assert.That(viewModel.GroupingProperties, Contains.Item("File Name"));
        }

        [Test]
        public void TestUseLogScaleYAxis_GetterSetter()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act & Assert - default is false
            Assert.That(viewModel.UseLogScaleYAxis, Is.False);

            // Act
            viewModel.UseLogScaleYAxis = true;

            // Assert
            Assert.That(viewModel.UseLogScaleYAxis, Is.True);

            // Act
            viewModel.UseLogScaleYAxis = false;

            // Assert
            Assert.That(viewModel.UseLogScaleYAxis, Is.False);
        }

        [Test]
        public void TestGroupingProperty_GetterSetter()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act & Assert - default is "None"
            Assert.That(viewModel.GroupingProperty, Is.EqualTo("None"));

            // Act
            viewModel.GroupingProperty = "Precursor Charge";

            // Assert
            Assert.That(viewModel.GroupingProperty, Is.EqualTo("Precursor Charge"));

            // Act
            viewModel.GroupingProperty = "File Name";

            // Assert
            Assert.That(viewModel.GroupingProperty, Is.EqualTo("File Name"));
        }

        [Test]
        public void TestGroupingProperty_NullValue_DefaultsToNone()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            viewModel.GroupingProperty = "Precursor Charge";

            // Act
            viewModel.GroupingProperty = null;

            // Assert - null should default to "None"
            Assert.That(viewModel.GroupingProperty, Is.EqualTo("None"));
        }

        [Test]
        public void TestAllowAmbiguousGroups_GetterSetter()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act & Assert - default is false
            Assert.That(viewModel.AllowAmbiguousGroups, Is.False);

            // Act
            viewModel.AllowAmbiguousGroups = true;

            // Assert
            Assert.That(viewModel.AllowAmbiguousGroups, Is.True);

            // Act
            viewModel.AllowAmbiguousGroups = false;

            // Assert
            Assert.That(viewModel.AllowAmbiguousGroups, Is.False);
        }

        [Test]
        public void TestMinRelativeCutoff_GetterSetter_ValidRange()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act & Assert - default is 0.0
            Assert.That(viewModel.MinRelativeCutoff, Is.EqualTo(0.0));

            // Act
            viewModel.MinRelativeCutoff = 50.0;

            // Assert
            Assert.That(viewModel.MinRelativeCutoff, Is.EqualTo(50.0));

            // Act
            viewModel.MinRelativeCutoff = 100.0;

            // Assert
            Assert.That(viewModel.MinRelativeCutoff, Is.EqualTo(100.0));
        }

        [Test]
        public void TestMinRelativeCutoff_BelowZero_ClampedToZero()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act
            viewModel.MinRelativeCutoff = -10.0;

            // Assert - should be clamped to 0
            Assert.That(viewModel.MinRelativeCutoff, Is.EqualTo(0.0));
        }

        [Test]
        public void TestMinRelativeCutoff_Above100_ClampedTo100()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act
            viewModel.MinRelativeCutoff = 150.0;

            // Assert - should be clamped to 100
            Assert.That(viewModel.MinRelativeCutoff, Is.EqualTo(100.0));
        }

        [Test]
        public void TestMinRelativeCutoff_ConstraintVsMaxRelativeCutoff()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            viewModel.MaxRelativeCutoff = 30.0;

            // Act - set MinRelativeCutoff above MaxRelativeCutoff
            viewModel.MinRelativeCutoff = 50.0;

            // Assert - should be clamped to MaxRelativeCutoff
            Assert.That(viewModel.MinRelativeCutoff, Is.EqualTo(30.0));
        }

        [Test]
        public void TestMaxRelativeCutoff_GetterSetter_ValidRange()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act & Assert - default is 100.0
            Assert.That(viewModel.MaxRelativeCutoff, Is.EqualTo(100.0));

            // Act
            viewModel.MaxRelativeCutoff = 50.0;

            // Assert
            Assert.That(viewModel.MaxRelativeCutoff, Is.EqualTo(50.0));

            // Act
            viewModel.MaxRelativeCutoff = 0.0;

            // Assert
            Assert.That(viewModel.MaxRelativeCutoff, Is.EqualTo(0.0));
        }

        [Test]
        public void TestMaxRelativeCutoff_BelowZero_ClampedToZero()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act
            viewModel.MaxRelativeCutoff = -10.0;

            // Assert - should be clamped to 0
            Assert.That(viewModel.MaxRelativeCutoff, Is.EqualTo(0.0));
        }

        [Test]
        public void TestMaxRelativeCutoff_Above100_ClampedTo100()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act
            viewModel.MaxRelativeCutoff = 150.0;

            // Assert - should be clamped to 100
            Assert.That(viewModel.MaxRelativeCutoff, Is.EqualTo(100.0));
        }

        [Test]
        public void TestMaxRelativeCutoff_ConstraintVsMinRelativeCutoff()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            viewModel.MinRelativeCutoff = 70.0;

            // Act - set MaxRelativeCutoff below MinRelativeCutoff
            viewModel.MaxRelativeCutoff = 50.0;

            // Assert - should be clamped to MinRelativeCutoff
            Assert.That(viewModel.MaxRelativeCutoff, Is.EqualTo(70.0));
        }

        [Test]
        public void TestNormalizeHistogramToFile_GetterSetter()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act & Assert - default is false
            Assert.That(viewModel.NormalizeHistogramToFile, Is.False);

            // Act
            viewModel.NormalizeHistogramToFile = true;

            // Assert
            Assert.That(viewModel.NormalizeHistogramToFile, Is.True);

            // Act
            viewModel.NormalizeHistogramToFile = false;

            // Assert
            Assert.That(viewModel.NormalizeHistogramToFile, Is.False);
        }

        [Test]
        public void TestDisplayFilteredOnly_GetterSetter()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;

            // Act & Assert - default is true
            Assert.That(viewModel.DisplayFilteredOnly, Is.True);

            // Act
            viewModel.DisplayFilteredOnly = false;

            // Assert
            Assert.That(viewModel.DisplayFilteredOnly, Is.False);

            // Act
            viewModel.DisplayFilteredOnly = true;

            // Assert
            Assert.That(viewModel.DisplayFilteredOnly, Is.True);
        }

        [Test]
        public void TestGetParameters_ReturnsClone()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            viewModel.UseLogScaleYAxis = true;
            viewModel.GroupingProperty = "Precursor Charge";
            viewModel.AllowAmbiguousGroups = true;
            viewModel.MinRelativeCutoff = 5.0;
            viewModel.MaxRelativeCutoff = 95.0;
            viewModel.NormalizeHistogramToFile = true;
            viewModel.DisplayFilteredOnly = false;

            // Act
            var parameters = viewModel.GetParameters();

            // Assert - values match
            Assert.That(parameters.UseLogScaleYAxis, Is.EqualTo(viewModel.UseLogScaleYAxis));
            Assert.That(parameters.GroupingProperty, Is.EqualTo(viewModel.GroupingProperty));
            Assert.That(parameters.AllowAmbiguousGroups, Is.EqualTo(viewModel.AllowAmbiguousGroups));
            Assert.That(parameters.MinRelativeCutoff, Is.EqualTo(viewModel.MinRelativeCutoff));
            Assert.That(parameters.MaxRelativeCutoff, Is.EqualTo(viewModel.MaxRelativeCutoff));
            Assert.That(parameters.NormalizeHistogramToFile, Is.EqualTo(viewModel.NormalizeHistogramToFile));
            Assert.That(parameters.DisplayFilteredOnly, Is.EqualTo(viewModel.DisplayFilteredOnly));

            // Assert - returned object is a clone (modifying it doesn't affect the ViewModel)
            parameters.UseLogScaleYAxis = !viewModel.UseLogScaleYAxis;
            Assert.That(viewModel.UseLogScaleYAxis, Is.Not.EqualTo(parameters.UseLogScaleYAxis));
        }

        [Test]
        public void TestLoadFromSnapshot_UpdatesAllProperties()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            var snapshot = new PlotModelStatParameters
            {
                UseLogScaleYAxis = true,
                GroupingProperty = "File Name",
                AllowAmbiguousGroups = true,
                MinRelativeCutoff = 10.0,
                MaxRelativeCutoff = 90.0,
                NormalizeHistogramToFile = true,
                DisplayFilteredOnly = false
            };

            // Act
            viewModel.LoadFromSnapshot(snapshot);

            // Assert
            Assert.That(viewModel.UseLogScaleYAxis, Is.EqualTo(snapshot.UseLogScaleYAxis));
            Assert.That(viewModel.GroupingProperty, Is.EqualTo(snapshot.GroupingProperty));
            Assert.That(viewModel.AllowAmbiguousGroups, Is.EqualTo(snapshot.AllowAmbiguousGroups));
            Assert.That(viewModel.MinRelativeCutoff, Is.EqualTo(snapshot.MinRelativeCutoff));
            Assert.That(viewModel.MaxRelativeCutoff, Is.EqualTo(snapshot.MaxRelativeCutoff));
            Assert.That(viewModel.NormalizeHistogramToFile, Is.EqualTo(snapshot.NormalizeHistogramToFile));
            Assert.That(viewModel.DisplayFilteredOnly, Is.EqualTo(snapshot.DisplayFilteredOnly));
        }

        [Test]
        public void TestLoadFromSnapshot_OriginalSnapshotUnaffected()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            var snapshot = new PlotModelStatParameters
            {
                UseLogScaleYAxis = true,
                GroupingProperty = "File Name",
                AllowAmbiguousGroups = true,
                MinRelativeCutoff = 10.0,
                MaxRelativeCutoff = 90.0,
                NormalizeHistogramToFile = true,
                DisplayFilteredOnly = false
            };

            // Act
            viewModel.LoadFromSnapshot(snapshot);

            // Modify ViewModel
            viewModel.UseLogScaleYAxis = false;
            viewModel.GroupingProperty = "None";
            viewModel.AllowAmbiguousGroups = false;
            viewModel.MinRelativeCutoff = 0.0;
            viewModel.MaxRelativeCutoff = 100.0;
            viewModel.NormalizeHistogramToFile = false;
            viewModel.DisplayFilteredOnly = true;

            // Assert - snapshot is not affected
            Assert.That(snapshot.UseLogScaleYAxis, Is.True);
            Assert.That(snapshot.GroupingProperty, Is.EqualTo("File Name"));
            Assert.That(snapshot.AllowAmbiguousGroups, Is.True);
            Assert.That(snapshot.MinRelativeCutoff, Is.EqualTo(10.0));
            Assert.That(snapshot.MaxRelativeCutoff, Is.EqualTo(90.0));
            Assert.That(snapshot.NormalizeHistogramToFile, Is.True);
            Assert.That(snapshot.DisplayFilteredOnly, Is.False);
        }

        [Test]
        public void TestPropertyChanged_UseLogScaleYAxis()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            bool eventRaised = false;
            string propertyName = null;

            viewModel.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == nameof(PlotModelStatParametersViewModel.UseLogScaleYAxis))
                {
                    eventRaised = true;
                    propertyName = e.PropertyName;
                }
            };

            // Act
            viewModel.UseLogScaleYAxis = true;

            // Assert
            Assert.That(eventRaised, Is.True);
            Assert.That(propertyName, Is.EqualTo(nameof(PlotModelStatParametersViewModel.UseLogScaleYAxis)));
        }

        [Test]
        public void TestPropertyChanged_GroupingProperty()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            bool eventRaised = false;
            string propertyName = null;

            viewModel.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == nameof(PlotModelStatParametersViewModel.GroupingProperty))
                {
                    eventRaised = true;
                    propertyName = e.PropertyName;
                }
            };

            // Act
            viewModel.GroupingProperty = "Precursor Charge";

            // Assert
            Assert.That(eventRaised, Is.True);
            Assert.That(propertyName, Is.EqualTo(nameof(PlotModelStatParametersViewModel.GroupingProperty)));
        }

        [Test]
        public void TestPropertyChanged_AllowAmbiguousGroups()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            bool eventRaised = false;
            string propertyName = null;

            viewModel.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == nameof(PlotModelStatParametersViewModel.AllowAmbiguousGroups))
                {
                    eventRaised = true;
                    propertyName = e.PropertyName;
                }
            };

            // Act
            viewModel.AllowAmbiguousGroups = true;

            // Assert
            Assert.That(eventRaised, Is.True);
            Assert.That(propertyName, Is.EqualTo(nameof(PlotModelStatParametersViewModel.AllowAmbiguousGroups)));
        }

        [Test]
        public void TestPropertyChanged_MinRelativeCutoff()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            bool eventRaised = false;
            string propertyName = null;

            viewModel.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == nameof(PlotModelStatParametersViewModel.MinRelativeCutoff))
                {
                    eventRaised = true;
                    propertyName = e.PropertyName;
                }
            };

            // Act
            viewModel.MinRelativeCutoff = 50.0;

            // Assert
            Assert.That(eventRaised, Is.True);
            Assert.That(propertyName, Is.EqualTo(nameof(PlotModelStatParametersViewModel.MinRelativeCutoff)));
        }

        [Test]
        public void TestPropertyChanged_MaxRelativeCutoff()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            bool eventRaised = false;
            string propertyName = null;

            viewModel.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == nameof(PlotModelStatParametersViewModel.MaxRelativeCutoff))
                {
                    eventRaised = true;
                    propertyName = e.PropertyName;
                }
            };

            // Act
            viewModel.MaxRelativeCutoff = 50.0;

            // Assert
            Assert.That(eventRaised, Is.True);
            Assert.That(propertyName, Is.EqualTo(nameof(PlotModelStatParametersViewModel.MaxRelativeCutoff)));
        }

        [Test]
        public void TestPropertyChanged_NormalizeHistogramToFile()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            bool eventRaised = false;
            string propertyName = null;

            viewModel.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == nameof(PlotModelStatParametersViewModel.NormalizeHistogramToFile))
                {
                    eventRaised = true;
                    propertyName = e.PropertyName;
                }
            };

            // Act
            viewModel.NormalizeHistogramToFile = true;

            // Assert
            Assert.That(eventRaised, Is.True);
            Assert.That(propertyName, Is.EqualTo(nameof(PlotModelStatParametersViewModel.NormalizeHistogramToFile)));
        }

        [Test]
        public void TestPropertyChanged_DisplayFilteredOnly()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            bool eventRaised = false;
            string propertyName = null;

            viewModel.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == nameof(PlotModelStatParametersViewModel.DisplayFilteredOnly))
                {
                    eventRaised = true;
                    propertyName = e.PropertyName;
                }
            };

            // Act
            viewModel.DisplayFilteredOnly = false;

            // Assert
            Assert.That(eventRaised, Is.True);
            Assert.That(propertyName, Is.EqualTo(nameof(PlotModelStatParametersViewModel.DisplayFilteredOnly)));
        }

        [Test]
        public void TestPropertyChanged_LoadFromSnapshot_RaisesAllProperties()
        {
            // Arrange
            var viewModel = PlotModelStatParametersViewModel.Instance;
            var snapshot = new PlotModelStatParameters
            {
                UseLogScaleYAxis = true,
                GroupingProperty = "File Name",
                AllowAmbiguousGroups = true,
                MinRelativeCutoff = 10.0,
                MaxRelativeCutoff = 90.0,
                NormalizeHistogramToFile = true,
                DisplayFilteredOnly = false
            };

            var propertiesChanged = new System.Collections.Generic.HashSet<string>();

            viewModel.PropertyChanged += (s, e) =>
            {
                propertiesChanged.Add(e.PropertyName);
            };

            // Act
            viewModel.LoadFromSnapshot(snapshot);

            // Assert
            Assert.That(propertiesChanged, Contains.Item(nameof(PlotModelStatParametersViewModel.UseLogScaleYAxis)));
            Assert.That(propertiesChanged, Contains.Item(nameof(PlotModelStatParametersViewModel.GroupingProperty)));
            Assert.That(propertiesChanged, Contains.Item(nameof(PlotModelStatParametersViewModel.AllowAmbiguousGroups)));
            Assert.That(propertiesChanged, Contains.Item(nameof(PlotModelStatParametersViewModel.MinRelativeCutoff)));
            Assert.That(propertiesChanged, Contains.Item(nameof(PlotModelStatParametersViewModel.MaxRelativeCutoff)));
            Assert.That(propertiesChanged, Contains.Item(nameof(PlotModelStatParametersViewModel.NormalizeHistogramToFile)));
            Assert.That(propertiesChanged, Contains.Item(nameof(PlotModelStatParametersViewModel.DisplayFilteredOnly)));
        }

        #endregion
    }
}
