using EngineLayer;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class SelectedPrecursorsTests
    {
        private string tempFilePath;

        [SetUp]
        public void Setup()
        {
            tempFilePath = Path.GetTempFileName();
        }

        [TearDown]
        public void TearDown()
        {
            if (File.Exists(tempFilePath))
            {
                File.Delete(tempFilePath);
            }
        }

        [Test]
        public void PrecursorInfo_Constructor_ValidValues_PropertiesSet()
        {
            // Arrange & Act
            var precursorInfo = new SelectedPrecursors.PrecursorInfo(500.5, 2, 10.5, 11.5);

            // Assert
            Assert.That(precursorInfo.Mz, Is.EqualTo(500.5));
            Assert.That(precursorInfo.Charge, Is.EqualTo(2));
            Assert.That(precursorInfo.RtStartInMinutes, Is.EqualTo(10.5));
            Assert.That(precursorInfo.RtEndInMinutes, Is.EqualTo(11.5));
        }

        [Test]
        public void PrecursorInfo_ToString_ReturnsTabSeparatedValues()
        {
            // Arrange
            var precursorInfo = new SelectedPrecursors.PrecursorInfo(500.5, 2, 10.5, 11.5);
            
            // Act
            string result = precursorInfo.ToString();
            
            // Assert
            Assert.That(result, Is.EqualTo("500.5\t2\t10.5\t11.5"));
        }

        [Test]
        public void PrecursorInfo_ToPrecursor_ConvertsToPrecursor()
        {
            // Arrange
            var precursorInfo = new SelectedPrecursors.PrecursorInfo(500.5, 2, 10.5, 11.5);
            
            // Act
            var precursor = precursorInfo.ToPrecursor();
            
            // Assert
            Assert.That(precursor.MonoisotopicPeakMz, Is.EqualTo(500.5));
            Assert.That(precursor.Charge, Is.EqualTo(2));
        }

        [Test]
        public void ReadSelectedPrecursors_FileNotFound_ReturnsEmptyListAndError()
        {
            // Arrange
            string nonExistentFilePath = "nonexistent.tsv";
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(nonExistentFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Is.EqualTo("Selected precursors file not found!"));
        }

        [Test]
        public void ReadSelectedPrecursors_MalformedLine_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t2\t10.5" // Missing RT End column
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("Expected 4 cells, but found 3"));
        }

        [Test]
        public void ReadSelectedPrecursors_InvalidMz_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "invalid\t2\t10.5\t11.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("m/z value on line 2 is not a valid number"));
        }

        [Test]
        public void ReadSelectedPrecursors_InvalidCharge_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\tinvalid\t10.5\t11.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("charge on line 2 is not a valid integer"));
        }

        [Test]
        public void ReadSelectedPrecursors_InvalidRtStart_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t2\tinvalid\t11.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("RT Start value on line 2 is not a valid number"));
        }

        [Test]
        public void ReadSelectedPrecursors_InvalidRtEnd_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t2\t10.5\tinvalid"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("RT End value on line 2 is not a valid number"));
        }

        [Test]
        public void ReadSelectedPrecursors_ZeroMz_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "0\t2\t10.5\t11.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("m/z value on line 2 must be greater than zero"));
        }

        [Test]
        public void ReadSelectedPrecursors_ZeroCharge_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t0\t10.5\t11.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("charge on line 2 can not be zero"));
        }

        [Test]
        public void ReadSelectedPrecursors_NegativeRtStart_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t2\t-1.0\t11.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("RT Start value on line 2 cannot be negative"));
        }

        [Test]
        public void ReadSelectedPrecursors_RtEndLessThanRtStart_ReturnsEmptyListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t2\t11.5\t10.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Is.Empty);
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("RT End value must be greater than or equal to RT Start"));
        }

        [Test]
        public void ReadSelectedPrecursors_DuplicatePrecursors_ReturnsListAndError()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t2\t10.5\t11.5\n" +
                "500.5\t2\t12.5\t13.5" // Duplicate m/z and charge
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(result, Has.Count.EqualTo(2));
            Assert.That(errors, Has.Count.EqualTo(1));
            Assert.That(errors[0], Contains.Substring("Duplicate precursors found"));
        }

        [Test]
        public void ReadSelectedPrecursors_ValidFile_ReturnsCorrectPrecursors()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t2\t10.5\t11.5\n" +
                "600.6\t3\t12.5\t13.5\n" +
                "700.7\t4\t14.5\t15.5"
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(errors, Is.Empty);
            Assert.That(result, Has.Count.EqualTo(3));
            
            Assert.That(result[0].Mz, Is.EqualTo(500.5));
            Assert.That(result[0].Charge, Is.EqualTo(2));
            Assert.That(result[0].RtStartInMinutes, Is.EqualTo(10.5));
            Assert.That(result[0].RtEndInMinutes, Is.EqualTo(11.5));
            
            Assert.That(result[1].Mz, Is.EqualTo(600.6));
            Assert.That(result[1].Charge, Is.EqualTo(3));
            Assert.That(result[1].RtStartInMinutes, Is.EqualTo(12.5));
            Assert.That(result[1].RtEndInMinutes, Is.EqualTo(13.5));
            
            Assert.That(result[2].Mz, Is.EqualTo(700.7));
            Assert.That(result[2].Charge, Is.EqualTo(4));
            Assert.That(result[2].RtStartInMinutes, Is.EqualTo(14.5));
            Assert.That(result[2].RtEndInMinutes, Is.EqualTo(15.5));
        }

        [Test]
        public void WriteSelectedPrecursorsToFile_EmptyList_WritesHeaderOnly()
        {
            // Arrange
            var precursors = new List<SelectedPrecursors.PrecursorInfo>();
            string outputDir = Path.GetDirectoryName(tempFilePath);
            string outputFileName = Path.GetFileName(tempFilePath);
            
            // Act
            string filePath = SelectedPrecursors.WriteSelectedPrecursorsToFile(precursors, outputDir, outputFileName);
            string fileContent = File.ReadAllText(filePath);
            
            // Assert
            Assert.That(filePath, Is.EqualTo(tempFilePath));
            Assert.That(fileContent.Trim(), Is.EqualTo("m/z\tCharge\tRT Start (min)\tRT End (min)"));
        }

        [Test]
        public void WriteSelectedPrecursorsToFile_ValidPrecursors_WritesCorrectly()
        {
            // Arrange
            var precursors = new List<SelectedPrecursors.PrecursorInfo>
            {
                new SelectedPrecursors.PrecursorInfo(500.5, 2, 10.5, 11.5),
                new SelectedPrecursors.PrecursorInfo(600.6, 3, 12.5, 13.5),
                new SelectedPrecursors.PrecursorInfo(700.7, 4, 14.5, 15.5)
            };
            string outputDir = Path.GetDirectoryName(tempFilePath);
            string outputFileName = Path.GetFileName(tempFilePath);
            
            // Act
            string filePath = SelectedPrecursors.WriteSelectedPrecursorsToFile(precursors, outputDir, outputFileName);
            string[] lines = File.ReadAllLines(filePath);
            
            // Assert
            Assert.That(filePath, Is.EqualTo(tempFilePath));
            Assert.That(lines.Length, Is.EqualTo(4)); // Header + 3 precursors
            Assert.That(lines[0], Is.EqualTo("m/z\tCharge\tRT Start (min)\tRT End (min)"));
            Assert.That(lines[1], Is.EqualTo("500.5\t2\t10.5\t11.5"));
            Assert.That(lines[2], Is.EqualTo("600.6\t3\t12.5\t13.5"));
            Assert.That(lines[3], Is.EqualTo("700.7\t4\t14.5\t15.5"));
        }

        [Test]
        public void GetErrorsInSelectedPrecursors_NullList_ReturnsError()
        {
            // Act
            string error = SelectedPrecursors.GetErrorsInSelectedPrecursors(null);
            
            // Assert
            Assert.That(error, Is.EqualTo("No precursors defined!"));
        }

        [Test]
        public void GetErrorsInSelectedPrecursors_EmptyList_ReturnsError()
        {
            // Arrange
            var precursors = new List<SelectedPrecursors.PrecursorInfo>();
            
            // Act
            string error = SelectedPrecursors.GetErrorsInSelectedPrecursors(precursors);
            
            // Assert
            Assert.That(error, Is.EqualTo("No precursors defined!"));
        }

        [Test]
        public void GetErrorsInSelectedPrecursors_OverlappingRtWindows_ReturnsError()
        {
            // Arrange
            var precursors = new List<SelectedPrecursors.PrecursorInfo>
            {
                new SelectedPrecursors.PrecursorInfo(500.5, 2, 10.5, 12.5),
                new SelectedPrecursors.PrecursorInfo(500.5, 2, 12.0, 13.5) // Overlapping with first precursor
            };
            
            // Act
            string error = SelectedPrecursors.GetErrorsInSelectedPrecursors(precursors);
            
            // Assert
            Assert.That(error, Contains.Substring("Overlapping retention time windows found"));
            Assert.That(error, Contains.Substring("m/z 500.5, charge 2"));
        }

        [Test]
        public void GetErrorsInSelectedPrecursors_InvalidRtRange_ReturnsError()
        {
            // Arrange
            var precursors = new List<SelectedPrecursors.PrecursorInfo>
            {
                new SelectedPrecursors.PrecursorInfo(500.5, 2, 12.5, 12.5) // Equal RT start and end
            };
            
            // Act
            string error = SelectedPrecursors.GetErrorsInSelectedPrecursors(precursors);
            
            // Assert
            Assert.That(error, Contains.Substring("Invalid retention time range"));
            Assert.That(error, Contains.Substring("RT Start (12.5) must be less than RT End (12.5)"));
        }

        [Test]
        public void GetErrorsInSelectedPrecursors_ValidPrecursors_ReturnsNull()
        {
            // Arrange
            var precursors = new List<SelectedPrecursors.PrecursorInfo>
            {
                new SelectedPrecursors.PrecursorInfo(500.5, 2, 10.5, 11.5),
                new SelectedPrecursors.PrecursorInfo(600.6, 3, 12.5, 13.5),
                // Different m/z, same charge
                new SelectedPrecursors.PrecursorInfo(700.7, 2, 14.5, 15.5),
                // Same m/z, different charge
                new SelectedPrecursors.PrecursorInfo(500.5, 3, 16.5, 17.5),
                // Same m/z and charge, non-overlapping RT windows
                new SelectedPrecursors.PrecursorInfo(800.8, 4, 18.5, 19.5),
                new SelectedPrecursors.PrecursorInfo(800.8, 4, 20.5, 21.5)
            };
            
            // Act
            string error = SelectedPrecursors.GetErrorsInSelectedPrecursors(precursors);
            
            // Assert
            Assert.That(error, Is.Null);
        }

        [Test]
        public void ReadSelectedPrecursors_CultureWithCommaDecimalSeparator_ParsesCorrectly()
        {
            // Arrange
            var originalCulture = CultureInfo.CurrentCulture;
            try
            {
                // Set culture with comma as decimal separator
                CultureInfo.CurrentCulture = new CultureInfo("fr-FR");
                
                File.WriteAllText(tempFilePath, 
                    "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                    "500.5\t2\t10.5\t11.5" // Using period as decimal separator in file
                );
                
                // Act
                var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
                
                // Assert
                Assert.That(errors, Is.Empty);
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0].Mz, Is.EqualTo(500.5));
                Assert.That(result[0].RtStartInMinutes, Is.EqualTo(10.5));
                Assert.That(result[0].RtEndInMinutes, Is.EqualTo(11.5));
            }
            finally
            {
                // Restore original culture
                CultureInfo.CurrentCulture = originalCulture;
            }
        }

        [Test]
        public void WriteSelectedPrecursorsToFile_CultureWithCommaDecimalSeparator_WritesWithPeriod()
        {
            // Arrange
            var originalCulture = CultureInfo.CurrentCulture;
            try
            {
                // Set culture with comma as decimal separator
                CultureInfo.CurrentCulture = new CultureInfo("fr-FR");
                
                var precursors = new List<SelectedPrecursors.PrecursorInfo>
                {
                    new SelectedPrecursors.PrecursorInfo(500.5, 2, 10.5, 11.5)
                };
                string outputDir = Path.GetDirectoryName(tempFilePath);
                string outputFileName = Path.GetFileName(tempFilePath);
                
                // Act
                string filePath = SelectedPrecursors.WriteSelectedPrecursorsToFile(precursors, outputDir, outputFileName);
                string[] lines = File.ReadAllLines(filePath);
                
                // Assert
                Assert.That(lines.Length, Is.EqualTo(2)); // Header + 1 precursor
                Assert.That(lines[1], Is.EqualTo("500.5\t2\t10.5\t11.5")); // Periods as decimal separators
            }
            finally
            {
                // Restore original culture
                CultureInfo.CurrentCulture = originalCulture;
            }
        }

        [Test]
        public void ReadSelectedPrecursors_NegativeCharge_IsAllowed()
        {
            // Arrange
            File.WriteAllText(tempFilePath, 
                "m/z\tCharge\tRT Start (min)\tRT End (min)\n" +
                "500.5\t-2\t10.5\t11.5" // Negative charge
            );
            
            // Act
            var result = SelectedPrecursors.ReadSelectedPrecursors(tempFilePath, out var errors);
            
            // Assert
            Assert.That(errors, Is.Empty);
            Assert.That(result, Has.Count.EqualTo(1));
            Assert.That(result[0].Mz, Is.EqualTo(500.5));
            Assert.That(result[0].Charge, Is.EqualTo(-2));
            Assert.That(result[0].RtStartInMinutes, Is.EqualTo(10.5));
            Assert.That(result[0].RtEndInMinutes, Is.EqualTo(11.5));
        }
    }
}