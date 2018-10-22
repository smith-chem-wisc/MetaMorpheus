using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using TaskLayer;

namespace Test
{
    class CalibrationTests
    {
        [Test]
        public static void ExperimentalDesignCalibrationTest()
        {
            CalibrationTask calibrationTask = new CalibrationTask();
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCalibration");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ExperimentalDesign.tsv");
            Directory.CreateDirectory(outputFolder);
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine("FileName\tCondition\tBiorep\tFraction\tTechrep");
                output.WriteLine("SmallCalibratible_Yeast" + "\t" + "condition" + "\t" + "1" + "\t" + "1" + "\t" + "1");
            }

            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
            var expDesignPath = Path.Combine(outputFolder, @"ExperimentalDesign.tsv");
            var expDesign = File.ReadAllLines(expDesignPath);
            Assert.That(expDesign[1].Contains("SmallCalibratible_Yeast-calib"));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"SmallCalibratible_Yeast-calib.mzML")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"SmallCalibratible_Yeast-calib.toml")));
            var lines = File.ReadAllLines(Path.Combine(outputFolder, @"SmallCalibratible_Yeast-calib.toml"));
            var tolerance = Regex.Match(lines[0], @"\d+\.\d*").Value;
            var tolerance1 = Regex.Match(lines[1], @"\d+\.\d*").Value;
            Assert.That(double.TryParse(tolerance, out double tol) == true);
            Assert.That(double.TryParse(tolerance1, out double tol1) == true);
            Assert.That(lines[0].Contains("PrecursorMassTolerance"));
            Assert.That(lines[1].Contains("ProductMassTolerance"));
            File.Delete(filePath);
            Directory.Delete(outputFolder, true);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }
    }
}
