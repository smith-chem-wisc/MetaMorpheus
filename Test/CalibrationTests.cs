using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using TaskLayer;

namespace Test
{
    class CalibrationTests
    {
        public static void RunCalibrationTest()
        {
            CalibrationTask calibrationTask = new CalibrationTask();
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCalibration");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.raw");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
        }
    }
}
