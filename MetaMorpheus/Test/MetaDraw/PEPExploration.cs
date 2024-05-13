using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using NUnit.Framework;
using TaskLayer;

namespace Test.MetaDraw
{
    internal class PEPExploration
    {
        //public static string DirectoryPath = @"B:\Users\Nic\Chimeras\PEPTesting";
        public static string DirectoryPath = @"B:\Users\Nic\Chimeras\TopDown_Analysis\Jurkat\SearchResults";
        public static string TwoRawFiles = @"B:\Users\Nic\Chimeras\PEPTesting\TwoRawFiles_SerializedPostSearhAnalysisTask.txt";
        public static string EntireTask = @"B:\Users\Nic\Chimeras\PEPTesting\FullJurkat_SerializedPostSearhAnalysisTask.txt";

        [Test]
        public static void TESTNAME()
        {
            // setup
            var task = ByteSerializer.ByteArrayFileToObject<PostSearchAnalysisTask>(EntireTask);
            GlobalVariables.AnalyteType = "Proteoform";
            string outputFolder = Path.Combine(DirectoryPath, "MetaMorpheus_NewPEP_NoNormNoMult");
            if (!Directory.Exists(outputFolder))
                Directory.CreateDirectory(outputFolder);
            task.Parameters.OutputFolder = outputFolder;


            task.Run();
        }

    }
}
