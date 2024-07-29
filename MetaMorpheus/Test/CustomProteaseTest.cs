using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics.ProteolyticDigestion;
using System.IO;
using NUnit.Framework;
using EngineLayer;
using Omics.Digestion;


namespace Test
{
    public static class CustomProteaseTest
    {
        [Test]
        public static void TestAddCustomProtease ()
        {
            string folderPath = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            string path = ((!string.IsNullOrWhiteSpace(folderPath) && AppDomain.CurrentDomain.BaseDirectory.Contains(folderPath) && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins")) ? Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus") : AppDomain.CurrentDomain.BaseDirectory);
            string path2 = Path.Combine(path, "ProteolyticDigestion", "proteases.tsv");

            string ProtDirectory = Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion");
            string customProteasePath = Path.Combine(ProtDirectory, @"CustomProtease.tsv");
            if (!File.Exists(customProteasePath))
            {
                File.Create(customProteasePath);
            }

            File.Copy(path2, customProteasePath, true);

            Protease p = new Protease("p", CleavageSpecificity.Full, "cusPsiNum", "cusPsiName", DigestionMotif.ParseDigestionMotifsFromString("A|"));
            string pString = "\np\tA|\t\t\tfull\tcusPsiNum\tcusPsiName\t\t\t\t\t";

            File.AppendAllText(customProteasePath, pString);
            Dictionary<string, Protease> dictionary = ProteaseDictionary.LoadProteaseDictionary(customProteasePath);

            string[] lines = File.ReadAllLines(customProteasePath);
        }
    }
}
