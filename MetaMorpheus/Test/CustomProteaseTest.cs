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
using MzLibUtil;


namespace Test
{
    public static class CustomProteaseTest
    {
        [Test]
        public static void TestAddCustomProtease ()
        {
            string ProtDirectory = Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion");
            string customProteasePath = Path.Combine(ProtDirectory, @"CustomProtease.tsv");

            string[] pString = new string[] { "p\tA|\t\t\tfull\tcusPsiNum\tcusPsiName\t\t\t\t\t" };
            File.WriteAllLines(customProteasePath, pString);

            GlobalVariables.LoadCustomProtease();
            var dictionary = ProteaseDictionary.Dictionary;
            Assert.That(ProteaseDictionary.Dictionary.ContainsKey("p"));
            Assert.That(ProteaseDictionary.Dictionary["p"].PsiMsName == "cusPsiName");

            File.Delete(customProteasePath);
        }
    }
}
