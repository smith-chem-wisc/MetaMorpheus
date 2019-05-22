using Chemistry;
using Proteomics.AminoAcidPolymer;
using System.IO;
using NUnit.Framework;
using EngineLayer;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public static class CustomAminoAcidsTest
    {
        [Test]
        public static void TestCustomAminoAcidReading()
        {
            string aminoAcidPath = Path.Combine(GlobalVariables.DataDir, @"CustomAminoAcids", @"CustomAminoAcids.txt");

            //Manually add an entry to it
            List<string> lines = new List<string>(File.ReadAllLines(aminoAcidPath));
            lines.Add("fake\tf\t60\tC5");
            File.WriteAllLines(aminoAcidPath, lines);

            GlobalVariables.RefreshAminoAcidDictionary(); //read the file

            //test that we read the new amino acid
            Assert.IsTrue(Residue.TryGetResidue('f', out Residue r));
            Assert.IsTrue(r.MonoisotopicMass.Equals(60));

            //now crash it intentionally with an invalid character
            lines.Add("evenFaker\tX\t72\tC6");
            File.WriteAllLines(aminoAcidPath, lines);
            try
            {
                GlobalVariables.RefreshAminoAcidDictionary(); //read the file
                //we're trying to crash it, so if we didn't, we failed :/
                Assert.IsTrue(false);
            }
            catch(MetaMorpheusException)
            {
                //Yay we passed!
            }

            //now crash it intentionally with a bad chemical formula
            lines.RemoveAt(lines.Count - 1); //get rid of that last bad one
            lines.Add("theFakest\ta\t50\t");
            File.WriteAllLines(aminoAcidPath, lines);
            try
            {
                GlobalVariables.RefreshAminoAcidDictionary(); //read the file
                //we're trying to crash it, so if we didn't, we failed :/
                Assert.IsTrue(false);
            }
            catch (MetaMorpheusException)
            {
                //Yay we passed!
            }

            //Delete so it doesn't crash the next time
            File.Delete(aminoAcidPath);
        }
    }
}
