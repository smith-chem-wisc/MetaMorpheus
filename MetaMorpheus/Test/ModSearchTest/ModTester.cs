using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.IO;
using System.Runtime.CompilerServices;
using TaskLayer;
using UsefulProteomicsDatabases;
using Nett;
using EngineLayer.GlycoSearch;
using FlashLFQ;
using SpectralAveraging;
using NUnit.Framework.Legacy;
using Omics.Modifications;
using ThermoFisher.CommonCore.BackgroundSubtraction;
using Easy.Common.Extensions;
using iText.IO.Font.Otf;
using static Nett.TomlObjectFactory;
using Omics.SpectrumMatch;
using TopDownProteomics;
using MzLibUtil;
using EngineLayer.ModSearch;
using Org.BouncyCastle.Asn1.BC;
using ThermoFisher.CommonCore.Data;
using SearchType = EngineLayer.ModSearch.SearchType;

namespace Test.ModSearchTest
{
    public class ModTester
    {
        // ModBox tester
        [Test]
        public static void ModBoxTest()
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            Modification mod = new Modification(_originalId: "TestMod0", _modificationType: "regularMod", _monoisotopicMass: 100.0, _target: motif, _locationRestriction: "Anywhere");
            Modification mod2 = new Modification(_originalId: "TestMod1", _modificationType: "O-Glycosylation", _monoisotopicMass: 100.0, _target: motif, _locationRestriction: "Anywhere");
            Modification mod3 = new Modification(_originalId: "TestMod2", _modificationType: "N-Glycosylation", _monoisotopicMass: 100.0, _target: motif, _locationRestriction: "Anywhere");
            ModBox modBox = new ModBox(new List<Modification> { mod, mod2, mod3 }, [0, 1, 2]);

            Assert.NotNull(modBox.ChildBoxes);
            Assert.AreEqual(modBox.ChildBoxes.Count(), 8);  //null, Mod1, Mod2, Mod3, Mod1+Mod2, Mod1+Mod3, Mod2+Mod3, Mod1+Mod2+Mod3
            Assert.AreEqual(modBox.ModIds.Length, 3);
            Assert.AreEqual(modBox.NumberOfMods, 3);
            Assert.AreEqual(modBox.Mass, 300.0);
            Assert.AreEqual(modBox.TargetDecoy, false);
            Assert.AreEqual(modBox.IsChild, false);
            Assert.AreEqual(modBox.SugarKind.Length, 11);
            Assert.IsTrue(modBox.IsChild == false); // The parent box should not be labeled as child
            Assert.IsTrue(modBox.ChildBoxes.Where(p => p.IsChild).Count() == 8); // all of the childBox should be labeled as child

            // In this modBox, there are 6 child boxes
            ModBox modBox_duplicate = new ModBox(new List<Modification> { mod, mod2, mod2 }, [0, 1, 1]);
            Assert.AreEqual(modBox_duplicate.ChildBoxes.Count(), 6);//null, Mod1, Mod2, Mod1+Mod2,  Mod2+Mod2, Mod1+Mod2+Mod2
        }

        // ModDatabaseReader tester
        [Test]
        public static void ModDatabaseReaderTest()
        {
            //Description: This test checks the ModDatabaseReader functionality with a sample mod box database file.
            //In this database, there are 10 glycans, 2 regular modifications.
            //We will generate two reader, one is default reader, the other is O-glycan search reader.

            // The default reader only read regular mod, so the mod number is 2.
            ModDatabaseReader reader = new ModDatabaseReader(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Glycan_Mods", @"General_Mods", @"ModBox Database.txt"));
            
            //Check the mod number and type
            Assert.AreEqual(reader.GlobalGlycans, null);
            Assert.AreEqual(reader.GlobalRegularMods.Count, 2);
            Assert.AreEqual(ModDatabaseReader.ModDictionary.Count(), 166); // There are total 166 modifications in the MetaMorpheus

            //Check the default setting
            Assert.AreEqual(reader.MaxModNum, 3);
            Assert.AreEqual(reader.MaxGlycanNum,3);
            Assert.AreEqual(reader.SearchType, SearchType.RegularModSearch);
            Assert.IsTrue(reader.ToGenerateIons == false);

            ModDatabaseReader reader_2 = new ModDatabaseReader(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Glycan_Mods", @"General_Mods", @"ModBox Database.txt"), 2,3,false,SearchType.O_GlycanSearch);

            //Check the mod number and type
            Assert.AreEqual(reader_2.GlobalGlycans.Count, 10);
            Assert.AreEqual(reader_2.GlobalRegularMods.Count, 2);

            //check the setting
            Assert.AreEqual(reader_2.MaxModNum, 3);
            Assert.AreEqual(reader_2.MaxGlycanNum, 2);
            Assert.AreEqual(reader_2.SearchType, SearchType.O_GlycanSearch);
            Assert.IsTrue(reader_2.ToGenerateIons == false);

            Assert.AreEqual(reader_2.ModBoxes.Where(p => p.IsChild ).Count(), 0); // The box in the reader will be parent box, and they should not be labeled as child

            // Exception test: maxGlycanNum larger than maxModNum
            try
            {
                ModDatabaseReader reader_3 = new ModDatabaseReader(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Glycan_Mods", @"General_Mods", @"ModBox Database.txt"), maxModNum: 2, maxGlycanNum: 3, searchType: SearchType.O_GlycanSearch);
            }
            catch (ArgumentException ex)
            {
                Assert.AreEqual("The maximum number of glycans cannot be greater than the maximum number of total modifications.", ex.Message);
            }
        }

        [Test]
        public static void SugarKindTester()
        {
            // Description: This test checks the SugarKind functionality in the ModBox class.
            // The SugarKind is used to store the sugar composition of the glycan modifications in this ModBox, used for O-glycan search.

            ModDatabaseReader reader = new ModDatabaseReader(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Glycan_Mods", @"General_Mods", @"ModBox Database.txt"), searchType: SearchType.O_GlycanSearch);
            ModBox modBox = reader.ModBoxes[94]; // N1 + H1N1 + Oxidatino on M
            Assert.AreEqual(modBox.SugarKind.Length, 11);
            CollectionAssert.AreEqual(modBox.SugarKind, new int[] { 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

            ModBox modBox1 = new ModBox(new List<Modification> { reader.GlobalMods[0], reader.GlobalMods[0], reader.GlobalMods[0] }, [0, 0, 0]); // box1 = N1 + N1 + N1
            ModBox modBox2 = new ModBox(new List<Modification> { reader.GlobalMods[0], reader.GlobalMods[1], reader.GlobalMods[2] }, [0, 1, 2]); // box2 = N1+ H1N1 + N2

            Assert.AreEqual(modBox1.SugarKind.Length, 11);
            CollectionAssert.AreEqual(modBox1.SugarKind, new int[] { 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0 }); // Total = 3*N
            Assert.AreEqual(modBox2.SugarKind.Length, 11);
            CollectionAssert.AreEqual(modBox2.SugarKind, new int[] { 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0 }); // Total = 1*H + 4*N 
        }


        [Test]
        public static void ModBoxCombinationTest()
        {
            // Description: This test checks the ModBox combination functionality in the ModDatabaseReader class.
            // The ModBox combination is used to generate all possible combinations of modifications in the ModBox class.
            // In this database, there are 10 glycans, 2 regular modifications.

            ModDatabaseReader reader = new ModDatabaseReader(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Glycan_Mods", @"General_Mods", @"ModBox Database.txt"), searchType: SearchType.O_GlycanSearch);
            reader.BuildModBoxes();
            Assert.IsTrue(!reader.ModBoxes.Where(p => p.IsChild == true).Any()); // The parent box should not be labeled as child
            Assert.AreEqual(reader.ModBoxes.Where(p => p.ModIds.Count() == 1).Count(), 12); // 12
            Assert.AreEqual(reader.ModBoxes.Where(p => p.ModIds.Count() == 2).Count(), 78); // C12_2 + 12 = 78
            Assert.AreEqual(reader.ModBoxes.Where(p => p.ModIds.Count() == 3).Count(), 364); // 12
        }


    }   
}
