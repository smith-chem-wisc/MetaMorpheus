using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using TaskLayer;
using MzLibUtil;
using Nett;
using NUnit.Framework.Legacy;
using Omics.Modifications;

namespace Test
{
    [TestFixture]
    public class TestNGlyco
    {
        [Test]
        public static void TestNGlycoPsmsHeader()
        {
            List<string> headerTerms = new()
            {
                "File Name",
                "Scan Number",
                "Scan Retention Time",
                "Precursor Scan Number",
                "Precursor MZ",
                "Precursor Charge",
                "Precursor Mass",
                "Protein Accession",
                "Organism Name",
                "Protein Name",
                "Missed Cleavages",
                "Start and End Residues In Protein",
                "Base Sequence",
                "Flanking Residues",
                "Full Sequence",
                "Number of Mods",
                "Peptide Monoisotopic Mass",
                "Score",
                "Rank",
                "Matched Ion Series",
                "Matched Ion Mass-To-Charge Ratios",
                "Matched Ion Mass Diff (Da)",
                "Matched Ion Mass Diff (Ppm)",
                "Matched Ion Intensities",
                "Matched Ion Counts",
                "Decoy/Contaminant/Target",
                "QValue",
                "PEP",
                "PEP_QValue",
                "Localization Score",
                "Yion Score",
                "DiagonosticIon Score",
                "Plausible Number Of Glycans",//Not used for N-Glyco
                "Total Glycosylation sites",//Not used for N-Glyco
                "GlycanMass",
                "Plausible GlycanComposition",
                "N-Glycan motif Check",//Not used for N-Glyco
                "R138/144",
                "Plausible GlycanStructure",
                "GlycanLocalizationLevel",
                "Localized Glycans with Peptide Site Specific Probability",
                "Localized Glycans with Protein Site Specific Probability",
                "All potential glycan localizations",//Not used for N-Glyco
                "All SiteSpecific Localization Probability",//Not used for N-Glyco
            };

            headerTerms = headerTerms.Select(p => p.ToLower()).ToList();

            string nglycoHeaderString = GlycoSpectralMatch.GetTabSepHeaderSingle() + GlycoSpectralMatch.GetTabSeperatedHeaderGlyco();
            List<string> nGlycoHeaderTerms = nglycoHeaderString.Split('\t').Select(p => p.ToLower()).ToList();

            CollectionAssert.AreEquivalent(headerTerms, nGlycoHeaderTerms);
        }

        [Test]
        public static void GlyTest_GetKindString()
        {
            byte[] kind = new byte[] {3, 4, 0, 0, 1, 0, 0, 0, 0, 0, 0 };
            string kindString = Glycan.GetKindString(kind);
            Assert.That(kindString, Is.EqualTo("H3N4F1"));
        }

        [Test]
        public static void GlyTest_ModificationSites()
        {
            PeptideWithSetModifications pep = new PeptideWithSetModifications("ELNPTPNVEVNVECR", null); 
            string[] motifs = new string[] { "Nxs", "Nxt"};
            var sites = GlycoSpectralMatch.GetPossibleModSites(pep, motifs).Select(p => p.Key).ToList();
            Assert.That(sites.Count() == 1 && sites[0] == 4);

            ModificationMotif.TryGetMotif("C", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Test of N", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.");
            var testN = new PeptideWithSetModifications("C[Common Fixed:Carbamidomethyl on C]N[Common Fixed:Test of N]SSDQPKL[Common Fixed:Carbamidomethyl on C]NLSGIETP", new Dictionary<string, Modification> { { "Carbamidomethyl on C", mod1 }, { "Test of N", mod2 } });
            var testSites = GlycoSpectralMatch.GetPossibleModSites(testN, motifs).Select(p => p.Key).ToList();
            Assert.That(testSites.Count() == 1 && testSites[0] == 11);


            var testC = new PeptideWithSetModifications("TELAAYLSC[Common Fixed:Carbamidomethyl on C]NATK", new Dictionary<string, Modification> { { "Carbamidomethyl on C", mod1 }});
            var testCSites = GlycoSpectralMatch.GetPossibleModSites(testC, motifs).Select(p => p.Key).ToList();
            Assert.That(testCSites.Count() == 1 && testCSites[0] == 11);
        }

        [Test]
        public static void GlyTest_GlyGetTheoreticalFragments()
        {
            Protein pep = new Protein("TKPREEQYNSTYR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 7);
            var aPeptideWithSetModifications = pep.Digest(digestionParams, new List<Modification>(), new List<Modification>());

            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = GlycoSpectralMatch.GetPossibleModSites(aPeptideWithSetModifications.Last(), motifs).Select(p => p.Key).ToList();
            Glycan glycan = Glycan.Struct2Glycan("(N(F)(N(H(H(N))(H(N)))))", 0).FirstOrDefault();

            CommonParameters commonParameters = new CommonParameters(deconvolutionMassTolerance: new PpmTolerance(20), trimMsMsPeaks: false);
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/Glyco_3383.mgf"); //"25170.mgf"
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();

            var glycopep = GlycoPeptides.GenerateGlycopeptide(sites[0], aPeptideWithSetModifications.Last(), glycan);
            List<Product> fragmentIons = new List<Product>();
            glycopep.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragmentIons);

            var glycanYIons = GlycoPeptides.GetGlycanYIons(listOfSortedms2Scans[0].PrecursorMass, glycan);
            var matchedGlycanYIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], glycanYIons, commonParameters);
            Assert.That(matchedGlycanYIons.Count, Is.EqualTo(14));

            //TO DO: The neutroloss is not annotated well.
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);

            var coreIons = GlycoPeptides.ScanGetTrimannosylCore(matchedFragmentIons, glycan);
            Assert.That(coreIons.Count, Is.EqualTo(6));
            var filter = GlycoPeptides.ScanTrimannosylCoreFilter(matchedFragmentIons, glycan);
            Assert.That(filter, Is.EqualTo(true));
            var NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths[0], true, false);
            var bestGlycans = GlycoPeptides.MatchBestGlycan(listOfSortedms2Scans[0], NGlycans.ToArray(), commonParameters).Where(p => p != null && p.Item2 >= 2).OrderByDescending(p => p.Item2).Take(100).OrderBy(p => p.Item3).ToArray(); ;

        }

        [Test]
        public static void GlyTest_RunTask()
        {
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);

            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\Q9C0Y4.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\yeast_glycan_25170.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask>
                {
                    db
                }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        [Test]
        public static void GlyTest_AIETD()
        {
            Protein pep = new Protein("TNSSFIQGFVDHVKEDCDR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 19);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var fixedModifications = new List<Modification>() { mod2 };
            var aPeptideWithSetModifications = pep.Digest(digestionParams, fixedModifications, new List<Modification>());

            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = GlycoSpectralMatch.GetPossibleModSites(aPeptideWithSetModifications.Last(), motifs).Select(p => p.Key).ToList();
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(H(H(H)))(H(H(H(H(H))))))))", 0).FirstOrDefault();

            Tolerance tolerance = new PpmTolerance(20);
            CommonParameters commonParameters = new CommonParameters(doPrecursorDeconvolution:false, trimMsMsPeaks:false, dissociationType:DissociationType.EThcD, productMassTolerance: tolerance);
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/11901_AIETD.mgf"); 
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();

            var XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(commonParameters.PrecursorMassTolerance.Value);
            var precusorMatched = XLPrecusorSearchMode.Accepts(aPeptideWithSetModifications.Last().MonoisotopicMass + (double)glycan.Mass/1E5, listOfSortedms2Scans[0].PrecursorMass);
            Assert.That(precusorMatched, Is.EqualTo(0));

            var glycopep = GlycoPeptides.GenerateGlycopeptide(sites[0], aPeptideWithSetModifications.Last(), glycan);
            List<Product> fragmentIons = new List<Product>();
            glycopep.Fragment(DissociationType.EThcD, FragmentationTerminus.Both, fragmentIons);
               
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);

            using (StreamWriter output = new StreamWriter(Path.Combine(TestContext.CurrentContext.TestDirectory, "11091_NGlyco_AIETD.tsv")))
            {
                foreach (var product in fragmentIons)
                {
                    output.WriteLine(product.Annotation + "\t" + ((double)glycan.Mass / 1E5 - product.NeutralLoss).ToString() + "\t" + product.NeutralMass.ToString());
                }
            }

        }

        [Test]
        public static void GlyTest_OxoniumIons()
        {
            CommonParameters commonParameters = new CommonParameters();
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/yeast_glycan_25170.mgf");
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();
            //Tips: Using debug mode to check the number of oxoniumIons, in this case will be 7.
            MassDiffAcceptor massDiffAcceptor = new SinglePpmAroundZeroSearchMode(20);
            var oxoinumIonsExist = GlycoPeptides.ScanOxoniumIonFilter(listOfSortedms2Scans[0], massDiffAcceptor);
            Assert.That(oxoinumIonsExist.Where(p => p > 0).Count(), Is.EqualTo(9));
        }

        [Test]
        public static void GlyTest_DistinguishGlycans()
        {
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(H(H(H)))(H(H(H(H)))(H(H(H)))))))", 0).FirstOrDefault();
            Glycan glycan2 = Glycan.Struct2Glycan("(N(N(H(H(H))(H(H(H))(H(H(H(H(H)))))))))", 0).FirstOrDefault();

            var test = Glycan.Equals(glycan, glycan2);
            Assert.That(test, Is.EqualTo(true));

            //TO DO: Test the glycan ions. 
            Glycan glycan3 = Glycan.Struct2Glycan("(N(F)(N(H(H(N(H(N(H(N(H))))))(N(H(N(H(N(F)(H(G))))))))(H(N(H(N(H(N(H(A)))))))(N(F)(H(N(F)(H(N(H)(F))))))))))", 8086).FirstOrDefault();
            Glycan glycan4 = Glycan.Struct2Glycan("(N(F)(N(H(H(N(H(N(H(N(H))))))(N(H(N(H(N(F)(H(A))))))))(H(N(H(N(H(N(H(G)))))))(N(F)(H(N(F)(H(N(H)(F))))))))))", 8087).FirstOrDefault();
        }

        [Test]
        public static void GlyTest_BisectHexNAc()
        {
            //The node here is for check the structure of the glycan. 
            Node node = Glycan.Struct2Node("(N(N(H(N)(H(N)(N))(H(N(H))))))"); //This glycan has a bisect hexnac 
            Assert.That(node.LeftChild.LeftChild.MiddleChild != null);

            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(N)(H(N)(N))(H(N(H))))))", 0).FirstOrDefault();
            Assert.That(glycan.Ions.Count, Is.EqualTo(18));
        }

        [Test]
        public static void GlyTest_GlycanDecoy()
        {
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(N)(H(N)(N))(H(N(H))))))", 0).FirstOrDefault();
            var test = Glycan.BuildTargetDecoyGlycans(new Glycan[] { glycan });
            Assert.That(test.Last().Decoy, Is.EqualTo(true));
            foreach (var ion in test.Last().Ions)
            {
                Assert.That(ion.IonMass + ion.LossIonMass, Is.EqualTo(test.Last().Mass));
            }
        }

        [Test]
        public static void GlyTest_BinarySearch()
        {
            //This is just to test how binary search works.
            double[] array = new double[] { 3.44, 3.45, 4.55, 4.55, 4.55, 4.55, 4.55, 4.55, 4.55, 5.66 };
            double x = 3.43;
            double x1 = 3.44;
            double x2 = 3.441;
            double x3 = 3.45;
            double y = 4.44;
            double z = 5.67;
            double d = 4.55;
            double t = 4.56;
            double t1 = 5.66;
            var xid = GlycoPeptides.BinarySearchGetIndex(array, x);
            var xid1 = GlycoPeptides.BinarySearchGetIndex(array, x1);
            var xid2 = GlycoPeptides.BinarySearchGetIndex(array, x2);
            var xid3 = GlycoPeptides.BinarySearchGetIndex(array, x3);

            var yid = GlycoPeptides.BinarySearchGetIndex(array, y);
            var zid = GlycoPeptides.BinarySearchGetIndex(array, z);
            var did = GlycoPeptides.BinarySearchGetIndex(array, d);
            var tid = GlycoPeptides.BinarySearchGetIndex(array, t);
            var tid1 = GlycoPeptides.BinarySearchGetIndex(array, t1);

            Assert.That(xid, Is.EqualTo(0));
            Assert.That(yid, Is.EqualTo(2));
            Assert.That(zid, Is.EqualTo(10)); //Index out range
            Assert.That(did, Is.EqualTo(2));
            Assert.That(tid, Is.EqualTo(9));
        }         

        [Test]
        public static void GlyTest_NGlycanCompositionFragments()
        {
            var testKind = GlycanDatabase.String2Kind("HexNAc(3)Hex(4)Fuc(2)NeuAc(1)Xylose(1)");

            var ions_NotFucExtended = GlycanDatabase.NGlycanCompositionFragments(testKind);

            var ions_fucExtended = GlycanDatabase.NGlycanCompositionFragments(testKind, true);

            Assert.That(ions_fucExtended.Count >= ions_NotFucExtended.Count);
            Assert.That(ions_NotFucExtended.Count == 35);
            Assert.That(ions_fucExtended.Count == 43);


            var kind = GlycanDatabase.String2Kind("HexNAc(3)Hex(4)Fuc(2)NeuAc(1)");

            Glycan glycan = Glycan.Struct2Glycan("(N(F)(N(H(H)(H(N(F)(H(A)))))))", 0).FirstOrDefault();

            var ionMass = ions_NotFucExtended.Select(p => p.IonMass).ToList();

            var glycanIonmass = glycan.Ions.Select(p => p.IonMass).ToList();

            var overlap = glycanIonmass.Intersect(ionMass).Count();

            Assert.That(overlap == 15);

        }

        [Test]
        public static void LoadCustomMonosaccharides_RoundTripsThroughCompositionAndMass()
        {
            // Register HexA (hexuronic acid, residue mass 176.03209 Da). Verify the loader
            // makes it parseable in composition-format glycan files, that Kind[] gets the right
            // count at the custom slot, that GetMass uses the custom mass, and that
            // GetKindString round-trips the custom code.

            string tsv = string.Join(Environment.NewLine, new[]
            {
                "# comment line ignored",
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "",
                "HexA\tU\t176.03209\t\tHexuronic acid"
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                GlycanDatabase.LoadCustomMonosaccharides(path);

                Assert.That(Glycan.KindCapacity, Is.EqualTo(12)); // 11 built-ins + 1 custom
                Assert.That(Glycan.NameCharDic.ContainsKey("HexA"));
                Assert.That(Glycan.NameCharDic["HexA"].Item1, Is.EqualTo('U'));
                Assert.That(Glycan.NameCharDic["HexA"].Item2, Is.EqualTo(11));
                Assert.That(Glycan.CharMassDic['U'], Is.EqualTo(17603209));

                byte[] kind = GlycanDatabase.String2Kind("HexNAc(2)Hex(5)HexA(1)");
                Assert.That(kind.Length, Is.EqualTo(12));
                Assert.That(kind[0], Is.EqualTo(5));   // Hex
                Assert.That(kind[1], Is.EqualTo(2));   // HexNAc
                Assert.That(kind[11], Is.EqualTo(1));  // HexA

                int expectedMass = 2 * 20307937 + 5 * 16205282 + 1 * 17603209;
                Assert.That(Glycan.GetMass(kind), Is.EqualTo(expectedMass));

                // GetKindString writes built-ins first (H..K) then customs at their slot order.
                Assert.That(Glycan.GetKindString(kind), Is.EqualTo("H5N2U1"));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_CommentsBlanksAndHeaderRowAreSkipped()
        {
            string tsv = string.Join(Environment.NewLine, new[]
            {
                "# only one real entry should be registered",
                "",
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "# another comment",
                "Pent\tT\t132.04226\t\tGeneric pentose",
                ""
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                GlycanDatabase.LoadCustomMonosaccharides(path);

                Assert.That(Glycan.KindCapacity, Is.EqualTo(12));
                Assert.That(Glycan.NameCharDic.ContainsKey("Pent"));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_NameCollisionThrowsWithFileAndLineNumber()
        {
            // "HexNAc" is a built-in; redefining it must fail.
            string tsv = "HexNAc\tU\t999.99999\t\tWill collide";
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                var ex = Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.LoadCustomMonosaccharides(path));
                Assert.That(ex.Message, Does.Contain(Path.GetFileName(path)));
                Assert.That(ex.Message, Does.Contain("line 1"));
                Assert.That(ex.Message, Does.Contain("HexNAc"));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_CharCollisionThrows()
        {
            // 'N' is the built-in single-char code for HexNAc; reusing it must fail.
            string tsv = "Foo\tN\t100.0\t\tCollides with HexNAc";
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                var ex = Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.LoadCustomMonosaccharides(path));
                Assert.That(ex.Message, Does.Contain("'N'"));
                Assert.That(ex.Message, Does.Contain("line 1"));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_NonLetterCodeRejected()
        {
            string tsv = "Foo\t1\t100.0\t\tDigit is not a letter";
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                var ex = Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.LoadCustomMonosaccharides(path));
                Assert.That(ex.Message, Does.Contain("ASCII letter"));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_NonNumericMassThrows()
        {
            string tsv = "Foo\tU\tnot-a-number\t\tBad mass";
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                var ex = Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.LoadCustomMonosaccharides(path));
                Assert.That(ex.Message, Does.Contain("MonoisotopicMass"));
                Assert.That(ex.Message, Does.Contain("not-a-number"));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_DiagnosticIonsEmittedWhenPresent()
        {
            string tsv = "HexA\tU\t176.03209\t175.02482,157.01425\tHexuronic acid";
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                GlycanDatabase.LoadCustomMonosaccharides(path);

                byte[] kind = GlycanDatabase.String2Kind("HexNAc(2)HexA(1)");
                var glycan = new Glycan(null, Glycan.GetMass(kind), kind, null, false, "Nxs", GlycanType.N_glycan);

                int expectedIon1 = (int)Math.Round(175.02482 * 1E5);
                int expectedIon2 = (int)Math.Round(157.01425 * 1E5);
                Assert.That(glycan.GlycanDiagnosticIons.Contains(expectedIon1));
                Assert.That(glycan.GlycanDiagnosticIons.Contains(expectedIon2));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_CustomCodeAcceptedInStructureFormat()
        {
            // Register HexA. A structure-format glycan file using 'U' must now parse without the
            // "Unrecognized character" validator error (which would have fired before registration).

            string monoTsv = "HexA\tU\t176.03209\t\tHexuronic acid";
            string monoPath = Path.GetTempFileName();
            string glycanPath = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(monoPath, monoTsv);
                GlycanDatabase.LoadCustomMonosaccharides(monoPath);

                File.WriteAllText(glycanPath, "(N(H)(U))" + Environment.NewLine);
                var glycans = GlycanDatabase.LoadGlycan(glycanPath, false, false).ToList();

                Assert.That(glycans.Count, Is.EqualTo(2)); // Nxs + Nxt motif duplication
                int expectedMass = 1 * 20307937 + 1 * 16205282 + 1 * 17603209;
                Assert.That(glycans[0].Mass, Is.EqualTo(expectedMass));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(monoPath);
                File.Delete(glycanPath);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_TooFewColumnsThrowsWithFileAndLineNumber()
        {
            // Comments and a valid row first; the bad row at line 4 has only 2 tab-separated
            // columns. The loader requires at least Name / SingleCharCode / MonoisotopicMass.
            string tsv = string.Join(Environment.NewLine, new[]
            {
                "# Custom monosaccharide file with one malformed row",
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "HexA\tU\t176.03209\t\tHexuronic acid",
                "Pent\tT"  // line 4 -- missing MonoisotopicMass column
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                var ex = Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.LoadCustomMonosaccharides(path));

                Assert.That(ex.Message, Does.Contain(Path.GetFileName(path)));
                Assert.That(ex.Message, Does.Contain("line 4"));
                Assert.That(ex.Message, Does.Contain("Expected at least 3 tab-separated columns"));
                // The raw line content should be in the message so the user can find it.
                Assert.That(ex.Message, Does.Contain("Pent"));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_MultiCharCodeRejectedWithFileAndLineNumber()
        {
            // The SingleCharCode column must be exactly one character. A two-character value
            // ("UU") is caught BEFORE Glycan.RegisterCustomMonosaccharide is called, so the
            // error message comes from GlycanDatabase, not from the inner ArgumentException.
            string tsv = string.Join(Environment.NewLine, new[]
            {
                "Name\tSingleCharCode\tMonoisotopicMass",
                "HexA\tU\t176.03209",
                "Pent\tTT\t132.04226"  // line 3 -- two-char code is invalid
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                var ex = Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.LoadCustomMonosaccharides(path));

                Assert.That(ex.Message, Does.Contain(Path.GetFileName(path)));
                Assert.That(ex.Message, Does.Contain("line 3"));
                Assert.That(ex.Message, Does.Contain("SingleCharCode must be exactly one character"));
                Assert.That(ex.Message, Does.Contain("\"TT\""));
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadCustomMonosaccharides_NonNumericDiagnosticIonThrowsWithFileAndLineNumber()
        {
            // The diagnostic-ion column is a comma-separated list of decimal m/z values.
            // A non-numeric entry partway through the list should be caught and reported with
            // the offending entry quoted.
            string tsv = string.Join(Environment.NewLine, new[]
            {
                "# Custom monosaccharide file",
                "",
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "HexA\tU\t176.03209\t175.02482,oops,157.01425\tHas a bad diagnostic ion"  // line 4
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                var ex = Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.LoadCustomMonosaccharides(path));

                Assert.That(ex.Message, Does.Contain(Path.GetFileName(path)));
                Assert.That(ex.Message, Does.Contain("line 4"));
                Assert.That(ex.Message, Does.Contain("DiagnosticIonMasses"));
                Assert.That(ex.Message, Does.Contain("\"oops\""));
                // The two well-formed diagnostic ions surrounding "oops" should not appear in
                // the registry -- the throw must happen before RegisterCustomMonosaccharide.
                Assert.That(Glycan.NameCharDic.ContainsKey("HexA"), Is.False);
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
                File.Delete(path);
            }
        }
    }
}
