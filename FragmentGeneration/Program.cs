using IndexSearchAndAnalyze;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MetaMorpheus;
using Proteomics;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace FragmentGeneration
{
    internal class Program
    {

        public static string unimodLocation = @"unimod_tables.xml";
        public static string psimodLocation = @"PSI-MOD.obo.xml";
        public static string elementsLocation = @"elements.dat";
        public static string uniprotLocation = @"ptmlist.txt";

        private static ObservableCollection<XMLdb> xMLdblist = new ObservableCollection<XMLdb>();
        private static ObservableCollection<RawDataAndResults> rawDataAndResultslist = new ObservableCollection<RawDataAndResults>();

        private static void Main(string[] args)
        {
            Console.WriteLine("Loading amino acid masses...");

            Console.WriteLine("Starting program");

            Console.WriteLine("Loading modification and element databases...");
            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);

            var unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
            var uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);

            Console.WriteLine("Loading amino acid masses...");
            AminoAcidMasses.LoadAminoAcidMasses();

            if (args[0].Equals("mouse"))
            {
                xMLdblist.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-mouse-reviewed-12-23-2016.xml"));
            }
            else
            {
                xMLdblist.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-human-reviewed-12-15-2016.xml"));
            }

            //xMLdblist.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-622-12-15-2016.xml"));
            //xMLdblist.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-P62263-12-16-2016.xml"));
            //xMLdblist.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-O60832-12-20-2016.xml"));

            ObservableCollection<ModList> collectionOfModLists = new ObservableCollection<ModList>();

            collectionOfModLists.Add(new ModList("p.txt", true, false, false));
            collectionOfModLists.Add(new ModList("v.txt", false, true, false));
            collectionOfModLists.Add(new ModList("f.txt", false, false, true));

            bool doFDRanalysis = true;

            Console.WriteLine("Reading modifications...");
            List<MorpheusModification> variableModifications = collectionOfModLists.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = collectionOfModLists.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = collectionOfModLists.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            GenerateModsFromStrings(xMLdblist.Select(b => b.FileName).ToList(), localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            Console.WriteLine("Unknown modifications (can add these to the localize list):");
            foreach (var hadsfm in unidentifiedModStrings)
                Console.WriteLine("\t" + hadsfm);

            Console.WriteLine("loading proteins...");
            var proteinList = xMLdblist.SelectMany(b => b.getProteins(doFDRanalysis, identifiedModsInXML)).ToList();
            Console.WriteLine("loaded!");
            Console.WriteLine("Proteins: " + proteinList.Count());
            Console.WriteLine("Max Mods: " + proteinList.Select(b => b.OneBasedPossibleLocalizedModifications.Count).Max());
            Console.WriteLine("Total Amino Acids: " + proteinList.Select(b => b.Length).Sum());

            var protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];

            List<CompactPeptide> peptideIndex;
            Dictionary<float, List<int>> fragmentIndexDict;

            Indices.GetPeptideAndFragmentIndices(out peptideIndex, out fragmentIndexDict, xMLdblist, collectionOfModLists, doFDRanalysis, variableModifications, fixedModifications, localizeableModifications, proteinList, protease);

            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

            List<SearchMode> searchModes = new List<SearchMode>();

            searchModes.Add(new SearchMode("open", (double a) => { return true; }));
            //searchModes.Add(new SearchMode("greaterThan-187", (double a) => { return a > -187; }));
            //searchModes.Add(new SearchMode("greaterThan-500", (double a) => { return a > -500; }));
            //searchModes.Add(new SearchMode("greaterThan-6500", (double a) => { return a > -6500; }));

            double[] exclude = Exclusions.PopulateExcludeList();

            double tolExclude = 0.0025;

            searchModes.Add(new SearchMode("greaterThan-187withExclusions", (double a) =>
            {
                return a > -187 && Exclusions.DoNotExclude(a, tolExclude, exclude);
            }));



            //    searchModes.Add(new SearchMode("greaterThan-500withExclusions", (double a) =>
            //    {
            //        return a > -500 && DoNotExclude(a, tolExclude);
            //    }));
            //    searchModes.Add(new SearchMode("greaterThan-6500withExclusions", (double a) =>
            //    {
            //        return a > -6500 && DoNotExclude(a, tolExclude);
            //    }));

            //    searchModes.Add(new SearchMode("500daltonWideSearch", (double a) => { return a > -500 && a < 500; }));
            //    searchModes.Add(new SearchMode("oldRangeWide", (double a) => { return a > -121.5 && a < 220.5; }));
            //    searchModes.Add(new SearchMode("oldRangeComb", (double a) => { return a > -121.5 && a < 220.5 && (a % 1 > 0.9 || (a % 1 < 0.2 && a % 1 > -0.1) || a % 1 < -0.8); }));
            searchModes.Add(new SearchMode("withinhalfAdaltonOfZero", (double a) => { return a > -0.5 && a < 0.5; }));
            //    double tol = 0.0075;
            //    searchModes.Add(new SearchMode("someNotches", (double a) =>
            //{
            //    return ((a > -0.5 && a < 0.5) ||
            //            (a > -17.026549 - tol && a < -17.026549 + tol) ||
            //            (a > 0.984016 - tol && a < 0.984016 + tol) ||
            //            (a > 15.994915 - tol && a < 15.994915 + tol) ||
            //            (a > 21.981943 - tol && a < 21.981943 + tol) ||
            //            (a > 31.989829 - tol && a < 31.989829 + tol) ||
            //            (a > 57.021464 - tol && a < 57.021464 + tol) ||
            //            (a > 79.966331 - tol && a < 79.966331 + tol) ||
            //            (a > -18.010565 - tol && a < -18.010565 + tol));
            //}));

            var allPsms = new List<NewPsm>[searchModes.Count];
            for (int j = 0; j < searchModes.Count; j++)
                allPsms[j] = new List<NewPsm>();

            List<string> dataFiles;

            if (args[0].Equals("mouse"))
            {
                dataFiles = new List<string>() {              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac1_9uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac2_9p5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac3_9p5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac4_8uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac5_4uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac6_5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac7_5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac8_9p5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac1_9uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac2_9uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac3_6uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac4_6uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac5_4uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac6_5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac7_6uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac8_9p5uL-Calibrated.mzML" ,
                                                          @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-30-13_CAST_Frac9_9p5uL-Calibrated.mzML"  };
            }
            else
            {
                dataFiles = new List<string>() { @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac1-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac2-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac3-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac4-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac5-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac6-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac7_120430090151-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac8_120430121912-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac9-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac10-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac11-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac12-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac13-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac14-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac15-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac16-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac17-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac18-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac19-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac20-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac21-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac22-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac23-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac24-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac25-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac26-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac27-Calibrated.mzML" ,
                                                              @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac28-Calibrated.mzML" };
            }

            foreach (var fileName in dataFiles)
                rawDataAndResultslist.Add(new RawDataAndResults(fileName, null, null));

            string output_folder = Path.Combine(Path.GetDirectoryName(dataFiles[0]), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));

            if (!Directory.Exists(output_folder))
                Directory.CreateDirectory(output_folder);

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();

            double fragmentTolerance = 0.01;

            for (int spectraFileIndex = 0; spectraFileIndex < dataFiles.Count; spectraFileIndex++)
            {
                var origDataFile = dataFiles[spectraFileIndex];
                Console.WriteLine("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile);
                Console.WriteLine("Opening spectra file...");
                myMsDataFile.Open();
                Console.WriteLine("Finished opening spectra file " + Path.GetFileName(origDataFile));

                SearchParams searchParams = new SearchParams(myMsDataFile, spectraFileIndex, peptideIndex, keys, fragmentIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, fragmentTolerance, protease, searchModes);
                SearchEngine searchEngine = new SearchEngine(searchParams);
                SearchResults searchResults = (SearchResults)searchEngine.Run();
                List<NewPsm>[] newPsms = searchResults.newPsms;

                Console.WriteLine(searchResults);

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(newPsms[i]);

                AnalysisParams analysisParams = new AnalysisParams(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, fragmentTolerance, unimodDeserialized, uniprotDeseralized, (MyNewTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => WriteToTabDelimitedTextFileWithDecoys(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s));
                AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);
                AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

                Console.WriteLine(analysisResults);

            }
            if (dataFiles.Count > 1)
            {

                AnalysisParams analysisParams = new AnalysisParams(allPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, null, fragmentTolerance, unimodDeserialized, uniprotDeseralized, (MyNewTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, "aggregate" + s), (List<NewPsmWithFDR> h, string s) => WriteToTabDelimitedTextFileWithDecoys(h, output_folder, "aggregate" + s));
                AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);
                AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();
            }

            Console.WriteLine("All Done!");
            Console.Read();
        }

        private static void WriteTree(MyNewTreeStructure myTreeStructure, string output_folder, string fileName)
        {
            using (StreamWriter output = new StreamWriter(Path.Combine(output_folder, fileName + ".mytsv")))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tNtermLocFrac\tCtermLocFrac\tUniprot");
                foreach (Bin bin in myTreeStructure.finalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + (double)b.Value / bin.CountTarget))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + (double)b.Value / bin.CountTarget))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.NlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }
        }



        public static void WriteToTabDelimitedTextFileWithDecoys(List<NewPsmWithFDR> items, string output_folder, string fileName)
        {
            Console.WriteLine("Writing psms");
            using (StreamWriter output = new StreamWriter(Path.Combine(output_folder, fileName + ".psmtsv")))
            {
                output.WriteLine("Spectrum File\tScan Number\tRetention Time\tPrecursor MZ\tPrecursor Charge\tPrecursor Intensity\tExperimental Peaks\tTotal Intensity\tPrecursor Mass\tScoreFromSearch\tPreviousAminoAcid\tSequence\tNextAminoAcid\tnumVariableMods\tStart Residue\tEnd Residue\tPeptide\tMissed Cleavages\tPeptide Mass\tProtein\tMass Diff(Da)\tMatched Fragments\tMatched Counts\tLocalized Scores\tImprovement\tImprovment Residue\tImprovement Terminus\tDecoy\tCumulative Target\tCumulative Decoy\tQ-value");
                for (int i = 0; i < items.Count; i++)
                    output.WriteLine(items[i].ToString());
            }
        }

        private static void GenerateModsFromStrings(List<string> listOfXMLdbs, List<MorpheusModification> modsKnown, out Dictionary<string, List<MorpheusModification>> modsToLocalize, out HashSet<string> modsInXMLtoTrim)
        {
            modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
            var modsInXML = ProteomeDatabaseReader.ReadXMLmodifications(listOfXMLdbs);
            modsInXMLtoTrim = new HashSet<string>(modsInXML);
            foreach (var knownMod in modsKnown)
                if (modsInXML.Contains(knownMod.NameInXML))
                {
                    if (modsToLocalize.ContainsKey(knownMod.NameInXML))
                        modsToLocalize[knownMod.NameInXML].Add(knownMod);
                    else
                        modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification>() { knownMod });
                    modsInXMLtoTrim.Remove(knownMod.NameInXML);
                }
        }
    }
}