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
        public static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized;
        public static UsefulProteomicsDatabases.Generated.obo psimodDeserialized;
        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

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
            unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
            psimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(psimodLocation);
            uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);

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
            searchModes.Add(new SearchMode("greaterThan-187", (double a) => { return a > -187; }));
            searchModes.Add(new SearchMode("greaterThan-500", (double a) => { return a > -500; }));
            searchModes.Add(new SearchMode("greaterThan-6500", (double a) => { return a > -6500; }));

            exclude = PopulateExcludeList();

            double tolExclude = 0.0025;

            searchModes.Add(new SearchMode("greaterThan-187withExclusions", (double a) =>
            {
                return a > -187 && DoNotExclude(a, tolExclude);
            }));



            searchModes.Add(new SearchMode("greaterThan-500withExclusions", (double a) =>
            {
                return a > -500 && DoNotExclude(a, tolExclude);
            }));
            searchModes.Add(new SearchMode("greaterThan-6500withExclusions", (double a) =>
            {
                return a > -6500 && DoNotExclude(a, tolExclude);
            }));

            searchModes.Add(new SearchMode("500daltonWideSearch", (double a) => { return a > -500 && a < 500; }));
            searchModes.Add(new SearchMode("oldRangeWide", (double a) => { return a > -121.5 && a < 220.5; }));
            searchModes.Add(new SearchMode("oldRangeComb", (double a) => { return a > -121.5 && a < 220.5 && (a % 1 > 0.9 || (a % 1 < 0.2 && a % 1 > -0.1) || a % 1 < -0.8); }));
            searchModes.Add(new SearchMode("withinhalfAdaltonOfZero", (double a) => { return a > -0.5 && a < 0.5; }));
            double tol = 0.0075;
            searchModes.Add(new SearchMode("someNotches", (double a) =>
        {
            return ((a > -0.5 && a < 0.5) ||
                    (a > -17.026549 - tol && a < -17.026549 + tol) ||
                    (a > 0.984016 - tol && a < 0.984016 + tol) ||
                    (a > 15.994915 - tol && a < 15.994915 + tol) ||
                    (a > 21.981943 - tol && a < 21.981943 + tol) ||
                    (a > 31.989829 - tol && a < 31.989829 + tol) ||
                    (a > 57.021464 - tol && a < 57.021464 + tol) ||
                    (a > 79.966331 - tol && a < 79.966331 + tol) ||
                    (a > -18.010565 - tol && a < -18.010565 + tol));
        }));

            var allPsms = new List<PSMwithPeptide>[searchModes.Count];
            for (int j = 0; j < searchModes.Count; j++)
                allPsms[j] = new List<PSMwithPeptide>();

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
                Console.WriteLine("Finished opening spectra file " + origDataFile);

                NewPsm[][] newPsms = SearchEngine.Search(myMsDataFile, spectraFileIndex, peptideIndex, keys, fragmentIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, fragmentTolerance, protease, searchModes);

                AddObservedPeptidesToDictionary(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease);
                fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(compactPeptideToProteinPeptideMatching);

                Console.WriteLine("Finished Search " + origDataFile);

                Console.WriteLine("Doing localization and ion match analysis.. " + origDataFile);

                for (int j = 0; j < searchModes.Count; j++)
                {
                    PSMwithPeptide[] psmsWithPeptides = new PSMwithPeptide[newPsms.Length];

                    Parallel.ForEach(Partitioner.Create(0, newPsms.Length), fff =>
                    {
                        for (int i = fff.Item1; i < fff.Item2; i++)
                        {
                            if (newPsms[i] != null)
                            {
                                var huh = newPsms[i][j];
                                if (huh != null && huh.ScoreFromSearch >= 1)
                                    psmsWithPeptides[i] = new PSMwithPeptide(huh, fullSequenceToProteinSingleMatch[huh.peptide], fragmentTolerance, myMsDataFile);
                            }
                        }
                    });

                    var orderedPsms = psmsWithPeptides.Where(b => b != null).OrderByDescending(b => b.ScoreFromSearch);

                    if (doFDRanalysis)
                    {
                        var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsms);
                        var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.QValue <= 0.01)).ToList();
                        if (limitedpsms_with_fdr.Where(b => !b.isDecoy).Count() > 0)
                            MyAnalysis(limitedpsms_with_fdr, Path.Combine(output_folder, Path.GetFileNameWithoutExtension(origDataFile) + searchModes[j].FileNameAddition + ".mytsv"));
                        WriteToTabDelimitedTextFileWithDecoys(orderedPsmsWithFDR, Path.Combine(output_folder, Path.GetFileNameWithoutExtension(origDataFile) + searchModes[j].FileNameAddition + ".psmtsv"));
                    }
                    else
                    {
                        WriteToTabDelimitedTextFile(orderedPsms, Path.Combine(output_folder, Path.GetFileNameWithoutExtension(origDataFile) + searchModes[j].FileNameAddition + ".psmtsv"));
                    }

                    if (dataFiles.Count > 1)
                    {
                        allPsms[j].AddRange(orderedPsms);
                    }
                }
            }
            if (dataFiles.Count > 1)
            {
                fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(compactPeptideToProteinPeptideMatching);

                for (int j = 0; j < searchModes.Count; j++)
                {
                    var orderedPsms = allPsms[j].OrderByDescending(b => b.ScoreFromSearch).Where(b => b.ScoreFromSearch >= 1).ToList();

                    Reassign(orderedPsms, fullSequenceToProteinSingleMatch);

                    if (doFDRanalysis)
                    {
                        var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsms);
                        var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.QValue <= 0.01)).ToList();
                        if (limitedpsms_with_fdr.Where(b => !b.isDecoy).Count() > 0)
                            MyAnalysis(limitedpsms_with_fdr, Path.Combine(output_folder, "aggregate" + searchModes[j].FileNameAddition + ".mytsv"));
                        WriteToTabDelimitedTextFileWithDecoys(orderedPsmsWithFDR, Path.Combine(output_folder, "aggregate" + searchModes[j].FileNameAddition + ".psmtsv"));
                    }
                    else
                    {
                        WriteToTabDelimitedTextFile(orderedPsms, Path.Combine(output_folder, "aggregate" + searchModes[j].FileNameAddition + ".psmtsv"));
                    }
                }
            }

            Console.WriteLine("All Done!");
            Console.Read();
        }

        private static double[] PopulateExcludeList()
        {
            // Do not exclude Lysine + Anything
            // Do not exclude Lysine, Arginine, Glycine, Asparagine, Alanine, Methionine

            // Exclude -76.134779 and -48.128629 - these are real reversed phosphorylations

            List<double> exclude = new List<double>();

            exclude.Add(-76.134779);
            exclude.Add(-48.128629);
            HashSet<AminoAcid> doNotExcludeEvenCombos = new HashSet<AminoAcid>() { AminoAcid.GetResidue('K') };
            HashSet<AminoAcid> doNotExclude = new HashSet<AminoAcid>() {
                AminoAcid.GetResidue('K'),
                AminoAcid.GetResidue('R'),
                AminoAcid.GetResidue('G'),
                AminoAcid.GetResidue('N'),
                AminoAcid.GetResidue('A'),
                AminoAcid.GetResidue('M'),
            };

            for (char c = 'A'; c <= 'Z'; c++)
            {
                AminoAcid residue;
                if (AminoAcid.TryGetResidue(c, out residue))
                {
                    if (!doNotExclude.Contains(residue))
                        exclude.Add(residue.MonoisotopicMass);
                    for (char cc = 'A'; cc <= 'Z'; cc++)
                    {
                        AminoAcid residueCC;
                        if (AminoAcid.TryGetResidue(cc, out residueCC))
                        {
                            if (!doNotExcludeEvenCombos.Contains(residueCC))
                                exclude.Add(residue.MonoisotopicMass + residueCC.MonoisotopicMass);
                        }
                    }
                }
            }
            return exclude.ToArray();
        }

        private static double[] exclude;
        private static bool DoNotExclude(double a, double tolExclude)
        {

            foreach (var heh in exclude)
                if (Math.Abs(heh - a) > tolExclude)
                    return false;
            return true;

        }

        private static void Reassign(List<PSMwithPeptide> orderedPsms, Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch)
        {
            foreach (var huhuh in orderedPsms)
            {
                huhuh.Reassign(fullSequenceToProteinSingleMatch);
            }
        }

        private static void AddObservedPeptidesToDictionary(NewPsm[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease)
        {
            foreach (var ah in newPsms)
            {
                if (ah != null)
                    foreach (var fhh in ah)
                    {
                        if (fhh != null && !fullSequenceToProteinPeptideMatching.ContainsKey(fhh.peptide))
                            fullSequenceToProteinPeptideMatching.Add(fhh.peptide, new HashSet<PeptideWithSetModifications>());
                    }
            }

            foreach (var protein in proteinList)
                foreach (var peptide in protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable).ToList())
                {
                    if (peptide.Length == 1 || peptide.Length > 252)
                        continue;
                    peptide.SetFixedModifications(fixedModifications);
                    var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, 4098, 3, localizeableModifications).ToList();
                    foreach (var yyy in ListOfModifiedPeptides)
                    {
                        HashSet<PeptideWithSetModifications> v;
                        if (fullSequenceToProteinPeptideMatching.TryGetValue(new CompactPeptide(yyy, variableModifications, localizeableModifications), out v))
                        {
                            v.Add(yyy);
                        }
                    }
                }
        }

        private static Dictionary<CompactPeptide, PeptideWithSetModifications> GetSingleMatchDictionary(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // Right now very stupid, add the first decoy one, and if no decoy, add the first one
            Dictionary<CompactPeptide, PeptideWithSetModifications> outDict = new Dictionary<CompactPeptide, PeptideWithSetModifications>();
            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                bool sawDecoy = false;
                foreach (var entry in kvp.Value)
                {
                    if (entry.protein.isDecoy)
                    {
                        outDict[kvp.Key] = entry;
                        sawDecoy = true;
                        break;
                    }
                }
                if (sawDecoy == false)
                    outDict[kvp.Key] = kvp.Value.First();
            }
            return outDict;
        }

        internal static byte PeptideTypeFromString(string peptideDescription)
        {
            switch (peptideDescription)
            {
                case "full":
                    return 0;

                case "full:M cleaved":
                    return 1;

                case "peptide start":
                    return 2;

                case "peptide end":
                    return 3;

                case "propeptide start":
                    return 4;

                case "propeptide end":
                    return 5;

                case "chain start":
                    return 6;

                case "chain end":
                    return 7;
            }
            return byte.MaxValue;
        }

        private static string GetStringDescriptionFromByte(byte peptideType)
        {
            switch (peptideType)
            {
                case 0:
                    return "full";

                case 1:
                    return "full: M cleaved";

                case 2:
                    return "peptide start";

                case 3:
                    return "peptide end";

                case 4:
                    return "propeptide start";

                case 5:
                    return "propeptide end";

                case 6:
                    return "chain start";

                case 7:
                    return "chain end";
            }
            return null;
        }

        private static Protein getProteinFromOffset(List<Protein> hm, int key)
        {
            var loc = Array.BinarySearch(hm.Select(b => b.offset).ToArray(), key);
            if (loc < 0)
                loc = ~loc;
            if (loc <= 0)
                return hm[0];
            if (loc >= hm.Count)
                return hm.Last();
            if (key >= hm[loc].offset)
                return hm[loc];
            return hm[loc - 1];
        }

        private static void MyAnalysis(List<NewPsmWithFDR> limitedpsms_with_fdr, string filepath)
        {
            MyNewTreeStructure myTreeStructure = new MyNewTreeStructure();
            myTreeStructure.GenerateBins(limitedpsms_with_fdr, 0.003);
            myTreeStructure.AddToBins(limitedpsms_with_fdr);

            Console.WriteLine("Identifying bins...");
            MyAnalysisClass.IdentifyUnimodBins(myTreeStructure, 0.003);
            MyAnalysisClass.IdentifyUniprotBins(myTreeStructure, 0.003);
            MyAnalysisClass.IdentifyAA(myTreeStructure, 0.003);

            Console.WriteLine("Identifying combos...");
            MyAnalysisClass.IdentifyCombos(myTreeStructure, 0.003);

            Console.WriteLine("Extracting residues from localizeable...");
            MyAnalysisClass.IdentifyResidues(myTreeStructure);

            Console.WriteLine("Identifying mods in common...");
            MyAnalysisClass.IdentifyMods(myTreeStructure);

            Console.WriteLine("Identifying AAs in common...");
            MyAnalysisClass.IdentifyAAsInCommon(myTreeStructure);

            Console.WriteLine("Identifying mine...");
            MyAnalysisClass.IdentifyMine(myTreeStructure, 0.003);

            Console.WriteLine("Writing...");
            using (StreamWriter output = new StreamWriter(filepath))
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
                        + "\t" + (bin.Count == 0 ? double.NaN : bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + b.Value / bin.CountTarget))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + b.Value / bin.CountTarget))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : bin.NlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : bin.ClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }

            Console.WriteLine("Done with my analysis analysis.. ");
        }

        private static List<NewPsmWithFDR> DoFalseDiscoveryRateAnalysis(IEnumerable<PSMwithPeptide> items)
        {
            List<NewPsmWithFDR> ids = new List<NewPsmWithFDR>();

            int cumulative_target = 0;
            int cumulative_decoy = 0;
            foreach (PSMwithPeptide item in items)
            {
                var isDecoy = item.isDecoy;
                if (isDecoy)
                    cumulative_decoy++;
                else
                    cumulative_target++;
                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                ids.Add(new NewPsmWithFDR(item, cumulative_target, cumulative_decoy, temp_q_value));
            }

            double min_q_value = double.PositiveInfinity;
            for (int i = ids.Count - 1; i >= 0; i--)
            {
                NewPsmWithFDR id = ids[i];
                if (id.QValue > min_q_value)
                    id.QValue = min_q_value;
                else if (id.QValue < min_q_value)
                    min_q_value = id.QValue;
            }

            return ids;
        }

        public static void WriteToTabDelimitedTextFileWithDecoys(List<NewPsmWithFDR> items, string filepath)
        {
            Console.WriteLine("Writing psms");
            using (StreamWriter output = new StreamWriter(filepath))
            {
                output.WriteLine("Spectrum File\tScan Number\tRetention Time\tPrecursor MZ\tPrecursor Charge\tPrecursor Intensity\tExperimental Peaks\tTotal Intensity\tPrecursor Mass\tScoreFromSearch\tPreviousAminoAcid\tSequence\tNextAminoAcid\tnumVariableMods\tStart Residue\tEnd Residue\tPeptide\tMissed Cleavages\tPeptide Mass\tProtein\tMass Diff(Da)\tMatched Fragments\tMatched Counts\tLocalized Scores\tImprovement\tImprovment Residue\tImprovement Terminus\tDecoy\tCumulative Target\tCumulative Decoy\tQ-value");
                for (int i = 0; i < items.Count; i++)
                    output.WriteLine(items[i].ToString());
            }
        }

        public static void WriteToTabDelimitedTextFile(IEnumerable<PSMwithPeptide> items, string filepath)
        {
            using (StreamWriter output = new StreamWriter(filepath))
            {
                output.WriteLine("Spectrum File\tScan Number\tRetention Time\tPrecursor MZ\tPrecursor Charge\tPrecursor Intensity\tExperimental Peaks\tTotal Intensity\tPrecursor Mass\tScoreFromSearch\tPreviousAminoAcid\tSequence\tNextAminoAcid\tnumVariableMods\tStart Residue\tEnd Residue\tPeptide\tMissed Cleavages\tPeptide Mass\tProtein\tMass Diff(Da)\tMatched Fragments\tMatched Counts\tLocalized Scores\tImprovement\tImprovment Residue\tImprovement Terminus\tDecoy");
                foreach (var i in items)
                    output.WriteLine(i.ToString());
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