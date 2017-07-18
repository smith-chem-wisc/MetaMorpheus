using EngineLayer;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace MyBenchmarkk
{
    internal class Program
    {

        #region Private Methods

        private static void Main(string[] args)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Restart();

            string fileName = "C:\\Users\\stepa\\Data\\PaperData\\uniprot-mouse-reviewed-7-13-2017.xml";
            bool generateDecoys = true;
            var Title = GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ?
                "MetaMorpheus: Not a release version" :
                "MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion;
            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalTaskLevelSettings.AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(modFile));
            GlobalTaskLevelSettings.AddMods(GlobalEngineLevelSettings.UniprotDeseralized.OfType<ModificationWithLocation>());
            var ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            var ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            IEnumerable<Modification> localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            bool isContaminant = false;
            List<Protein> proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateDecoys, localizeableModifications, isContaminant, new List<string>(), out Dictionary<string, Modification> um);
            int totalProteins = proteinList.Count;
            Protease protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            int maximumMissedCleavages = 0;
            int? minPeptideLength = 5;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 0;
            int max_mods_for_peptide = 0;
            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;

            Console.WriteLine("Load time: " + stopWatch.Elapsed);
            Console.WriteLine();
            stopWatch.Restart();

            var lockObject = new object();
            int globalSeenUnique = 0;
            Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            {
                int seenUnique = 0;
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            seenUnique++;
                        }
                    }
                }
                lock (lockObject)
                {
                    globalSeenUnique += seenUnique;
                }
            });
            stopWatch.Stop();

            Console.WriteLine("Loop1 time: " + stopWatch.Elapsed);
            Console.WriteLine("globalSeen: " + globalSeenUnique);
            stopWatch.Restart();

            lockObject = new object();
            globalSeenUnique = 0;
            var observed_sequences = new HashSet<string>();
            Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            {
                int seenUnique = 0;
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            var hc = yyy.Sequence;
                            var observed = observed_sequences.Contains(hc);
                            if (observed)
                                continue;
                            lock (observed_sequences)
                            {
                                observed = observed_sequences.Contains(hc);
                                if (observed)
                                    continue;
                                observed_sequences.Add(hc);
                            }
                            seenUnique++;
                        }
                    }
                }
                lock (lockObject)

                {
                    globalSeenUnique += seenUnique;
                }
            });
            stopWatch.Stop();

            Console.WriteLine("Loop2 time: " + stopWatch.Elapsed);
            Console.WriteLine("globalSeenUnique: " + globalSeenUnique);
            stopWatch.Restart();

            lockObject = new object();
            globalSeenUnique = 0;
            var compactPeptides = new HashSet<CompactPeptide>();
            Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            {
                int seenUnique = 0;
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, initiatorMethionineBehavior, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            var hc = yyy.CompactPeptide;
                            var observed = compactPeptides.Contains(hc);
                            if (observed)
                                continue;
                            lock (compactPeptides)
                            {
                                observed = compactPeptides.Contains(hc);
                                if (observed)
                                    continue;
                                compactPeptides.Add(hc);
                            }
                            seenUnique++;
                        }
                    }
                }
                lock (lockObject)

                {
                    globalSeenUnique += seenUnique;
                }
            });
            stopWatch.Stop();

            Console.WriteLine("Loop3 time: " + stopWatch.Elapsed);
            Console.WriteLine("globalSeenUnique: " + globalSeenUnique);
            stopWatch.Restart();

            stopWatch.Stop();
            Console.WriteLine(stopWatch.Elapsed);
        }

        #endregion Private Methods

    }
}