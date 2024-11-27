#nullable enable
using EngineLayer;
using Omics;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Proteomics;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace TaskLayer
{

    /// <summary>
    /// Writes a modification annotated database for BioPolymers. 
    /// This class handles the conversion between BioPolymers and the database writer.
    /// </summary>
    public static class PrunedDatabaseWriter
    {
        public static event EventHandler<SingleFileEventArgs> FinishedWritingFileHandler = null!;

        public static async Task WriteDataAsync<T>(string outputPath, List<T> bioPolymers)
            where T : IBioPolymer
        {
            // Simulate writing data to the database
            await Task.Run(() =>
            {
                WriteData(outputPath, bioPolymers);
            });
        }

        public static void WriteData<T>(string outputPath, List<T> bioPolymers, List<string>? nestedIds = null)
            where T: IBioPolymer
        {
            if (bioPolymers.Select(p => p.GetType()).Distinct().Count() != 1)
                throw new ArgumentException("All bioPolymers must be of the same type.");

            switch (bioPolymers.First())
            {
                case Protein:
                    ProteinDbWriter.WriteXmlDatabase([], bioPolymers.Cast<Protein>().ToList(), outputPath);
                    break;
                case RNA:
                    ProteinDbWriter.WriteXmlDatabase([], bioPolymers.Cast<RNA>().ToList(), outputPath);
                    break;
            }

            FinishedWritingFileHandler?.Invoke(outputPath, new SingleFileEventArgs(outputPath, nestedIds));
        }

        /// <summary>
        /// Returns three categories of modification types to write to the database. 
        ///  See <see cref="SearchParameters.ModsToWriteSelection"/>
        /// </summary>
        /// <param name="modTypesToWrite">
        ///     Key is modification type
        ///     Value is integer 0, 1, 2 and 3 interpreted as:
        ///         0:   Do not Write
        ///         1:   Write if in DB and Observed
        ///         2:   Write if in DB
        ///         3:   Write if Observed
        ///   </param>
        /// <returns></returns>
        public static (HashSet<Modification> modificationsToWriteIfBoth, HashSet<Modification> modificationsToWriteIfInDatabase, 
            HashSet<Modification> modificationsToWriteIfObserved) GetModificationsToWrite(Dictionary<string, int> modTypesToWrite)
        {
            var modificationsToWriteIfBoth = new HashSet<Modification>();
            var modificationsToWriteIfInDatabase = new HashSet<Modification>();
            var modificationsToWriteIfObserved = new HashSet<Modification>();

            foreach (var modType in modTypesToWrite)
            {
                foreach (Modification mod in GlobalVariables.AllModsKnown.Where(b => b.ModificationType.Equals(modType.Key)))
                {
                    if (modType.Value == 1) // Write if observed and in database
                    {
                        modificationsToWriteIfBoth.Add(mod);
                    }
                    if (modType.Value == 2) // Write if in database
                    {
                        modificationsToWriteIfInDatabase.Add(mod);
                    }
                    if (modType.Value == 3) // Write if observed
                    {
                        modificationsToWriteIfObserved.Add(mod);
                    }
                }
            }

            return (modificationsToWriteIfBoth, modificationsToWriteIfInDatabase, modificationsToWriteIfObserved);
        }
    }
}
