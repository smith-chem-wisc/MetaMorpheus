using EngineLayer;
using Omics;
using Omics.Modifications;
using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    internal class PrunedDatabaseWriter
    {
        // Constructor
        public PrunedDatabaseWriter()
        {

        }

        // Methods
        public async Task WriteDataAsync(string outputPath, List<IBioPolymer> bioPolymers)
        {
            // Simulate writing data to the database
            await Task.Run(() =>
            {
                WriteData(outputPath, bioPolymers);
            });
        }

        public void WriteData(string outputPath, List<IBioPolymer> bioPolymers)
        {
            foreach (var typeGroupedBioPolymers in bioPolymers.GroupBy(p => p.GetType()))
            {
                if (typeGroupedBioPolymers.Key == typeof(Protein))
                    ProteinDbWriter.WriteXmlDatabase([], typeGroupedBioPolymers.Cast<Protein>().ToList(), outputPath);
                else if (typeGroupedBioPolymers.Key == typeof(RNA))
                    ProteinDbWriter.WriteXmlDatabase([], typeGroupedBioPolymers.Cast<RNA>().ToList(), outputPath);
            }
        }

        public static (HashSet<Modification> modificationsToWriteIfBoth, HashSet<Modification> modificationsToWriteIfInDatabase, 
            HashSet<Modification> modificationsToWriteIfObserved) GetModificationsToWrite(Dictionary<string, int> modsToWrite)
        {
            var modificationsToWriteIfBoth = new HashSet<Modification>();
            var modificationsToWriteIfInDatabase = new HashSet<Modification>();
            var modificationsToWriteIfObserved = new HashSet<Modification>();

            foreach (var modType in modsToWrite)
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
