using IndexSearchAndAnalyze;
using MetaMorpheus;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;

namespace FragmentGeneration
{
    internal class Indices
    {
        public static void GetPeptideAndFragmentIndices(out List<CompactPeptide> peptideIndex, out Dictionary<float, List<int>> fragmentIndexDict, ObservableCollection<XMLdb> xMLdblist, ObservableCollection<ModList> collectionOfModLists, bool doFDRanalysis, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, List<Protein> hm, Protease protease)
        {
            #region Index file names

            string folderName = Path.GetDirectoryName(xMLdblist.First().FileName);
            StringBuilder indexFileSB = new StringBuilder();
            foreach (var heh in xMLdblist)
                indexFileSB.Append(Path.GetFileNameWithoutExtension(heh.FileName));
            if (doFDRanalysis)
                indexFileSB.Append("-WithDecoys");
            if (collectionOfModLists.Where(b => b.Fixed).Count() > 0)
            {
                indexFileSB.Append("-fixed");
                foreach (var heh in collectionOfModLists.Where(b => b.Fixed))
                    indexFileSB.Append("-" + Path.GetFileNameWithoutExtension(heh.FileName));
            }
            if (collectionOfModLists.Where(b => b.Fixed).Count() > 0)
            {
                indexFileSB.Append("-variable");
                foreach (var heh in collectionOfModLists.Where(b => b.Variable))
                    indexFileSB.Append("-" + Path.GetFileNameWithoutExtension(heh.FileName));
            }
            if (collectionOfModLists.Where(b => b.Localize).Count() > 0)
            {
                indexFileSB.Append("-localize");
                foreach (var heh in collectionOfModLists.Where(b => b.Localize))
                    indexFileSB.Append("-" + Path.GetFileNameWithoutExtension(heh.FileName));
            }

            string peptideIndexFile = Path.Combine(folderName, indexFileSB.ToString() + "-peptideIndex.ind");
            string fragmentIndexFile = Path.Combine(folderName, indexFileSB.ToString() + "-fragmentIndex.ind");

            #endregion Index file names

            if (!File.Exists(peptideIndexFile) || !File.Exists(fragmentIndexFile))
            {
                Console.WriteLine("Generating indices...");

                //IndexParams indexParams = new IndexParams(hm, variableModifications, fixedModifications, localizeableModifications, protease);
                //IndexEngine indexEngine = new IndexEngine(indexParams);
                //IndexResults indexResults = (IndexResults)indexEngine.Run();
                //peptideIndex = indexResults.peptideIndex;
                //fragmentIndexDict = indexResults.fragmentIndexDict;
                peptideIndex = null;
                fragmentIndexDict = null;

                //Console.WriteLine("Writing peptide index...");
                //writePeptideIndex(peptideIndex, peptideIndexFile);
                //Console.WriteLine("Writing fragment index...");
                //writeFragmentIndexNetSerializer(fragmentIndexDict, fragmentIndexFile);
                //Console.WriteLine("Done Writing fragment index");
            }
            else
            {
                Console.WriteLine("Reading peptide index...");
                peptideIndex = readPeptideIndex(peptideIndexFile);
                Console.WriteLine("Reading fragment index...");
                fragmentIndexDict = readFragmentIndexNetSerializer(fragmentIndexFile);
            }
        }

        private static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            foreach (var ok in type.Assembly.GetTypes().Where(t => t.IsSubclassOf(type)))
                yield return ok;
            yield return type;
        }

        internal static Dictionary<float, List<int>> readFragmentIndexNetSerializer(string fragmentIndexFile)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();

            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            Dictionary<float, List<int>> newPerson;
            using (var file = File.OpenRead(fragmentIndexFile))
                newPerson = (Dictionary<float, List<int>>)ser.Deserialize(file);

            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = string.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                ts.Hours, ts.Minutes, ts.Seconds,
                ts.Milliseconds / 10);
            Console.WriteLine("Time to read fragment index with netSerializer: " + elapsedTime);

            return newPerson;
        }

        internal static void writeFragmentIndexNetSerializer(Dictionary<float, List<int>> fragmentIndex, string fragmentIndexFile)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();

            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);

            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = string.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                ts.Hours, ts.Minutes, ts.Seconds,
                ts.Milliseconds / 10);
            Console.WriteLine("Time to write fragment index with netserializer: " + elapsedTime);
        }

        internal static void writePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFile))
            {
                ser.Serialize(file, peptideIndex);
            }
        }

        internal static List<CompactPeptide> readPeptideIndex(string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);
            List<CompactPeptide> newPerson;
            using (var file = File.OpenRead(peptideIndexFile))
            {
                newPerson = (List<CompactPeptide>)ser.Deserialize(file);
            }

            return newPerson;
        }
    }
}