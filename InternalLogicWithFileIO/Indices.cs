using InternalLogic;
using OldInternalLogic;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;

namespace InternalLogicWithFileIO
{
    internal class Indices
    {
        public static void GetPeptideAndFragmentIndices(out List<CompactPeptide> peptideIndex, out Dictionary<float, List<int>> fragmentIndexDict, List<ModListForSearch> collectionOfModLists, bool doFDRanalysis, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, List<Protein> hm, Protease protease, AllTasksParams po, string output_folder)
        {
            #region Index file names

            string folderName = output_folder;
            StringBuilder indexFileSB = new StringBuilder();
            foreach (var heh in po.xMLdblist)
                indexFileSB.Append(Path.GetFileNameWithoutExtension(heh));
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
                po.output("Generating indices...");

                IndexParams indexParams = new IndexParams(hm, variableModifications, fixedModifications, localizeableModifications, protease, po);
                IndexEngine indexEngine = new IndexEngine(indexParams);
                IndexResults indexResults = (IndexResults)indexEngine.Run();
                peptideIndex = indexResults.peptideIndex;
                fragmentIndexDict = indexResults.fragmentIndexDict;

                po.output("Writing peptide index...");
                writePeptideIndex(peptideIndex, peptideIndexFile);
                po.output("Writing fragment index...");
                writeFragmentIndexNetSerializer(fragmentIndexDict, fragmentIndexFile, po);
                po.output("Done Writing fragment index");
            }
            else
            {
                po.output("Reading peptide index...");
                peptideIndex = readPeptideIndex(peptideIndexFile);
                po.output("Reading fragment index...");
                fragmentIndexDict = readFragmentIndexNetSerializer(fragmentIndexFile, po);
            }
        }

        private static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            foreach (var ok in type.Assembly.GetTypes().Where(t => t.IsSubclassOf(type)))
                yield return ok;
            yield return type;
        }

        internal static Dictionary<float, List<int>> readFragmentIndexNetSerializer(string fragmentIndexFile, AllTasksParams po)
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
            po.output("Time to read fragment index with netSerializer: " + elapsedTime);

            return newPerson;
        }

        internal static void writeFragmentIndexNetSerializer(Dictionary<float, List<int>> fragmentIndex, string fragmentIndexFile, AllTasksParams po)
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
            po.output("Time to write fragment index with netserializer: " + elapsedTime);
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