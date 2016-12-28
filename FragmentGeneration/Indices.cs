using MetaMorpheus;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
                GenerateIndices(hm, variableModifications, fixedModifications, out peptideIndex, out fragmentIndexDict, localizeableModifications, protease);
                Console.WriteLine("Writing peptide index...");
                writePeptideIndex(peptideIndex, peptideIndexFile);
                Console.WriteLine("Writing fragment index...");
                writeFragmentIndexNetSerializer(fragmentIndexDict, fragmentIndexFile);
                Console.WriteLine("Done Writing fragment index");
            }
            else
            {
                Console.WriteLine("Reading peptide index...");
                peptideIndex = readPeptideIndex(peptideIndexFile);
                Console.WriteLine("Reading fragment index...");
                fragmentIndexDict = readFragmentIndexNetSerializer(fragmentIndexFile);
            }
        }

        internal static void GenerateIndices(List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, out List<CompactPeptide> myDictionaryOut, out Dictionary<float, List<int>> myFragmentDictionaryOut, List<MorpheusModification> localizeableModifications, Protease protease)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();

            var myDictionary = new List<CompactPeptide>();
            var myFragmentDictionary = new Dictionary<float, List<int>>(100000);
            int numProteins = 0;
            int totalProteins = proteinList.Count();
            HashSet<string> level3_observed = new HashSet<string>();
            HashSet<string> level4_observed = new HashSet<string>();

            var lp = new List<ProductType>() { ProductType.b, ProductType.y };
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<float, List<int>> myInnerDictionary = new Dictionary<float, List<int>>(100000);
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length == 1 || peptide.Length > 252)
                            continue;

                        if (peptide.OneBasedPossibleLocalizedModifications.Count == 0)
                        {
                            lock (level3_observed)
                            {
                                var hc = peptide.BaseLeucineSequence;
                                var observed = level3_observed.Contains(hc);
                                if (observed)
                                    continue;
                                level3_observed.Add(hc);
                            }
                        }

                        peptide.SetFixedModifications(fixedModifications);

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, 4098, 3, localizeableModifications).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            if (peptide.OneBasedPossibleLocalizedModifications.Count > 0)
                            {
                                lock (level4_observed)
                                {
                                    var hc = yyy.Sequence;
                                    var observed = level4_observed.Contains(hc);
                                    if (observed)
                                        continue;
                                    level4_observed.Add(hc);
                                }
                            }

                            var ps = new CompactPeptide(yyy, variableModifications, localizeableModifications);

                            //var ps = new CompactPeptide(yyy.OneBasedStartResidueInProtein + yyy.protein.offset,
                            //                                          (byte)yyy.Length,
                            //                                          varMod1Type,
                            //                                          varMod1Loc,
                            //                                          varMod2Type,
                            //                                          varMod2Loc,
                            //                                          varMod3Type,
                            //                                          varMod3Loc,
                            //                                          protein.isDecoy,
                            //                                          (float)yyy.MonoisotopicMass,
                            //                                          Program.PeptideTypeFromString(yyy.PeptideDescription));

                            //Console.WriteLine(yyy.Sequence);
                            int index;
                            lock (myDictionary)
                            {
                                index = myDictionary.Count;
                                myDictionary.Add(ps);
                            }

                            foreach (var huhu in yyy.FastUnsortedProductMasses(lp))
                            {
                                //if (yyy.BaseSequence.Equals("RL"))
                                //Console.WriteLine("Adding " + index + " to " + huhu);
                                float rounded = (float)Math.Round(huhu, 3);
                                List<int> value;
                                if (myInnerDictionary.TryGetValue(rounded, out value))
                                    value.Add(index);
                                else
                                    myInnerDictionary.Add(rounded, new List<int>() { index });
                            }
                        }
                    }
                    numProteins++;
                    if (numProteins % 100 == 0)
                        Console.WriteLine("Proteins: " + numProteins + " / " + totalProteins);
                }
                lock (myFragmentDictionary)
                {
                    foreach (var huhu in myInnerDictionary)
                    {
                        List<int> value;
                        foreach (var hhhh in huhu.Value)
                        {
                            if (myFragmentDictionary.TryGetValue(huhu.Key, out value))
                                value.Add(hhhh);
                            else
                                myFragmentDictionary.Add(huhu.Key, new List<int>() { hhhh });
                        }
                    }
                }
            });

            //bool[] vetted = new bool[myDictionary.Count];
            //Dictionary<int, HashSet<int>> AmbiguityGroups = new Dictionary<int, HashSet<int>>();
            //foreach (var djfkjdf in myFragmentDictionary.OrderByDescending(b => b.Key))
            //{
            //    //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //        //Console.WriteLine("Considering " + djfkjdf.Key);
            //    var mass = djfkjdf.Key;
            //    var listOfPeptides = djfkjdf.Value;

            //    if (listOfPeptides.Count == 1)
            //    {
            //        //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //            //Console.WriteLine("vetting  " + listOfPeptides.First());
            //        vetted[listOfPeptides.First()] = true;
            //        continue;
            //    }

            //    HashSet<int> nonVettedGroup = new HashSet<int>();
            //    // At least two peptides for this fragment... ambiguous!
            //    foreach (var pep in listOfPeptides)
            //    {
            //        if (vetted[pep])
            //            continue;

            //        // Non-vetted peptide
            //        nonVettedGroup.Add(pep);
            //    }
            //    if (nonVettedGroup.Count == 0)
            //    {
            //        continue;
            //    }
            //    if (nonVettedGroup.Count == 1)
            //    {
            //        //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //            //Console.WriteLine("vetting  " + nonVettedGroup.First());
            //        vetted[nonVettedGroup.First()] = true;
            //        continue;
            //    }
            //    // At least two non-vetted peptides here!

            //    //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //    //    Console.WriteLine("At least two non-vetted peptides here! " + string.Join(",", nonVettedGroup));

            //    // At this point, want to add ALL the ambiguity groups:
            //    // For every peptide in nonVettedGroup. But nothing else.
            //    //
            //    // At this point, want to trim ambiguity groups:
            //    //
            //    // Nonvetted group = [1,2,3]
            //    // Case A:
            //    // Existing ambiguity group of 1: Do Not Exist. -> Add 2 and 3 to the ambiguity group
            //    // Existing ambiguity group of 1: [2] -> do nothing
            //    // Existing ambiguity group of 1: [4] -> Make empty and vet for 1! Also, remove 1 from ambiguity group of 4.
            //    foreach (var pep in nonVettedGroup)
            //    {
            //        //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //        //    Console.WriteLine("  Considering " + pep);

            //        HashSet<int> AmbiguityGroup;
            //        if (AmbiguityGroups.TryGetValue(pep, out AmbiguityGroup))
            //        {
            //            //foreach (var sdfdsf in AmbiguityGroup)
            //            //{
            //            //    //Console.WriteLine("   Considering ambiguity group of " + sdfdsf);
            //            //    var ag = AmbiguityGroups[sdfdsf];
            //            //    if (nonVettedGroup.Contains(sdfdsf))
            //            //    {
            //            //        if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //            //            Console.WriteLine("   nonVettedGroup contains " + sdfdsf);
            //            //    }
            //            //    else
            //            //    {
            //            //        if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //            //            Console.WriteLine("   nonVettedGroup does not contain " + sdfdsf);
            //            //    }

            //            //}
            //            //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //            //    Console.WriteLine("Ambiguity group was: " + string.Join(",", AmbiguityGroup));
            //            AmbiguityGroup.IntersectWith(nonVettedGroup);
            //            //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //            //    Console.WriteLine("Ambiguity group now: " + string.Join(",", AmbiguityGroup));
            //            if (AmbiguityGroup.Count == 0)
            //            {
            //                vetted[pep] = true;
            //                //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //                //    Console.WriteLine("vetting  " + pep + " because now AmbiguityGroup of " + pep + " has no members");
            //            }
            //        }
            //    }

            //    nonVettedGroup = new HashSet<int>();
            //    // At least two peptides for this fragment... ambiguous!
            //    foreach (var pep in listOfPeptides)
            //    {
            //        if (vetted[pep])
            //            continue;

            //        // Non-vetted peptide
            //        nonVettedGroup.Add(pep);
            //    }
            //    if (nonVettedGroup.Count == 0)
            //    {
            //        continue;
            //    }
            //    if (nonVettedGroup.Count == 1)
            //    {
            //      //  if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //        //    Console.WriteLine("vetting  " + nonVettedGroup.First());
            //        vetted[nonVettedGroup.First()] = true;
            //        continue;
            //    }

            //    foreach (var pep in nonVettedGroup)
            //    {
            //        HashSet<int> AmbiguityGroup;
            //        if (!AmbiguityGroups.TryGetValue(pep, out AmbiguityGroup))
            //        {
            //            AmbiguityGroups.Add(pep, new HashSet<int>(nonVettedGroup.Where(b => b != pep)));
            //            //if (djfkjdf.Value.Contains(113) || djfkjdf.Value.Contains(114))
            //            //    Console.WriteLine("New ambiguity group of " + pep + ": " + string.Join(",", AmbiguityGroups[pep]));
            //        }
            //    }
            //    //Console.WriteLine("AmbiguityGroups:");
            //    //foreach (var a in AmbiguityGroups.Where(b => b.Value.Count > 0))
            //    //    Console.WriteLine(" AmbiguityGroup of " + a.Key + ": " + string.Join(",", string.Join(",", a.Value)));

            //}

            //foreach (var a in AmbiguityGroups.Where(b => b.Value != null && b.Value.Count > 0).ToList())
            //{
            //    AmbiguityGroups[a.Key] = new HashSet<int>(a.Value.Where(b => Math.Abs(myDictionary[b].mass - myDictionary[a.Key].mass) < 0.001 && !vetted[b]));
            //}

            //foreach (var a in AmbiguityGroups.Where(b => b.Value != null && b.Value.Count > 0).OrderByDescending(b => b.Value.Count))
            //{
            //    Console.WriteLine(a.Key + " of mass " + myDictionary[a.Key].mass.ToString() + " ambiguous with: " + string.Join(", ", a.Value.Select(b => b.ToString() + " of mass " + myDictionary[b].mass)));

            //    var compactPeptide = myDictionary[a.Key];
            //    PeptideWithSetModifications pepSet = GetPeptideWithSetModifications(compactPeptide, hm, fixedModifications, localizeableModifications, variableModifications);
            //    Console.WriteLine(pepSet.ExtendedSequence);
            //    foreach (var huhu in a.Value)
            //    {
            //        var compactPeptideInGroup = myDictionary[huhu];
            //        PeptideWithSetModifications pepSetInGroup = GetPeptideWithSetModifications(compactPeptideInGroup, hm, fixedModifications, localizeableModifications, variableModifications);
            //        Console.WriteLine(pepSetInGroup.ExtendedSequence);

            //    }
            //}

            Console.WriteLine("finished generating peptide index");
            myDictionaryOut = myDictionary;
            myFragmentDictionaryOut = myFragmentDictionary;

            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = string.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                ts.Hours, ts.Minutes, ts.Seconds,
                ts.Milliseconds / 10);
            Console.WriteLine("Time to generate indices: " + elapsedTime);
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