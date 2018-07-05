using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Neo
{
    public static class NeoExport
    {
        public static string Path;
        public static string Folder;
        public static CommonParameters CommonParameters;

        public static void ExportAll(List<NeoPsm> psms, Ms2ScanWithSpecificMass[] spectra, string databaseFileName)
        {
            Folder = DateTime.Now.ToString("yyyy-MM-dd_hh-mm-ss");
            Path = "";
            string[] temp = databaseFileName.Split('\\').ToArray();
            for (int i = 0; i < temp.Length - 1; i++)
                Path += temp[i] + '\\';

            Directory.CreateDirectory(Path + Folder);
            ExportCandidates(psms, spectra, Path, CommonParameters);
            // ExportFullFASTA(psms, databaseFileName, path);
            ExportFASTAAppendix(psms, databaseFileName, Path);

            //  ExportFilteredFusionPeptideAppendix(psms, databaseFileName, path);
        }

        public static void ExportCandidates(List<NeoPsm> psms, Ms2ScanWithSpecificMass[] spectra, string path, CommonParameters commonParameters)
        {
            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "ExportedFusionCandidatesAll.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                //double progress = 0;
                Parallel.ForEach(psms, new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (psm) =>
                {
                    Ms2ScanWithSpecificMass spectrum = spectra[psm.ScanNumber];
                    //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                    string allPossibleSequences = "";
                    string mostProbableSequences = "";
                    int indexOfFirstProbableSequence = -1;
                    string mostProbableParents = "";
                    string allPossibleParents = "";
                    FusionType probableType = psm.FusionType;
                    for (int fc = 0; fc < psm.Candidates.Count; fc++)
                    {
                        FusionCandidate fusionCandidate = psm.Candidates[fc]; //need fc for indexOfFirstProbableSequence
                        char[] tempArray = fusionCandidate.Seq.ToCharArray();
                        //if most probable, add to all and most
                        if (fusionCandidate.FusionType.Equals(probableType))
                        {
                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Length; i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.JunctionIndexes)
                                {
                                    if (junction == i)
                                    {
                                        mostProbableSequences += "-";
                                        allPossibleSequences += "-";
                                    }
                                }
                            }
                            mostProbableSequences += "|";
                            allPossibleSequences += "|";

                            //record parents
                            string tempParents = "";
                            //switch (probableType)
                            //{
                            //    case FusionCandidate.FusionType.TL:
                            //        if (fusionCandidate.translatedParents.Count == 0)
                            //            throw new Exception();
                            //        tempParents += GenerateParentOutput(fusionCandidate.translatedParents, new List<CisParent>(), new List<TransParent>());
                            //        break;

                            //    case FusionCandidate.FusionType.NC:
                            //    case FusionCandidate.FusionType.RC:
                            //        if (fusionCandidate.cisParents.Count == 0)
                            //            throw new Exception();
                            //        tempParents += GenerateParentOutput(new List<TranslatedParent>(), fusionCandidate.cisParents, new List<TransParent>());
                            //        break;

                            //    default: //if trans
                            //        if (fusionCandidate.transParents.Count == 0)
                            //            throw new Exception();
                            //        tempParents += GenerateParentOutput(new List<TranslatedParent>(), new List<CisParent>(), fusionCandidate.transParents);
                            //        break;
                            //}
                            mostProbableParents += tempParents;
                            allPossibleParents += tempParents;
                        }
                        else //not most probable, only add it to allPossibleSequences
                        {
                            //record sequences
                            for (int i = 0; i < tempArray.Length; i++)
                            {
                                allPossibleSequences += tempArray[i];
                                if (fusionCandidate.JunctionIndexes.Contains(i))
                                {
                                    allPossibleSequences += "-";
                                }
                            }
                            allPossibleSequences += "|";
                            //record parents
                            //UNCOMMENT   allPossibleParents += GenerateParentOutput(fusionCandidate.translatedParents, fusionCandidate.cisParents, fusionCandidate.transParents);
                        }
                        /*        foreach(ParentInfo PI in fusionCandidate.parentInfo)
                                {
                                    parents += PI.accession + "_" + PI.parentType.ToString() + "_" + PI.seqFound + "|";
                                }*/
                    }

                    allPossibleSequences = allPossibleSequences.Substring(0, allPossibleSequences.Length - 1); //remove last "|"
                    mostProbableSequences = mostProbableSequences.Substring(0, mostProbableSequences.Length - 1); //remove last "|"

                    string ambiguity = "";
                    NeoFindAmbiguity.FindIons(psm.Candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                            //      mutableError_message += e;
                    bool[] foundIons = psm.Candidates[indexOfFirstProbableSequence].FoundIons;
                    char[] firstSeq = psm.Candidates[indexOfFirstProbableSequence].Seq.ToCharArray();
                    //   if(foundIons.Length==firstSeq.Length) //prevent crashing if something went wrong
                    // {
                    bool ambiguous = false;
                    for (int i = 1; i < foundIons.Length; i++)
                    {
                        if (foundIons[i]) //if found
                        {
                            ambiguity += firstSeq[i - 1]; //add aa
                            if (ambiguous) //if it is part of an ambiguous sequence
                            {
                                ambiguity += ")";
                                ambiguous = false; //no longer ambiguous
                            }
                        }
                        else
                        {
                            if (!ambiguous)
                            {
                                ambiguous = true;
                                ambiguity += "(";
                            }
                            ambiguity += firstSeq[i - 1];
                        }
                    }
                    ambiguity += firstSeq[foundIons.Length - 1];
                    if (ambiguous)
                    {
                        ambiguity += ")";
                    }
                    string potentialFalsePositives = "";
                    //foreach (Variant v in psm.variants)
                    //{
                    //    potentialFalsePositives += v.id + "_" + v.start + "-" + (v.start + v.peptideLength - 1) + "(" + v.pepSeq + ")" + v.varType + "|";
                    //}
                    if (potentialFalsePositives.Length > 0)
                    {
                        potentialFalsePositives = potentialFalsePositives.Substring(0, potentialFalsePositives.Length - 1); //remove last |
                        if (potentialFalsePositives.Length > 30000)
                            potentialFalsePositives = potentialFalsePositives.Substring(0, 30000);
                    }
                    //workarounds for excel. Actual limit is 32767, but that doesn't seem to work
                    if (mostProbableParents.Length > 30000)
                    {
                        mostProbableParents = mostProbableParents.Substring(0, 30000);
                    }
                    if (allPossibleParents.Length > 30000)
                    {
                        allPossibleParents = allPossibleParents.Substring(0, 30000);
                    }
                    int score = 0;
                    foreach (bool b in psm.Candidates[indexOfFirstProbableSequence].FoundIons)
                    {
                        if (b)
                        {
                            score++;
                        }
                    }
                    lock (file)
                    {
                        file.WriteLine(psm.ScanNumber.ToString() + '\t' + psm.ExpMass.ToString() + '\t' + psm.NInfo.Seq + '\t' + psm.NInfo.Score + '\t' + psm.CInfo.Seq + '\t' + psm.CInfo.Score + '\t' + psm.Candidates[indexOfFirstProbableSequence].Seq + '\t' + ambiguity + '\t' + psm.FusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.Candidates.Count.ToString() + '\t' + potentialFalsePositives + '\t' + score);
                        //progress++;
                        //this.worker.ReportProgress((int)(progress / psms.Count() * 100));
                    }
                });
            }

            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "ExportedFusionCandidatesTL.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                Parallel.ForEach(psms, new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (psm) =>
                {
                    if (psm.Candidates.Any(x => x.FusionType == FusionType.TL))
                    {
                        Ms2ScanWithSpecificMass spectrum = spectra[psm.ScanNumber];
                        //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                        string allPossibleSequences = "";
                        string mostProbableSequences = "";
                        int indexOfFirstProbableSequence = -1;
                        string mostProbableParents = "";
                        string allPossibleParents = "";
                        FusionType probableType = psm.FusionType;
                        for (int fc = 0; fc < psm.Candidates.Count; fc++)
                        {
                            FusionCandidate fusionCandidate = psm.Candidates[fc]; //need fc for indexOfFirstProbableSequence
                            if (fusionCandidate.FusionType != FusionType.TL)
                                continue;
                            char[] tempArray = fusionCandidate.Seq.ToCharArray();
                            //if most probable, add to all and most

                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Length; i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.JunctionIndexes)
                                {
                                    if (junction == i)
                                    {
                                        mostProbableSequences += "-";
                                        allPossibleSequences += "-";
                                    }
                                }
                            }
                            mostProbableSequences += "|";
                            allPossibleSequences += "|";

                            //record parents
                            string tempParents = "";

                            mostProbableParents += tempParents;
                            allPossibleParents += tempParents;
                        }

                        allPossibleSequences = allPossibleSequences.Substring(0, allPossibleSequences.Length - 1); //remove last "|"
                        mostProbableSequences = mostProbableSequences.Substring(0, mostProbableSequences.Length - 1); //remove last "|"

                        string ambiguity = "";
                        NeoFindAmbiguity.FindIons(psm.Candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                                //      mutableError_message += e;
                        bool[] foundIons = psm.Candidates[indexOfFirstProbableSequence].FoundIons;
                        char[] firstSeq = psm.Candidates[indexOfFirstProbableSequence].Seq.ToCharArray();

                        bool ambiguous = false;
                        for (int i = 1; i < foundIons.Length; i++)
                        {
                            if (foundIons[i]) //if found
                            {
                                ambiguity += firstSeq[i - 1]; //add aa
                                if (ambiguous) //if it is part of an ambiguous sequence
                                {
                                    ambiguity += ")";
                                    ambiguous = false; //no longer ambiguous
                                }
                            }
                            else
                            {
                                if (!ambiguous)
                                {
                                    ambiguous = true;
                                    ambiguity += "(";
                                }
                                ambiguity += firstSeq[i - 1];
                            }
                        }
                        ambiguity += firstSeq[foundIons.Length - 1];
                        if (ambiguous)
                        {
                            ambiguity += ")";
                        }
                        string potentialFalsePositives = "";

                        if (potentialFalsePositives.Length > 0)
                        {
                            potentialFalsePositives = potentialFalsePositives.Substring(0, potentialFalsePositives.Length - 1); //remove last |
                            if (potentialFalsePositives.Length > 30000)
                            {
                                potentialFalsePositives = potentialFalsePositives.Substring(0, 30000);
                            }
                        }
                        //workarounds for excel. Actual limit is 32767, but that doesn't seem to work
                        if (mostProbableParents.Length > 30000)
                        {
                            mostProbableParents = mostProbableParents.Substring(0, 30000);
                        }
                        if (allPossibleParents.Length > 30000)
                        {
                            allPossibleParents = allPossibleParents.Substring(0, 30000);
                        }
                        int score = 0;
                        foreach (bool b in psm.Candidates[indexOfFirstProbableSequence].FoundIons)
                            if (b)
                                score++;
                        lock (file)
                        {
                            file.WriteLine(psm.ScanNumber.ToString() + '\t' + psm.ExpMass.ToString() + '\t' + psm.NInfo.Seq + '\t' + psm.NInfo.Score + '\t' + psm.CInfo.Seq + '\t' + psm.CInfo.Score + '\t' + psm.Candidates[indexOfFirstProbableSequence].Seq + '\t' + ambiguity + '\t' + psm.FusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.Candidates.Count.ToString() + '\t' + potentialFalsePositives + '\t' + score);
                            //progress++;
                            //this.worker.ReportProgress((int)(progress / psms.Count() * 100));
                        }
                    }
                });
            }

            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "ExportedFusionCandidatesNC.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                Parallel.ForEach(psms, new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (psm) =>
                {
                    if (psm.Candidates.Any(x => x.FusionType == FusionType.NC))
                    {
                        Ms2ScanWithSpecificMass spectrum = spectra[psm.ScanNumber];
                        //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                        string allPossibleSequences = "";
                        string mostProbableSequences = "";
                        int indexOfFirstProbableSequence = -1;
                        string mostProbableParents = "";
                        string allPossibleParents = "";
                        FusionType probableType = psm.FusionType;
                        for (int fc = 0; fc < psm.Candidates.Count; fc++)
                        {
                            FusionCandidate fusionCandidate = psm.Candidates[fc]; //need fc for indexOfFirstProbableSequence
                            if (fusionCandidate.FusionType != FusionType.NC)
                            {
                                continue;
                            }
                            char[] tempArray = fusionCandidate.Seq.ToCharArray();
                            //if most probable, add to all and most

                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Length; i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.JunctionIndexes)
                                {
                                    if (junction == i)
                                    {
                                        mostProbableSequences += "-";
                                        allPossibleSequences += "-";
                                    }
                                }
                            }
                            mostProbableSequences += "|";
                            allPossibleSequences += "|";

                            //record parents
                            string tempParents = "";

                            mostProbableParents += tempParents;
                            allPossibleParents += tempParents;
                        }

                        allPossibleSequences = allPossibleSequences.Substring(0, allPossibleSequences.Length - 1); //remove last "|"
                        mostProbableSequences = mostProbableSequences.Substring(0, mostProbableSequences.Length - 1); //remove last "|"

                        string ambiguity = "";
                        NeoFindAmbiguity.FindIons(psm.Candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                                //      mutableError_message += e;
                        bool[] foundIons = psm.Candidates[indexOfFirstProbableSequence].FoundIons;
                        char[] firstSeq = psm.Candidates[indexOfFirstProbableSequence].Seq.ToCharArray();

                        bool ambiguous = false;
                        for (int i = 1; i < foundIons.Length; i++)
                        {
                            if (foundIons[i]) //if found
                            {
                                ambiguity += firstSeq[i - 1]; //add aa
                                if (ambiguous) //if it is part of an ambiguous sequence
                                {
                                    ambiguity += ")";
                                    ambiguous = false; //no longer ambiguous
                                }
                            }
                            else
                            {
                                if (!ambiguous)
                                {
                                    ambiguous = true;
                                    ambiguity += "(";
                                }
                                ambiguity += firstSeq[i - 1];
                            }
                        }
                        ambiguity += firstSeq[foundIons.Length - 1];
                        if (ambiguous)
                            ambiguity += ")";
                        string potentialFalsePositives = "";

                        if (potentialFalsePositives.Length > 0)
                        {
                            potentialFalsePositives = potentialFalsePositives.Substring(0, potentialFalsePositives.Length - 1); //remove last |
                            if (potentialFalsePositives.Length > 30000)
                            {
                                potentialFalsePositives = potentialFalsePositives.Substring(0, 30000);
                            }
                        }
                        //workarounds for excel. Actual limit is 32767, but that doesn't seem to work
                        if (mostProbableParents.Length > 30000)
                        {
                            mostProbableParents = mostProbableParents.Substring(0, 30000);
                        }
                        if (allPossibleParents.Length > 30000)
                        {
                            allPossibleParents = allPossibleParents.Substring(0, 30000);
                        }
                        int score = 0;
                        foreach (bool b in psm.Candidates[indexOfFirstProbableSequence].FoundIons)
                        {
                            if (b)
                            {
                                score++;
                            }
                        }
                        lock (file)
                        {
                            file.WriteLine(psm.ScanNumber.ToString() + '\t' + psm.ExpMass.ToString() + '\t' + psm.NInfo.Seq + '\t' + psm.NInfo.Score + '\t' + psm.CInfo.Seq + '\t' + psm.CInfo.Score + '\t' + psm.Candidates[indexOfFirstProbableSequence].Seq + '\t' + ambiguity + '\t' + psm.FusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.Candidates.Count.ToString() + '\t' + potentialFalsePositives + '\t' + score);
                            //progress++;
                            //this.worker.ReportProgress((int)(progress / psms.Count() * 100));
                        }
                    }
                });
            }

            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "ExportedFusionCandidatesTS.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                Parallel.ForEach(psms, new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (psm) =>
                {
                    if (psm.Candidates.Any(x => x.FusionType == FusionType.TS))
                    {
                        Ms2ScanWithSpecificMass spectrum = spectra[psm.ScanNumber];
                        //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                        string allPossibleSequences = "";
                        string mostProbableSequences = "";
                        int indexOfFirstProbableSequence = -1;
                        string mostProbableParents = "";
                        string allPossibleParents = "";
                        FusionType probableType = psm.FusionType;
                        for (int fc = 0; fc < psm.Candidates.Count; fc++)
                        {
                            FusionCandidate fusionCandidate = psm.Candidates[fc]; //need fc for indexOfFirstProbableSequence
                            if (fusionCandidate.FusionType != FusionType.TS)
                                continue;
                            char[] tempArray = fusionCandidate.Seq.ToCharArray();
                            //if most probable, add to all and most

                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Length; i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.JunctionIndexes)
                                {
                                    if (junction == i)
                                    {
                                        mostProbableSequences += "-";
                                        allPossibleSequences += "-";
                                    }
                                }
                            }
                            mostProbableSequences += "|";
                            allPossibleSequences += "|";

                            //record parents
                            string tempParents = "";

                            mostProbableParents += tempParents;
                            allPossibleParents += tempParents;
                        }

                        allPossibleSequences = allPossibleSequences.Substring(0, allPossibleSequences.Length - 1); //remove last "|"
                        mostProbableSequences = mostProbableSequences.Substring(0, mostProbableSequences.Length - 1); //remove last "|"

                        string ambiguity = "";
                        NeoFindAmbiguity.FindIons(psm.Candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                                //      mutableError_message += e;
                        bool[] foundIons = psm.Candidates[indexOfFirstProbableSequence].FoundIons;
                        char[] firstSeq = psm.Candidates[indexOfFirstProbableSequence].Seq.ToCharArray();

                        bool ambiguous = false;
                        for (int i = 1; i < foundIons.Length; i++)
                        {
                            if (foundIons[i]) //if found
                            {
                                ambiguity += firstSeq[i - 1]; //add aa
                                if (ambiguous) //if it is part of an ambiguous sequence
                                {
                                    ambiguity += ")";
                                    ambiguous = false; //no longer ambiguous
                                }
                            }
                            else
                            {
                                if (!ambiguous)
                                {
                                    ambiguous = true;
                                    ambiguity += "(";
                                }
                                ambiguity += firstSeq[i - 1];
                            }
                        }
                        ambiguity += firstSeq[foundIons.Length - 1];
                        if (ambiguous)
                        {
                            ambiguity += ")";
                        }
                        string potentialFalsePositives = "";

                        if (potentialFalsePositives.Length > 0)
                        {
                            potentialFalsePositives = potentialFalsePositives.Substring(0, potentialFalsePositives.Length - 1); //remove last |
                            if (potentialFalsePositives.Length > 30000)
                                potentialFalsePositives = potentialFalsePositives.Substring(0, 30000);
                        }
                        //workarounds for excel. Actual limit is 32767, but that doesn't seem to work
                        if (mostProbableParents.Length > 30000)
                        {
                            mostProbableParents = mostProbableParents.Substring(0, 30000);
                        }
                        if (allPossibleParents.Length > 30000)
                        {
                            allPossibleParents = allPossibleParents.Substring(0, 30000);
                        }
                        int score = 0;
                        foreach (bool b in psm.Candidates[indexOfFirstProbableSequence].FoundIons)
                        {
                            if (b)
                            {
                                score++;
                            }
                        }
                        lock (file)
                        {
                            file.WriteLine(psm.ScanNumber.ToString() + '\t' + psm.ExpMass.ToString() + '\t' + psm.NInfo.Seq + '\t' + psm.NInfo.Score + '\t' + psm.CInfo.Seq + '\t' + psm.CInfo.Score + '\t' + psm.Candidates[indexOfFirstProbableSequence].Seq + '\t' + ambiguity + '\t' + psm.FusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.Candidates.Count.ToString() + '\t' + potentialFalsePositives + '\t' + score);
                            //progress++;
                            //this.worker.ReportProgress((int)(progress / psms.Count() * 100));
                        }
                    }
                });
            }
        }

        private static void ExportFullFASTA(List<NeoPsm> psms, string databaseFileName, string path)
        {
            //@"C:\Users\Zach Rolfs\Desktop\Chemistry\Smith Research\Fusion Peptides\Neo\Results\"
            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "FullFusionDatabase.fasta"))
            {
                //copy database
                string[] FASTARead = File.ReadAllLines(databaseFileName);
                foreach (string s in FASTARead)
                {
                    file.WriteLine(s);
                }

                //Write Fusion Peptide Candidates
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for FASTA output
                    string scan = psm.ScanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.Candidates)
                    {
                        file.WriteLine(">sp|" + scan + fusionCandidate.FusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.FusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.Seq;
                        for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file. It is unlikely to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Length)
                            {
                                file.WriteLine(seq.Substring(i, 60));
                            }
                            else
                            {
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                        }
                        fusNum++;
                    }
                }
            }
        }

        private static void ExportFASTAAppendix(List<NeoPsm> psms, string databaseFileName, string path)
        {
            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "FusionDatabaseAppendixAll.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.ScanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.Candidates)
                    {
                        file.WriteLine(">sp|" + scan + fusionCandidate.FusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.FusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.Seq;
                        for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Length)
                            {
                                file.WriteLine(seq.Substring(i, 60));
                            }
                            else
                            {
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "FusionDatabaseAppendixTL.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for FASTA output
                    string scan = psm.ScanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.Candidates)
                    {
                        if (fusionCandidate.FusionType != FusionType.TL)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.FusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.FusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.Seq;
                        for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Length)
                            {
                                file.WriteLine(seq.Substring(i, 60));
                            }
                            else
                            {
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "FusionDatabaseAppendixNC.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.ScanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.Candidates)
                    {
                        if (fusionCandidate.FusionType != FusionType.NC)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.FusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.FusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.Seq;
                        for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Length)
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "FusionDatabaseAppendixTS.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.ScanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.Candidates)
                    {
                        if (fusionCandidate.FusionType != FusionType.TS)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.FusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.FusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.Seq;
                        for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Length)
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "FusionDatabaseAppendixTop.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.ScanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.Candidates)
                    {
                        if (fusionCandidate.FusionType != psm.FusionType)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.FusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.FusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.Seq;
                        for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Length)
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }
        }

        private static void printPTMInfo(string ptmName, int position, StreamWriter file)
        {
            file.WriteLine("<feature type=" + '"' + "modified residue" + '"' + " description=" + '"' + ptmName + '"' + ">");
            file.WriteLine("<location>");
            file.WriteLine("<position position=" + '"' + position + '"' + " />");
            file.WriteLine("</location>");
            file.WriteLine("</feature>");
        }

        private static void ExportFilteredFusionPeptideAppendix(List<NeoPsm> psms, string databaseFileName, string path)
        {
            using (StreamWriter file = new StreamWriter(path + Folder + @"\" + Folder + "FilteredFusionDatabaseAppendix.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;
                    if (!psm.FusionType.Equals(FusionType.TL))
                    {
                        //the following code is for FASTA output
                        string scan = psm.ScanNumber.ToString();
                        foreach (FusionCandidate fusionCandidate in psm.Candidates)
                        {
                            file.WriteLine(">sp|" + scan + fusionCandidate.FusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.FusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                            string seq = fusionCandidate.Seq;
                            for (int i = 0; i < seq.Length; i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                            {
                                if ((i + 60) < seq.Length)
                                    file.WriteLine(seq.Substring(i, 60));
                                else
                                    file.WriteLine(seq.Substring(i, seq.Length - i));
                            }
                            fusNum++;
                        }
                    }
                }
            }
        }

        private static string GenerateParentOutput(List<TranslatedParent> translatedParents, List<CisParent> cisParents, List<TransParent> transParents)
        {
            string output = "";
            foreach (TranslatedParent tlp in translatedParents)
            {
                output += tlp.ID + "_" + tlp.Start + "-" + (tlp.Start + tlp.PeptideLength - 1) + "(" + tlp.Seq.Substring(tlp.Start, tlp.PeptideLength) + ")" + "|";
            }
            foreach (CisParent cp in cisParents)
            {
                foreach (int ns in cp.NStart)
                {
                    foreach (int cs in cp.CStart)
                    {
                        output += cp.ID + "_" + ns + "-" + (ns + cp.NLength - 1) + "(" + cp.Seq.Substring(ns, cp.NLength) + ")"
                            + "&" + cs + "-" + (cs + cp.CLength - 1) + "(" + cp.Seq.Substring(cs, cp.CLength) + ")" + "|";
                    }
                }
            }
            string nOutput = "";
            string cOutput = "";
            foreach (TransParent tp in transParents)
            {
                if (tp.Terminal.Equals(ParentInfo.Terminal.N))
                    foreach (int ts in tp.Start)
                        nOutput += tp.ID + "_" + ts + "-" + (ts + tp.PeptideLength - 1) + "(" + tp.Seq.Substring(ts, tp.PeptideLength) + ")" + "|";
                else
                    foreach (int ts in tp.Start)
                        cOutput += tp.ID + "_" + ts + "-" + (ts + tp.PeptideLength - 1) + "(" + tp.Seq.Substring(ts, tp.PeptideLength) + ")" + "|";
            }
            if (nOutput.Length > 0) //if there was a trans, append and delimit the two lists
                output += nOutput + "&&|" + cOutput;
            output = output.Substring(0, output.Length - 1); //remove last |

            return output;
        }
    }
}