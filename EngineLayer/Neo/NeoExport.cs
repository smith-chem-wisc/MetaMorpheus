using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Neo
{
    public static class NeoExport
    {
        public static string path;
        public static string folder;

        public static void ExportAll(List<NeoPsm> psms, Ms2ScanWithSpecificMass[] spectra, string databaseFileName)
        {
            folder = DateTime.Now.ToString("yyyy-MM-dd_hh-mm-ss");
            path = "";
            string[] temp = databaseFileName.Split('\\').ToArray();
            for (int i = 0; i < temp.Count() - 1; i++)
                path += temp[i] + '\\';
            
            Directory.CreateDirectory(path + folder);
            ExportCandidates(psms, spectra, path);
           // ExportFullFASTA(psms, databaseFileName, path);
            ExportFASTAAppendix(psms, databaseFileName, path);

          //  ExportFilteredFusionPeptideAppendix(psms, databaseFileName, path);
        }

        public static void ExportCandidates(List<NeoPsm> psms, Ms2ScanWithSpecificMass[] spectra, string path)
        {
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "ExportedFusionCandidatesAll.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                //double progress = 0;
                Parallel.ForEach(psms, (psm) =>
                {
                    Ms2ScanWithSpecificMass spectrum = spectra[psm.scanNumber];
                    //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                    string allPossibleSequences = "";
                    string mostProbableSequences = "";
                    int indexOfFirstProbableSequence = -1;
                    string mostProbableParents = "";
                    string allPossibleParents = "";
                    FusionCandidate.FusionType probableType = psm.fusionType;
                    for (int fc = 0; fc < psm.candidates.Count(); fc++)
                    {
                        FusionCandidate fusionCandidate = psm.candidates[fc]; //need fc for indexOfFirstProbableSequence
                        char[] tempArray = fusionCandidate.seq.ToCharArray();
                        //if most probable, add to all and most
                        if (fusionCandidate.fusionType.Equals(probableType))
                        {
                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Count(); i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.junctionIndexes)
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
                            for (int i = 0; i < tempArray.Count(); i++)
                            {
                                allPossibleSequences += tempArray[i];
                                if (fusionCandidate.junctionIndexes.Contains(i))
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
                    NeoFindAmbiguity.FindIons(psm.candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                                              //      mutableError_message += e;
                    bool[] foundIons = psm.candidates[indexOfFirstProbableSequence].foundIons;
                    char[] firstSeq = psm.candidates[indexOfFirstProbableSequence].seq.ToCharArray();
                    //   if(foundIons.Count()==firstSeq.Count()) //prevent crashing if something went wrong
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
                        ambiguity += ")";
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
                    foreach (bool b in psm.candidates[indexOfFirstProbableSequence].foundIons)
                        if (b)
                            score++;
                    lock (file)
                    {
                        file.WriteLine(psm.scanNumber.ToString() + '\t' + psm.expMass.ToString() + '\t' + psm.nInfo.seq + '\t' + psm.nInfo.score + '\t' + psm.cInfo.seq + '\t' + psm.cInfo.score + '\t' + psm.candidates[indexOfFirstProbableSequence].seq + '\t' + ambiguity + '\t' + psm.fusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.candidates.Count().ToString() + '\t' + potentialFalsePositives + '\t' + score);
                        //progress++;
                        //this.worker.ReportProgress((int)(progress / psms.Count() * 100));
                    }
                });
            }

            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "ExportedFusionCandidatesTL.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                double progress = 0;
                Parallel.ForEach(psms, (psm) =>
                {
                    if (psm.candidates.Any(x => x.fusionType == FusionCandidate.FusionType.TL))
                    {
                        Ms2ScanWithSpecificMass spectrum = spectra[psm.scanNumber];
                        //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                        string allPossibleSequences = "";
                        string mostProbableSequences = "";
                        int indexOfFirstProbableSequence = -1;
                        string mostProbableParents = "";
                        string allPossibleParents = "";
                        FusionCandidate.FusionType probableType = psm.fusionType;
                        for (int fc = 0; fc < psm.candidates.Count(); fc++)
                        {
                            FusionCandidate fusionCandidate = psm.candidates[fc]; //need fc for indexOfFirstProbableSequence
                            if (fusionCandidate.fusionType != FusionCandidate.FusionType.TL)
                                continue;
                            char[] tempArray = fusionCandidate.seq.ToCharArray();
                            //if most probable, add to all and most

                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Count(); i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.junctionIndexes)
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
                        NeoFindAmbiguity.FindIons(psm.candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                                                  //      mutableError_message += e;
                        bool[] foundIons = psm.candidates[indexOfFirstProbableSequence].foundIons;
                        char[] firstSeq = psm.candidates[indexOfFirstProbableSequence].seq.ToCharArray();

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
                        foreach (bool b in psm.candidates[indexOfFirstProbableSequence].foundIons)
                            if (b)
                                score++;
                        lock (file)
                        {
                            file.WriteLine(psm.scanNumber.ToString() + '\t' + psm.expMass.ToString() + '\t' + psm.nInfo.seq + '\t' + psm.nInfo.score + '\t' + psm.cInfo.seq + '\t' + psm.cInfo.score + '\t' + psm.candidates[indexOfFirstProbableSequence].seq + '\t' + ambiguity + '\t' + psm.fusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.candidates.Count().ToString() + '\t' + potentialFalsePositives + '\t' + score);
                            //progress++;
                            //this.worker.ReportProgress((int)(progress / psms.Count() * 100));
                        }
                    }
                });
            }


            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "ExportedFusionCandidatesNC.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                double progress = 0;
                Parallel.ForEach(psms, (psm) =>
                {
                    if (psm.candidates.Any(x => x.fusionType == FusionCandidate.FusionType.NC))
                    {
                        Ms2ScanWithSpecificMass spectrum = spectra[psm.scanNumber];
                        //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                        string allPossibleSequences = "";
                        string mostProbableSequences = "";
                        int indexOfFirstProbableSequence = -1;
                        string mostProbableParents = "";
                        string allPossibleParents = "";
                        FusionCandidate.FusionType probableType = psm.fusionType;
                        for (int fc = 0; fc < psm.candidates.Count(); fc++)
                        {
                            FusionCandidate fusionCandidate = psm.candidates[fc]; //need fc for indexOfFirstProbableSequence
                            if (fusionCandidate.fusionType != FusionCandidate.FusionType.NC)
                                continue;
                            char[] tempArray = fusionCandidate.seq.ToCharArray();
                            //if most probable, add to all and most

                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Count(); i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.junctionIndexes)
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
                        NeoFindAmbiguity.FindIons(psm.candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                                                  //      mutableError_message += e;
                        bool[] foundIons = psm.candidates[indexOfFirstProbableSequence].foundIons;
                        char[] firstSeq = psm.candidates[indexOfFirstProbableSequence].seq.ToCharArray();

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
                        foreach (bool b in psm.candidates[indexOfFirstProbableSequence].foundIons)
                            if (b)
                                score++;
                        lock (file)
                        {
                            file.WriteLine(psm.scanNumber.ToString() + '\t' + psm.expMass.ToString() + '\t' + psm.nInfo.seq + '\t' + psm.nInfo.score + '\t' + psm.cInfo.seq + '\t' + psm.cInfo.score + '\t' + psm.candidates[indexOfFirstProbableSequence].seq + '\t' + ambiguity + '\t' + psm.fusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.candidates.Count().ToString() + '\t' + potentialFalsePositives + '\t' + score);
                            //progress++;
                            //this.worker.ReportProgress((int)(progress / psms.Count() * 100));
                        }
                    }
                });
            }


            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "ExportedFusionCandidatesTS.txt"))
            {
                file.WriteLine("Scan" + '\t' + "ExperimentalMass" + '\t' + "OriginalNSequence" + '\t' + "OriginalNScore" + '\t' + "OriginalCSequence" + '\t' + "OriginalCScore" + '\t' + "SampleSequence" + '\t' + "Ambiguity" + '\t' + "ProbableType" + '\t' + "MostProbableSequenceJunctions" + '\t' + "MostProbableSequence(s)" + '\t' + "MostProbableParents" + '\t' + "AllPossibleSequenceJunctions" + '\t' + "AllPossibleSequence(s)" + '\t' + "AllPossibleParent(s)" + '\t' + "NumberOfPossibleSequences" + '\t' + "PotentialFalsePositives" + '\t' + "TotalScore");

                double progress = 0;
                Parallel.ForEach(psms, (psm) =>
                {
                    if (psm.candidates.Any(x => x.fusionType == FusionCandidate.FusionType.TS))
                    {
                        Ms2ScanWithSpecificMass spectrum = spectra[psm.scanNumber];
                        //printout the scan, the mass, the sequences with and without junctions, the number of potential sequences
                        string allPossibleSequences = "";
                        string mostProbableSequences = "";
                        int indexOfFirstProbableSequence = -1;
                        string mostProbableParents = "";
                        string allPossibleParents = "";
                        FusionCandidate.FusionType probableType = psm.fusionType;
                        for (int fc = 0; fc < psm.candidates.Count(); fc++)
                        {
                            FusionCandidate fusionCandidate = psm.candidates[fc]; //need fc for indexOfFirstProbableSequence
                            if (fusionCandidate.fusionType != FusionCandidate.FusionType.TS)
                                continue;
                            char[] tempArray = fusionCandidate.seq.ToCharArray();
                            //if most probable, add to all and most

                            //record sequences
                            if (indexOfFirstProbableSequence < 0)
                            {
                                indexOfFirstProbableSequence = fc;
                            }
                            for (int i = 0; i < tempArray.Count(); i++)
                            {
                                mostProbableSequences += tempArray[i];
                                allPossibleSequences += tempArray[i];
                                foreach (int junction in fusionCandidate.junctionIndexes)
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
                        NeoFindAmbiguity.FindIons(psm.candidates[indexOfFirstProbableSequence], psm, spectrum); //this should be carried over, but it's not...
                                                                                                                                  //      mutableError_message += e;
                        bool[] foundIons = psm.candidates[indexOfFirstProbableSequence].foundIons;
                        char[] firstSeq = psm.candidates[indexOfFirstProbableSequence].seq.ToCharArray();

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
                        foreach (bool b in psm.candidates[indexOfFirstProbableSequence].foundIons)
                            if (b)
                                score++;
                        lock (file)
                        {
                            file.WriteLine(psm.scanNumber.ToString() + '\t' + psm.expMass.ToString() + '\t' + psm.nInfo.seq + '\t' + psm.nInfo.score + '\t' + psm.cInfo.seq + '\t' + psm.cInfo.score + '\t' + psm.candidates[indexOfFirstProbableSequence].seq + '\t' + ambiguity + '\t' + psm.fusionType.ToString() + '\t' + mostProbableSequences + '\t' + mostProbableSequences.Replace("-", "") + '\t' + mostProbableParents + '\t' + allPossibleSequences + '\t' + allPossibleSequences.Replace("-", "") + '\t' + allPossibleParents + '\t' + psm.candidates.Count().ToString() + '\t' + potentialFalsePositives + '\t' + score);
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
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FullFusionDatabase.fasta"))
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
                    string scan = psm.scanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.candidates)
                    {
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unlikely to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }
        }

        private static void ExportFASTAAppendix(List<NeoPsm> psms, string databaseFileName, string path)
        {
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FusionDatabaseAppendixAll.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.scanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.candidates)
                    {
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FusionDatabaseAppendixTL.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //fusionPeptides = Regex.Replace(candidateRow[5].ToString(), "(\\(.*?\\))", "");//remove PTM annotations

                    //the following code is for FASTA output
                    string scan = psm.scanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.candidates)
                    {
                        if (fusionCandidate.fusionType != FusionCandidate.FusionType.TL)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FusionDatabaseAppendixNC.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.scanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.candidates)
                    {
                        if (fusionCandidate.fusionType != FusionCandidate.FusionType.NC)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FusionDatabaseAppendixTS.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.scanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.candidates)
                    {
                        if (fusionCandidate.fusionType != FusionCandidate.FusionType.TS)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
                                file.WriteLine(seq.Substring(i, 60));
                            else
                                file.WriteLine(seq.Substring(i, seq.Length - i));
                        }
                        fusNum++;
                    }
                }
            }

            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FusionDatabaseAppendixTop.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;

                    //the following code is for FASTA output
                    string scan = psm.scanNumber.ToString();
                    foreach (FusionCandidate fusionCandidate in psm.candidates)
                    {
                        if (fusionCandidate.fusionType != psm.fusionType)
                            continue;
                        file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                        string seq = fusionCandidate.seq;
                        for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                        {
                            if ((i + 60) < seq.Count())
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
            using (StreamWriter file = new StreamWriter(path + folder + @"\" + folder + "FilteredFusionDatabaseAppendix.fasta"))
            {
                foreach (NeoPsm psm in psms)
                {
                    int fusNum = 1;
                    if (!psm.fusionType.Equals(FusionCandidate.FusionType.TL))
                    {
                        //the following code is for FASTA output
                        string scan = psm.scanNumber.ToString();
                        foreach (FusionCandidate fusionCandidate in psm.candidates)
                        {
                            file.WriteLine(">sp|" + scan + fusionCandidate.fusionType.ToString() + fusNum + "| " + scan + "-" + fusNum + " Proposed " + fusionCandidate.fusionType.ToString() + " fusion peptide GN=FUS PE=1 SV=1");
                            string seq = fusionCandidate.seq;
                            for (int i = 0; i < seq.Count(); i += 60) //60 used as number of AAs per line in a FASTA file. It is unliekly to observe something this large currently, but is here as a precaution.
                            {
                                if ((i + 60) < seq.Count())
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
                output += tlp.id + "_" + tlp.start + "-" + (tlp.start + tlp.peptideLength - 1) + "(" + tlp.seq.Substring(tlp.start, tlp.peptideLength) + ")" + "|";
            }
            foreach (CisParent cp in cisParents)
            {
                foreach (int ns in cp.nStart)
                {
                    foreach (int cs in cp.cStart)
                    {
                        output += cp.id + "_" + ns + "-" + (ns + cp.nLength - 1) + "(" + cp.seq.Substring(ns, cp.nLength) + ")"
                            + "&" + cs + "-" + (cs + cp.cLength - 1) + "(" + cp.seq.Substring(cs, cp.cLength) + ")" + "|";
                    }
                }
            }
            string nOutput = "";
            string cOutput = "";
            foreach (TransParent tp in transParents)
            {
                if (tp.terminal.Equals(ParentInfo.terminal.N))
                    foreach (int ts in tp.start)
                        nOutput += tp.id + "_" + ts + "-" + (ts + tp.peptideLength - 1) + "(" + tp.seq.Substring(ts, tp.peptideLength) + ")" + "|";
                else
                    foreach (int ts in tp.start)
                        cOutput += tp.id + "_" + ts + "-" + (ts + tp.peptideLength - 1) + "(" + tp.seq.Substring(ts, tp.peptideLength) + ")" + "|";
            }
            if (nOutput.Length > 0) //if there was a trans, append and delimit the two lists
                output += nOutput + "&&|" + cOutput;
            output = output.Substring(0, output.Length - 1); //remove last |

            return output;
        }
    }
}
