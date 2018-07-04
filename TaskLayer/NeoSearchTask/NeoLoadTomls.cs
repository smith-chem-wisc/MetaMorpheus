using EngineLayer;
using Nett;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;

namespace TaskLayer
{
    public static class NeoLoadTomls
    {
        public static List<MetaMorpheusTask> LoadTomls(NeoSearchTask ye5) //tomls are located in EngineLayer//Neo//Data//TomlFiles
        {
            List<MetaMorpheusTask> novelCollection = new List<MetaMorpheusTask>();

            string defaultFolderPath = Path.Combine(GlobalVariables.DataDir, @"Neo", @"TomlFiles");

            // write TOML
            var tomlFileName = Path.Combine(defaultFolderPath, ye5.GetType().Name + "config.toml");
            Toml.WriteFile(ye5, tomlFileName, MetaMorpheusTask.tomlConfig);

            if (ye5.NeoParameters.Calibrate)
            {
                string caliFileName = "CalibrationTaskconfig.toml";
                caliFileName = Path.Combine(defaultFolderPath, caliFileName);
                UpdateTomls(tomlFileName, caliFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<CalibrationTask>(caliFileName, MetaMorpheusTask.tomlConfig);//FIXME, An item with the same key has already been added, dictionary in Toml.ReadFile
                novelCollection.Add(yeo);//multiple protease issue
            }

            if (ye5.NeoParameters.GPTMD)
            {
                string gptmdFileName = "GptmdTaskconfig.toml";
                gptmdFileName = Path.Combine(defaultFolderPath, gptmdFileName);
                UpdateTomls(tomlFileName, gptmdFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<GptmdTask>(gptmdFileName, MetaMorpheusTask.tomlConfig);
                novelCollection.Add(yeo);
            }

            if (ye5.NeoParameters.TargetSearch)
            {
                string targetFileName = "SearchTaskTargetconfig.toml";
                targetFileName = Path.Combine(defaultFolderPath, targetFileName);
                UpdateTomls(tomlFileName, targetFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                novelCollection.Add(yeo);
            }

            if (ye5.NeoParameters.DecoySearch)
            {
                string decoyFileName = "SearchTaskDecoyconfig.toml";
                decoyFileName = Path.Combine(defaultFolderPath, decoyFileName);
                UpdateTomls(tomlFileName, decoyFileName, ye5.CommonParameters, TerminusType.None, false);
                var yeo = Toml.ReadFile<SearchTask>(decoyFileName, MetaMorpheusTask.tomlConfig);
                novelCollection.Add(yeo);
            }

            var yeo5_1 = ye5.Clone();
            yeo5_1.NeoType = NeoSearchTask.NeoTaskType.AggregateTargetDecoyFiles;

            novelCollection.Add(yeo5_1);

            if (ye5.NeoParameters.SearchNTerminus)
            {
                string targetFileName = "SearchTaskNconfig.toml";
                targetFileName = Path.Combine(defaultFolderPath, targetFileName);
                UpdateTomls(tomlFileName, targetFileName, ye5.CommonParameters, TerminusType.N, false);
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                novelCollection.Add(yeo);
            }

            if (ye5.NeoParameters.SearchCTerminus)
            {
                string targetFileName = "SearchTaskCconfig.toml";
                targetFileName = Path.Combine(defaultFolderPath, targetFileName);
                UpdateTomls(tomlFileName, targetFileName, ye5.CommonParameters, TerminusType.C, false);
                var yeo = Toml.ReadFile<SearchTask>(targetFileName, MetaMorpheusTask.tomlConfig);
                novelCollection.Add(yeo);
            }

            var yeo5_2 = ye5.Clone();
            yeo5_2.NeoType = NeoSearchTask.NeoTaskType.GenerateSplicedPeptides;
            novelCollection.Add(yeo5_2);

            string cisFileName = "SearchTaskCisconfig.toml";
            cisFileName = Path.Combine(defaultFolderPath, cisFileName);
            UpdateTomls(tomlFileName, cisFileName, ye5.CommonParameters, TerminusType.None, true);
            var yeocis = Toml.ReadFile<SearchTask>(cisFileName, MetaMorpheusTask.tomlConfig);
            novelCollection.Add(yeocis);

            var yeo5_3 = ye5.Clone();
            yeo5_3.NeoType = NeoSearchTask.NeoTaskType.SearchTransDb;
            novelCollection.Add(yeo5_3);

            string transFileName = "SearchTaskTransconfig.toml";
            transFileName = Path.Combine(defaultFolderPath, transFileName);
            UpdateTomls(tomlFileName, transFileName, ye5.CommonParameters, TerminusType.None, true);
            var yeotrans = Toml.ReadFile<SearchTask>(transFileName, MetaMorpheusTask.tomlConfig);
            novelCollection.Add(yeotrans);

            var yeo5_4 = ye5.Clone();
            yeo5_4.NeoType = NeoSearchTask.NeoTaskType.AggregateNormalSplicedFiles;
            novelCollection.Add(yeo5_4);

            // DeleteTomlFile

            File.Delete(tomlFileName);

            return novelCollection;
        }

        private static void UpdateTomls(string tomlFileName, string fileName, CommonParameters ye5, TerminusType terminusType, bool spliceSearch)
        {
            string[] oldTomlLines = File.ReadAllLines(@fileName);
            List<string> newTomlLines = new List<string>();
            foreach (string line in oldTomlLines)
            {
                if (line.Contains("LocalizeAll") && terminusType.Equals(TerminusType.None))
                    newTomlLines.Add(GetCorrectValue("LocalizeAll", tomlFileName, line));
                else if (line.Contains("ListOfModsFixed"))
                    newTomlLines.Add(GetCorrectValue("ListOfModsFixed", tomlFileName, line));
                else if (line.Contains("ListOfModsVariable") && terminusType.Equals(TerminusType.None) && !spliceSearch)
                    newTomlLines.Add(GetCorrectValue("ListOfModsVariable", tomlFileName, line));
                else if (line.Contains("BIons"))
                {
                    if (terminusType.Equals(TerminusType.N) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("BIons", tomlFileName, line));
                    else
                        newTomlLines.Add("BIons = false");
                }
                else if (line.Contains("YIons"))
                {
                    if (terminusType.Equals(TerminusType.C) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("YIons", tomlFileName, line));
                    else
                        newTomlLines.Add("YIons = false");
                }
                else if (line.Contains("ZdotIons"))
                {
                    if (terminusType.Equals(TerminusType.C) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("ZdotIons", tomlFileName, line));
                    else
                        newTomlLines.Add("ZdotIons = false");
                }
                else if (line.Contains("CIons"))
                {
                    if (terminusType.Equals(TerminusType.N) || terminusType.Equals(TerminusType.None))
                        newTomlLines.Add(GetCorrectValue("CIons", tomlFileName, line));
                    else
                        newTomlLines.Add("CIons = false");
                }
                else if (line.Contains("ProductMassTolerance"))
                    newTomlLines.Add(GetCorrectValue("ProductMassTolerance", tomlFileName, line));
                else if (line.Contains("PrecursorMassTolerance"))
                    newTomlLines.Add(GetCorrectValue("PrecursorMassTolerance", tomlFileName, line));
                else if (line.Contains("MaxMissedCleavages"))
                    newTomlLines.Add(GetCorrectValue("MaxMissedCleavages", tomlFileName, line));
                else if (line.Contains("InitiatorMethionineBehavior"))
                    newTomlLines.Add(GetCorrectValue("InitiatorMethionineBehavior", tomlFileName, line));
                else if (line.Contains("MinPeptideLength") && !!terminusType.Equals(TerminusType.None))
                    newTomlLines.Add(GetCorrectValue("MinPeptideLength", tomlFileName, line));
                else if (line.Contains("MaxPeptideLength"))
                    newTomlLines.Add(GetCorrectValue("MaxPeptideLength", tomlFileName, line));
                else if (line.Contains("MaxModificationIsoforms"))
                    newTomlLines.Add(GetCorrectValue("MaxModificationIsoforms", tomlFileName, line));
                else if (line.Contains("MaxModsForPeptide"))
                    newTomlLines.Add(GetCorrectValue("MaxModsForPeptide", tomlFileName, line));
                else if (line.Contains("SemiProteaseDigestion"))
                    newTomlLines.Add(GetCorrectValue("SemiProteaseDigestion", tomlFileName, line));
                else if (line.Contains("TerminusTypeSemiProtease"))
                    newTomlLines.Add(GetCorrectValue("TerminusTypeSemiProtease", tomlFileName, line));
                else if (line.Contains("Protease") && terminusType.Equals(TerminusType.None) && !spliceSearch) //this must be last, else other names including protease will be overwritten and crash.
                    newTomlLines.Add(GetCorrectValue("Protease", tomlFileName, line));
                else
                    newTomlLines.Add(line);
            }
            using (StreamWriter file = new StreamWriter(fileName))
            {
                foreach (string line in newTomlLines)
                {
                    file.WriteLine(line);
                }
            }
        }

        private static string GetCorrectValue(string parameter, string tomlFileName, string oldLine)
        {
            string[] newTomlLines = File.ReadAllLines(@tomlFileName);
            foreach (string line in newTomlLines)
            {
                if (line.Contains(parameter))
                {
                    return line;
                }
            }
            return oldLine;
        }
    }
}