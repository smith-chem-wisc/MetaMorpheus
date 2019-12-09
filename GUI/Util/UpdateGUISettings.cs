using EngineLayer;
using Nett;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;

namespace MetaMorpheusGUI
{
    public static class UpdateGUISettings
    {
        public static GuiGlobalParams Params { get; internal set; }

        public static bool LoadGUISettings()
        {
            bool fileExists = File.Exists(Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"));
            try //catches toml read errors
            {
                Params = Toml.ReadFile<GuiGlobalParams>(Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"));
            }
            catch
            {
                Params = new GuiGlobalParams(); //create an instance
                Toml.WriteFile(Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
            }

            if (GlobalVariables.MetaMorpheusVersion.Contains("Not a release version"))
            {
                Params.AskAboutUpdating = false;
            }
            return fileExists;
        }


        public static bool UseNonSpecificRecommendedSettings()
        {
            bool useRecommendedSettings = false;
            //check with the user to update params
            if (Params.AskAboutNonSpecificParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Non-Specific Recommendations?",
                    "We recommend using the following parameters for non-specific searches:\n" +
                    "\t-Check 'Non-Specific Search' (Search Task Only)\n" +
                    "\t-Use '25' for 'Max Peptide Length'.\n\n" +
                    "Would you like to use these recommended settings?");

                if (results.UseSettings)
                {
                    useRecommendedSettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    Params.AskAboutNonSpecificParams = false;
                    Params.UseNonSpecificParams = results.UseSettings;

                    Toml.WriteFile(Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                }
            }
            else if (Params.UseNonSpecificParams) //user didn't want to check in, but wanted to update last time
            {
                useRecommendedSettings = true;
            }
            //else do nothing
            return useRecommendedSettings;
        }

        public static bool UseTopDownRecommendedSettings()
        {
            bool useRecommendedSettings = false;
            //check with the user to update params
            if (Params.AskAboutTopDownParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Top-Down Recommendations?",
                    "We recommend using the following parameters for top-down searches:\n" +
                        "\t-Uncheck 'Use Provided Precursor'\n" +
                        "\t-Use '60' for 'Deconvolution Max Assumed Charge State'\n" +
                        "\t-Uncheck 'Trim MS2 Peaks'\n" +
                        "\t-Uncheck all variable mods (Please use a GPTMD database instead)\n" +
                        "\t-SEARCH TASK ONLY: Check 'No Quantification'\n" +
                        "\t-SEARCH TASK ONLY: Check '1, 2, or 3 Missed Monoisotopic Peaks'\n" +
                        "\t-GPTMD TASK ONLY: Search for only acetylation, phosphorylation, and oxidation of M\n\n" +
                    "Would you like to use these recommended settings?");

                if (results.UseSettings)
                {
                    useRecommendedSettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    Params.AskAboutTopDownParams = false;
                    Params.UseTopDownParams = results.UseSettings;

                    Toml.WriteFile(Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                }
            }
            else if (Params.UseTopDownParams) //user didn't want to check in, but wanted to update last time
            {
                useRecommendedSettings = true;
            }
            //else do nothing
            return useRecommendedSettings;
        }

        public static bool UseArgCRecommendedSettings()
        {
            bool useRecommendedSettings = false;
            //check with the user to update params
            if (Params.AskAboutArgCParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Arg-C Recommendations?",
                    "We recommend using the following parameters for Arg-C searches:\n" +
                    "\t-Use 'trypsin'\n\n" +

                    "Would you like to use these recommended settings?");

                if (results.UseSettings)
                {
                    useRecommendedSettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    Params.AskAboutArgCParams = false;
                    Params.UseArgCParams = results.UseSettings;

                    Toml.WriteFile(Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                }
            }
            else if (Params.UseArgCParams) //user didn't want to check in, but wanted to update last time
            {
                useRecommendedSettings = true;
            }
            //else do nothing
            return useRecommendedSettings;
        }

        public static bool UseChymotrypsinRecommendedSettings()
        {
            bool useRecommendedSettings = false;
            //check with the user to update params
            if (Params.AskAboutChymotrypsinParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Chymotrypsin Recommendations?",
                    "We recommend using the following parameters for chymotrypsin searches:\n" +
                    "\t-Check 'Semi-Specific Search' (Search Task Only)\n" +
                    "\t-Use '50' for the 'Max Peptide Len' (Search Task Only)\n" +
                    "\t-Use '3' for 'Max Missed Cleavages'.\n\n" +
                    "Would you like to use these recommended settings?");

                if (results.UseSettings)
                {
                    useRecommendedSettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    Params.AskAboutChymotrypsinParams = false;
                    Params.UseChymotrypsinParams = results.UseSettings;

                    Toml.WriteFile(Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                }
            }
            else if (Params.UseChymotrypsinParams) //user didn't want to check in, but wanted to update last time
            {
                useRecommendedSettings = true;
            }
            //else do nothing
            return useRecommendedSettings;
        }

        public static bool UseElastaseRecommendedSettings()
        {
            bool useRecommendedSettings = false;
            //check with the user to update params
            if (Params.AskAboutElastaseParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Elastase Recommendations?",
                    "We recommend using the following parameters for elastase searches:\n" +
                    "\t-Check 'Semi-Specific Search' (Search Task Only)\n" +
                    "\t-Use '50' for the 'Max Peptide Len' (Search Task Only)\n" +
                    "\t-Use '16' for 'Max Missed Cleavages'.\n\n" +
                    "Would you like to use these recommended settings?");

                if (results.UseSettings)
                {
                    useRecommendedSettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    Params.AskAboutElastaseParams = false;
                    Params.UseElastaseParams = results.UseSettings;

                    Toml.WriteFile(Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                }
            }
            else if (Params.UseElastaseParams) //user didn't want to check in, but wanted to update last time
            {
                useRecommendedSettings = true;
            }
            //else do nothing
            return useRecommendedSettings;
        }

        public static bool UseSemiTrypsinRecommendedSettings()
        {
            bool useRecommendedSettings = false;
            //check with the user to update params
            if (Params.AskAboutSemiTrypsinParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Semi-Trypsin Recommendations?",
                    "We recommend using the following parameters for semi-trypsin searches:\n" +
                    "\t-Check 'Semi-Specific Search'\n" +
                    "\t-Use 'trypsin'\n\n" +
                    "Would you like to use these recommended settings?");

                if (results.UseSettings)
                {
                    useRecommendedSettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    Params.AskAboutSemiTrypsinParams = false;
                    Params.UseSemiTrypsinParams = results.UseSettings;

                    Toml.WriteFile(Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                }
            }
            else if (Params.UseSemiTrypsinParams) //user didn't want to check in, but wanted to update last time
            {
                useRecommendedSettings = true;
            }
            //else do nothing
            return useRecommendedSettings;
        }

        public static List<(string, string)> TopDownModsForGPTMD = new List<(string, string)>
        {
            ("Common Variable", "Oxidation on M"),
            ("Common Biological", "Acetylation on K"),
            ("Common Biological", "Acetylation on X"),          
            ("Common Biological", "Phosphorylation on S"),
            ("Common Biological", "Phosphorylation on T"),
            ("Common Biological", "Phosphorylation on Y"),
        };
    }
}
