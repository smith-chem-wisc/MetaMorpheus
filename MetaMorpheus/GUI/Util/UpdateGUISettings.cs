using System.Collections.Generic;

namespace MetaMorpheusGUI
{
    public static class UpdateGUISettings
    {
        public static bool UseNonSpecificRecommendedSettings()
        {
            bool useRecommendedSettings = false;
            //check with the user to update params
            if (GuiGlobalParamsViewModel.Instance.AskAboutNonSpecificParams)
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
                    GuiGlobalParamsViewModel.Instance.AskAboutNonSpecificParams = false;
                    GuiGlobalParamsViewModel.Instance.UseNonSpecificParams = results.UseSettings;
                }
            }
            else if (GuiGlobalParamsViewModel.Instance.UseNonSpecificParams) //user didn't want to check in, but wanted to update last time
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
            if (GuiGlobalParamsViewModel.Instance.AskAboutTopDownParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Top-Down Recommendations?",
                    "We recommend using the following parameters for top-down searches:\n" +
                        "\t-Uncheck 'Use Provided Precursor'\n" +
                        "\t-Use '60' for 'Deconvolution Max Assumed Charge State'\n" +
                        "\t-Uncheck 'Trim MS2 Peaks'\n" +
                        "\t-Uncheck all variable mods (Please use a GPTMD database instead)\n" +
                        "\t-GPTMD TASK ONLY: Search for only acetylation, phosphorylation, and oxidation of M\n" +
                        "\t-SEARCH TASK ONLY: Check 'No Quantification'\n" +
                        "\t-SEARCH TASK ONLY: Check '1, 2, or 3 Missed Monoisotopic Peaks'\n" +
                        "\t-SEARCH TASK ONLY: Check 'Internal Ions - Min Internal Length 10'\n" +
                    "Would you like to use these recommended settings?");

                if (results.UseSettings)
                {
                    useRecommendedSettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    GuiGlobalParamsViewModel.Instance.AskAboutTopDownParams = false;
                    GuiGlobalParamsViewModel.Instance.UseTopDownParams = results.UseSettings;
                }
            }
            else if (GuiGlobalParamsViewModel.Instance.UseTopDownParams) //user didn't want to check in, but wanted to update last time
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
            if (GuiGlobalParamsViewModel.Instance.AskAboutArgCParams)
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
                    GuiGlobalParamsViewModel.Instance.AskAboutArgCParams = false;
                    GuiGlobalParamsViewModel.Instance.UseArgCParams = results.UseSettings;
                }
            }
            else if (GuiGlobalParamsViewModel.Instance.UseArgCParams) //user didn't want to check in, but wanted to update last time
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
            if (GuiGlobalParamsViewModel.Instance.AskAboutChymotrypsinParams)
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
                    GuiGlobalParamsViewModel.Instance.AskAboutChymotrypsinParams = false;
                    GuiGlobalParamsViewModel.Instance.UseChymotrypsinParams = results.UseSettings;
                }
            }
            else if (GuiGlobalParamsViewModel.Instance.UseChymotrypsinParams) //user didn't want to check in, but wanted to update last time
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
            if (GuiGlobalParamsViewModel.Instance.AskAboutElastaseParams)
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
                    GuiGlobalParamsViewModel.Instance.AskAboutElastaseParams = false;
                    GuiGlobalParamsViewModel.Instance.UseElastaseParams = results.UseSettings;
                }
            }
            else if (GuiGlobalParamsViewModel.Instance.UseElastaseParams) //user didn't want to check in, but wanted to update last time
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
            if (GuiGlobalParamsViewModel.Instance.AskAboutSemiTrypsinParams)
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
                    GuiGlobalParamsViewModel.Instance.AskAboutSemiTrypsinParams = false;
                    GuiGlobalParamsViewModel.Instance.UseSemiTrypsinParams = results.UseSettings;
                }
            }
            else if (GuiGlobalParamsViewModel.Instance.UseSemiTrypsinParams) //user didn't want to check in, but wanted to update last time
            {
                useRecommendedSettings = true;
            }
            //else do nothing
            return useRecommendedSettings;
        }

        public static bool UseSpectralRecoveryMandatorySettings()
        {
            bool useMandatorySettings = false;
            //check with the user to update params
            if (GuiGlobalParamsViewModel.Instance.AskAboutSpectralRecoveryParams)
            {
                var results = ProteaseSpecificMsgBox.Show("Use Spectral Recovery Settings?",
                    "The following parameters are necessary for the Spectral Recovery algorithm:\n" +
                    "\t-Check 'Match between runs' (Search Task Only)\n" +
                    "\t-Check 'Write Spectral Library' (Search Task Only)\n" +
                    "\t-SEARCH TASK ONLY: Increase the maximum allowed modified isoforms to 4096'\n" +
                    "Would you like to use these settings?");

                if (results.UseSettings)
                {
                    useMandatorySettings = true;
                }
                //else do nothing

                //if they don't want to see this again, save the answer
                if (!results.AskAgain)
                {
                    GuiGlobalParamsViewModel.Instance.AskAboutSpectralRecoveryParams = false;
                    GuiGlobalParamsViewModel.Instance.UseSpectralRecoveryParams = results.UseSettings;
                }
            }
            else if (GuiGlobalParamsViewModel.Instance.UseSpectralRecoveryParams) //user didn't want to check in, but wanted to update last time
            {
                useMandatorySettings = true;
            }
            //else do nothing
            return useMandatorySettings;
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
