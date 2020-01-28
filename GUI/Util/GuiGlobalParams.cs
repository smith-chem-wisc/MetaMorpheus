namespace MetaMorpheusGUI
{
    public class GuiGlobalParams
    {
        public bool AskAboutUpdating { get; internal set; } = true;
        public bool AskBeforeExitingMetaMorpheus { get; internal set; } = true; 

        //Ask about protease-specific parameter recommendations
        public bool AskAboutTopDownParams { get; internal set; } = true;
        public bool AskAboutChymotrypsinParams { get; internal set; } = true;
        public bool AskAboutElastaseParams { get; internal set; } = true;
        public bool AskAboutNonSpecificParams { get; internal set; } = true;
        public bool AskAboutSemiTrypsinParams { get; internal set; } = true;
        public bool AskAboutArgCParams { get; internal set; } = true;

        //Use protease-specific parameter recommendations
        public bool UseTopDownParams { get; internal set; } = true;
        public bool UseChymotrypsinParams { get; internal set; } = true;
        public bool UseElastaseParams { get; internal set; } = true;
        public bool UseNonSpecificParams { get; internal set; } = true;
        public bool UseSemiTrypsinParams { get; internal set; } = true;
        public bool UseArgCParams { get; internal set; } = true;
    }
}