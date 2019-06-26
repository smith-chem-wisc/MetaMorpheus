namespace EngineLayer
{
    public class GlobalSettings : IGlobalSettings
    {
        public bool WriteExcelCompatibleTSVs { get; set; }
        public bool UserHasAgreedToThermoRawFileReaderLicence { get; set; }
    }
}