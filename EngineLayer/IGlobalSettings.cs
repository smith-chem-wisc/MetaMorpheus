namespace EngineLayer
{
    public interface IGlobalSettings
    {
        bool WriteExcelCompatibleTSVs { get; }
        bool UserHasAgreedToThermoRawFileReaderLicence { get; }
    }
}