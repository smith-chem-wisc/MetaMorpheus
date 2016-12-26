namespace MetaMorpheus
{
    public class RawDataAndResults
    {
        public bool Use { get; set; }
        public string FileName { get; private set; }
        public bool UseMzid { get; set; }
        public string mzidName { get; private set; }
        public bool UsePsmsTsv { get; set; }
        public string psmsTSVName { get; private set; }

        public RawDataAndResults(string FileName, string mzidName, string psmsTSVName)
        {
            this.FileName = FileName;
            if (FileName != null)
                Use = true;
            this.mzidName = mzidName;
            if (mzidName != null)
                UseMzid = true;
            this.psmsTSVName = psmsTSVName;
            if (psmsTSVName != null)
                UsePsmsTsv = true;
        }

        public void AddFilePath(string FileName)
        {
            this.FileName = FileName;
            Use = true;
        }

        public void AddMZID(string mzidName)
        {
            this.mzidName = mzidName;
            UseMzid = true;
        }

        public void AddTSV(string psmsTSVName)
        {
            this.psmsTSVName = psmsTSVName;
            UsePsmsTsv = true;
        }
    }
}