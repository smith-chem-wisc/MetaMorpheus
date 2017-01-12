using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class ModernSearchResults : MyResults
    {
        private int numMS2spectra;
        private int[] numMS2spectraMatched;

        public List<ModernSpectrumMatch>[] newPsms { get; private set; }

        public ModernSearchResults(List<ModernSpectrumMatch>[] newPsms, int numMS2spectra, int[] numMS2spectraMatched, ModernSearchEngine s) : base(s)
        {
            this.numMS2spectra = numMS2spectra;
            this.numMS2spectraMatched = numMS2spectraMatched;
            this.newPsms = newPsms;
        }

        protected override string GetStringForOutput()
        {
            var sp = (ModernSearchEngine)s;
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("\t\tTotal ms2 spectra seen: " + numMS2spectra);
            sb.Append("\t\t"+string.Join(Environment.NewLine, sp.searchModes.Zip(numMS2spectraMatched, (a, b) => a.FileNameAddition + " : " + b)));

            return sb.ToString();
        }
    }
}