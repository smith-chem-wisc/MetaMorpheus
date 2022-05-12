using System.Text;
using EngineLayer;
using FlashLFQ;

namespace TaskLayer
{
    public class MbrWriter
    {
        public MbrWriter()
        {

        }
        public static string TabSeparatedHeader
        {
            
            get
            {
                var sb = new StringBuilder();
                sb.Append(PeptideSpectralMatch.GetTabSeparatedHeader());
                sb.Append('\t');
                sb.Append(ChromatographicPeak.TabSeparatedHeader);
                return sb.ToString();
            }
        }

        public static string MbrPsmOutputString(ChromatographicPeak peak, PeptideSpectralMatch psm)
        {
            var sb = new StringBuilder();
            sb.Append(psm.ToString());
            sb.Append('\t');
            sb.Append(peak.ToString());
            return sb.ToString();
        }

    }
}
