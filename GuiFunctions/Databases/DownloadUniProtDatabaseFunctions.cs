using Proteomics.ProteolyticDigestion;
using System.Text;
using static Org.BouncyCastle.Bcpg.Attr.ImageAttrib;

namespace GuiFunctions.Databases
{
    /// <summary>
    /// We place functions from the GUI in the class to enable testing.
    /// </summary>

    public class DownloadUniProtDatabaseFunctions
    {
        public static string GetUniProtHtmlQueryString(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed)
        {
            StringBuilder htmlQueryString = new StringBuilder();
            htmlQueryString.Append("https://rest.uniprot.org/uniprotkb/stream?");

            if (!xmlFormat) //fasta
            {
                htmlQueryString.Append("format=fasta");
            }
            else
            {
                htmlQueryString.Append("format=xml");
            }

            if (isoforms && !xmlFormat)
            {
                htmlQueryString.Append("&includeIsoform=true");
            }

            if (reviewed)
            {
                htmlQueryString.Append("&reviewed=true");
            }

            if (compressed)
            {
                htmlQueryString.Append("&compressed=true");
            }
            htmlQueryString.Append("&query=%28%28proteome%3A" + proteomeID + "%29%29");
            return htmlQueryString.ToString();
        }

        public static string GetUniprotFilename(string uniprotProteomeId, bool reviewed, bool isoforms, bool xmlFormat, bool compressed)
        {
            string filename = "\\" + uniprotProteomeId;
            if (reviewed)
            {
                filename += "_reviewed";
            }
            else
            {
                filename += "_withUnreviewed";
            }

            if (isoforms && !xmlFormat)
            {
                filename += "_isoform";
            }

            if (xmlFormat)
            {
                filename += ".xml";
            }
            else
            {
                filename += ".fasta";
            }

            if (compressed)
            {
                filename += ".gz";
            }

            return filename;
        }
    }
}