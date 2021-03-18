using System.Text;

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
            htmlQueryString.Append("https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " ");

            if (!xmlFormat) //fasta
            {
                htmlQueryString.Append("&format=fasta");
            }
            else
            {
                htmlQueryString.Append("&format=xml");
            }

            if (reviewed)
            {
                htmlQueryString.Append("&reviewed:yes");
            }

            if (isoforms)
            {
                htmlQueryString.Append("&include:yes");
            }
            else
            {
                htmlQueryString.Append("&include:no");
            }

            if (compressed)
            {
                htmlQueryString.Append("&compress=yes");
            }
            else
            {
                htmlQueryString.Append("&compress=no");
            }

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

            if (isoforms)
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