using System.Linq;
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
            htmlQueryString.Append("https://rest.uniprot.org/uniprotkb/search?");

            string[] queryArray = new string[4];

            queryArray[0] = compressed ? "compressed=true" : "";
            queryArray[1] = xmlFormat ? "format=xml" : "format=fasta";
            queryArray[2] = isoforms && !xmlFormat ? "includeIsoform=true" : "";
            queryArray[3] = reviewed ? $"query=reviewed:true&proteome:{proteomeID}" : $"query=proteome:{proteomeID}";

            string[] queryArrayReduced = queryArray.Where(x => x != "").ToArray() ;

            htmlQueryString.Append(string.Join("&", queryArrayReduced));

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