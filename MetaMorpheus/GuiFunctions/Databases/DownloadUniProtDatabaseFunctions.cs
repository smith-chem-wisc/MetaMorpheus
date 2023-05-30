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

//if (compressed)
//{
//    htmlQueryString.Append("compressed=true");

//    if (!xmlFormat && reviewed && isoforms)
//    {
//        htmlQueryString.Append("&format=fasta&reviewed=true&includeIsoform=true");
//    }
//    else if (xmlFormat && reviewed && isoforms)
//    {
//        htmlQueryString.Append("&format=xml&reviewed=true");
//    }
//    else if (!xmlFormat && !reviewed && isoforms)
//    {

//    }
//}


//if (compressed)
//{
//    htmlQueryString.Append("compressed=true");
//}

//if (!xmlFormat && compressed) //fasta
//{
//    htmlQueryString.Append("&format=fasta");
//}
//else if (!xmlFormat && !compressed)
//{
//    htmlQueryString.Append("format=fasta");
//}
//else if (xmlFormat && compressed)
//{
//    htmlQueryString.Append("&format=xml");
//}
//else
//{
//    htmlQueryString.Append("format=xml");
//}

//if (isoforms && !xmlFormat) //Only the .fasta file can be isoforms
//{
//    htmlQueryString.Append("&includeIsoform=true");
//}

//// Query 

//htmlQueryString.Append("&query=");

//if (reviewed)
//{
//    htmlQueryString.Append("reviewed:true&");
//}

//htmlQueryString.Append("proteome:" + proteomeID);