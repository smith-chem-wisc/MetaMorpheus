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
            StringBuilder sb = new StringBuilder();
            sb.Append("https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " ");

            if (!xmlFormat) //fasta
            {
                sb.Append("&format=fasta");
                if (reviewed)
                {
                    sb.Append("&reviewed:yes");
                    if (isoforms)
                    {
                        sb.Append("&include:yes");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                    else//no isoforms
                    {
                        sb.Append("&include:no");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                }
                else //revied and unreviewed both
                {
                    //we don't need to append anything if we want both reviewed and unreviewed
                    if (isoforms)
                    {
                        sb.Append("&include:yes");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                    else//no isoforms
                    {
                        sb.Append("&include:no");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                }
            }
            else //xml -- currently isoforms are automatically included (I think)
            {
                sb.Append("&format=xml");
                if (reviewed)
                {
                    sb.Append("&reviewed:yes");
                    if (isoforms)
                    {
                        sb.Append("&include:yes");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                    else//no isoforms
                    {
                        sb.Append("&include:no");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                }
                else //revied and unreviewed both
                {
                    //we don't need to append anything if we want both reviewed and unreviewed
                    if (isoforms)
                    {
                        sb.Append("&include:yes");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                    else//no isoforms
                    {
                        sb.Append("&include:no");
                        if (compressed)
                        {
                            sb.Append("&compress=yes");
                        }
                        else //not compressed
                        {
                            sb.Append("&compress=no");
                        }
                    }
                }
            }

            return sb.ToString();
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