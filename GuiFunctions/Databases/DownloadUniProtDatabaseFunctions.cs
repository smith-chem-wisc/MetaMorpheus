using System;
using System.Collections.Generic;
using System.Text;
using static UsefulProteomicsDatabases.ProteinDbRetriever;

namespace GuiFunctions.Databases
{
    /// <summary>
    /// We place functions from the GUI in the class to enable testing.
    /// </summary>
    
    public class DownloadUniProtDatabaseFunctions
    {

        public static string GetUniProtHtmlQueryString(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed)
        {
            string htmlQueryString;
            if (!xmlFormat) //fasta
            {
                if (reviewed)
                {
                    if (isoforms)
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=yes" + "&format=fasta" + "&include:yes";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=no" + "&format=fasta" + "&include:yes";
                        }
                    }
                    else//no isoforms
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=yes" + "&format=fasta" + "&include:no";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=no" + "&format=fasta" + "&include:no";
                        }
                    }
                }
                else //revied and unreviewed both
                {
                    if (isoforms)
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=yes" + "&format=fasta" + "&include:yes";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=no" + "&format=fasta" + "&include:yes";
                        }
                    }
                    else//no isoforms
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=yes" + "&format=fasta" + "&include:no";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=no" + "&format=fasta" + "&include:no";
                        }
                    }
                }
            }
            else //xml -- currently isoforms are automatically included (I think)
            {
                if (reviewed)
                {
                    if (isoforms)
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=yes" + "&format=xml" + "&include:yes";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=no" + "&format=xml" + "&include:yes";
                        }
                    }
                    else//no isoforms
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=yes" + "&format=xml" + "&include:no";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:yes" + "&compress=no" + "&format=xml" + "&include:no";
                        }
                    }
                }
                else //revied and unreviewed both
                {
                    if (isoforms)
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=yes" + "&format=xml" + "&include:yes";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=no" + "&format=xml" + "&include:yes";
                        }
                    }
                    else//no isoforms
                    {
                        if (compressed)
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=yes" + "&format=xml" + "&include:no";
                        }
                        else //not compressed
                        {
                            htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=no" + "&format=xml" + "&include:no";
                        }
                    }
                }
            }

            return htmlQueryString;
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
