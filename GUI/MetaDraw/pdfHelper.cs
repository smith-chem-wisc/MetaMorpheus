using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using iTextSharp;
using System.Drawing;
using iTextSharp.text;
using System.IO;
using iTextSharp.text.pdf;

namespace MetaMorpheusGUI
{
    class pdfHelper
    {
        public static void saveToPDF(System.Drawing.Image img,double wid, double hei,string filename)
        {
            Document doc = new Document(new iTextSharp.text.Rectangle((float)wid,(float)hei));
            DirectoryInfo di = Directory.CreateDirectory(@".\PDF");
            PdfWriter.GetInstance(doc, new FileStream(Path.Combine(di.ToString(),filename), FileMode.Create));
            doc.Open();
            iTextSharp.text.Image pdfImage = iTextSharp.text.Image.GetInstance(img, System.Drawing.Imaging.ImageFormat.Bmp);
            pdfImage.SetAbsolutePosition(0, 0);
            doc.Add(pdfImage);
            doc.Close();
        }

        //from PDFSharp Sample
        
    }
}
