using iTextSharp.text;
using System.IO;
using iTextSharp.text.pdf;

namespace MetaMorpheusGUI
{
    class PdfWriter
    {
        public static void WriteToPdf(System.Drawing.Image img, double wid, double hei, string filePath)
        {
            Document doc = new Document(new Rectangle((float)wid, (float)hei));
            iTextSharp.text.pdf.PdfWriter.GetInstance(doc, new FileStream(filePath, FileMode.Create));
            doc.Open();
            Image pdfImage = Image.GetInstance(img, System.Drawing.Imaging.ImageFormat.Bmp);
            pdfImage.SetAbsolutePosition(0, 0);
            doc.Add(pdfImage);
            doc.Close();
        }
    }
}
