using iTextSharp.text;
using System.IO;
using iTextSharp.text.pdf;
using System.Windows;
using System.Windows.Media.Imaging;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    class CustomPdfWriter
    {
        public static void WriteToPdfMetaDraw(double wid, double hei, string pathToImage, string pathToLegend, string filePath)
        {
            Document doc = new Document(new Rectangle((float)wid, (float)hei));
            PdfWriter.GetInstance(doc, new FileStream(filePath, FileMode.Create));
            doc.Open();

            Image coverageMap = Image.GetInstance(pathToImage);
            coverageMap.ScalePercent(70f);
            coverageMap.Alignment = Image.ALIGN_CENTER;
            coverageMap.SetAbsolutePosition(((float)wid - coverageMap.ScaledWidth) / 2, ((float)hei - coverageMap.ScaledHeight) / 2);
            doc.Add(coverageMap);

            Image legend = Image.GetInstance(pathToLegend);
            legend.ScalePercent(80f);
            legend.Alignment = Image.ALIGN_CENTER;
            legend.SetAbsolutePosition(((float)wid - coverageMap.ScaledWidth) / 2 + 200, 10);
            doc.Add(legend);
            doc.Close();

            File.Delete(pathToImage);
            File.Delete(pathToLegend);
        }

        // renders a bitmap image and save as a png file
        public static void RenderImage(int width, int height, string path, System.Windows.Controls.Canvas canvas)
        {
            canvas.Measure(new Size(width, height));
            canvas.Arrange(new Rect(new Size(width, height)));

            RenderTargetBitmap renderBitmap = new RenderTargetBitmap(width, height, 96, 96, PixelFormats.Pbgra32);

            renderBitmap.Render(canvas);
            PngBitmapEncoder encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(renderBitmap));

            using (FileStream file = File.Create(path))
            {
                encoder.Save(file);
            }
        }
    }
}
