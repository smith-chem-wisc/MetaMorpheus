using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace MetaMorpheusGUI
{
    class BaseDraw
    {
        /// <summary>
        /// Draw the line seperator @ top
        /// </summary>
        /// <param name="cav">Canvas to where draws the line</param>
        /// <param name="topLoc">Location of the starting point on top (exp: 0,0)</param>
        public static void topSplittingDrawing(Canvas cav, Point topLoc, Color clr, string footnote)
        {
            double x = topLoc.X, y = topLoc.Y;
            Polyline bot = new Polyline();
            bot.Points = new PointCollection() { new Point(x + 10, y), new Point(x, y + 10), new Point(x, y + 40) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 2;
            cav.Children.Add(bot);
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
            TextBlock tb = new TextBlock();
            tb.Foreground = new SolidColorBrush(clr);
            tb.Text = footnote;
            tb.FontSize = 10;
            Canvas.SetTop(tb, y - 10);
            Canvas.SetLeft(tb, x + 10);
            Canvas.SetZIndex(tb, 2);
            cav.Children.Add(tb);
        }

        /// <summary>
        /// Draw the line seperator @ bottom
        /// </summary>
        /// <param name="cav">Canvas to where draws the line</param>
        /// <param name="botLoc">Location of the starting point on bottom (exp: 0,50)</param>
        public static void botSplittingDrawing(Canvas cav, Point botLoc, Color clr, string footnote)
        {
            double x = botLoc.X, y = botLoc.Y;
            Polyline bot = new Polyline();
            bot.Points = new PointCollection() { new Point(x - 10, y), new Point(x, y - 10), new Point(x, y - 40) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 2;
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
            cav.Children.Add(bot);
            TextBlock tb = new TextBlock();
            tb.Foreground = new SolidColorBrush(clr);
            tb.Text = footnote;
            tb.FontSize = 10;
            Canvas.SetTop(tb, y - 8);
            Canvas.SetLeft(tb, x - 22);
            Canvas.SetZIndex(tb, 2);
            cav.Children.Add(tb);
        }

        /// <summary>
        /// Create text blocks on canvas
        /// </summary>
        /// <param name="cav"> Canvas Board </param>
        /// <param name="loc"> Provate the (x,y) coordinates for textblock</param>
        /// <param name="txt"> Message for textblock</param>
        /// <returns> the width of current addup</returns>
        public static void txtDrawing(Canvas cav, Point loc, string txt, Brush clr)
        {
            TextBlock tb = new TextBlock();
            tb.Foreground = clr;
            tb.Text = txt;
            tb.Height = 30;
            tb.FontSize = 25;
            tb.FontWeight = FontWeights.Bold;
            tb.FontFamily = new FontFamily("Courier New"); // monospaced font

            Canvas.SetTop(tb, loc.Y);
            Canvas.SetLeft(tb, loc.X);
            Panel.SetZIndex(tb, 2); //lower priority
            cav.Children.Add(tb);
            cav.UpdateLayout();
        }

        /// <summary>
        /// Display texts surounded by circle
        /// </summary>
        /// <param name="cav"></param>
        /// <param name="loc"></param>
        /// <param name="txt"></param>
        /// <param name="clr"></param>
        /// <returns></returns>
        public static void circledTxtDraw(Canvas cav, Point loc, SolidColorBrush clr)
        {
            Ellipse circle = new Ellipse()
            {
                Width = 24,
                Height = 24,
                Stroke = clr,
                StrokeThickness = 1,
                Fill = clr,
                Opacity = 0.7
            };
            Canvas.SetLeft(circle, loc.X);
            Canvas.SetTop(circle, loc.Y);
            Panel.SetZIndex(circle, 1);
            cav.Children.Add(circle);
        }

        /// <summary>
        /// Clear canvas board
        /// </summary>
        /// <param name="cav">board to clear</param>
        public static void clearCanvas(Canvas cav)
        {
            cav.Children.Clear();
        }
    }
}
