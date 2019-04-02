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
            bot.Points = new PointCollection() { new Point(x + 10, y + 10), new Point(x, y + 10), new Point(x, y + 24) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 1;
            cav.Children.Add(bot);
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
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
            bot.Points = new PointCollection() { new Point(x - 10, y - 10), new Point(x, y - 10), new Point(x, y - 24) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 1;
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
            cav.Children.Add(bot);
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
            tb.FontFamily = new FontFamily("Arial"); // monospaced font

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

        public static void DrawCrosslinker(Canvas cav, Point alphaBotLoc, Point betaBotLoc, Color clr)
        {
            Polyline bot = new Polyline();
            double distance = (betaBotLoc.Y - alphaBotLoc.Y) / 2;
            bot.Points = new PointCollection() { alphaBotLoc, new Point(alphaBotLoc.X, alphaBotLoc.Y + distance), new Point(betaBotLoc.X, alphaBotLoc.Y + distance), betaBotLoc };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 2;
            cav.Children.Add(bot);
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