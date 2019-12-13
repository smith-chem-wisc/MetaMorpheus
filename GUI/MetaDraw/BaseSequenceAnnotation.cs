using Proteomics;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;
using System;
using System.Windows.Input;
using System.Windows.Media.Animation;

namespace MetaMorpheusGUI
{
    class BaseSequenceAnnotation
    {
        /// <summary>
        /// Draw the line seperator @ top
        /// </summary>
        /// <param name="cav">Canvas to where draws the line</param>
        /// <param name="topLoc">Location of the starting point on top (exp: 0,0)</param>
        public static void DrawCTerminalIon(Canvas cav, Point topLoc, Color clr, string annotationString)
        {
            double x = topLoc.X, y = topLoc.Y;
            Polyline annotation = new Polyline();
            annotation.Points = new PointCollection() { new Point(x + 10, y + 10), new Point(x, y + 10), new Point(x, y + 24) };
            annotation.Stroke = new SolidColorBrush(clr);
            annotation.StrokeThickness = 1.5;

            annotation.ToolTip = annotationString;

            cav.Children.Add(annotation);
            Canvas.SetZIndex(annotation, 1); //on top of any other things in canvas
        }

        /// <summary>
        /// Draw the line seperator @ bottom
        /// </summary>
        /// <param name="cav">Canvas to where draws the line</param>
        /// <param name="botLoc">Location of the starting point on bottom (exp: 0,50)</param>
        public static void DrawNTerminalIon(Canvas cav, Point botLoc, Color clr, string annotationString)
        {
            double x = botLoc.X, y = botLoc.Y;
            Polyline annotation = new Polyline();
            annotation.Points = new PointCollection() { new Point(x - 10, y - 10), new Point(x, y - 10), new Point(x, y - 24) };
            annotation.Stroke = new SolidColorBrush(clr);
            annotation.StrokeThickness = 1.5;

            annotation.ToolTip = annotationString;

            Canvas.SetZIndex(annotation, 1); //on top of any other things in canvas
            cav.Children.Add(annotation);
        }

        /// <summary>
        /// Create text blocks on canvas
        /// </summary>
        /// <param name="cav"> Canvas Board </param>
        /// <param name="loc"> Provate the (x,y) coordinates for textblock</param>
        /// <param name="txt"> Message for textblock</param>
        /// <returns> the width of current addup</returns>
        public static TextBlock DrawText(Canvas cav, Point loc, string txt, Brush clr)
        {
            TextBlock text = new TextBlock();
            text.Foreground = clr;
            text.Text = txt;
            //tb.Height = 30;
            text.Width = 20;
            text.FontSize = 25;
            text.FontWeight = FontWeights.Medium;
            text.FontFamily = new FontFamily("Arial");
            text.TextAlignment = TextAlignment.Center;
            text.IsHitTestVisible = false;

            Canvas.SetTop(text, loc.Y);
            Canvas.SetLeft(text, loc.X);
            Panel.SetZIndex(text, 2); //lower priority
            cav.Children.Add(text);
            cav.UpdateLayout();

            return text;
        }

        /// <summary>
        /// Display texts surounded by circle
        /// </summary>
        /// <param name="cav"></param>
        /// <param name="loc"></param>
        /// <param name="txt"></param>
        /// <param name="color"></param>
        /// <returns></returns>
        public static void DrawModification(Canvas cav, Point loc, SolidColorBrush color, Modification mod, bool isAmbiguous)
        {
            Rectangle square = new Rectangle()
            {
                Width = 21,
                Height = 22,
                Stroke = Brushes.Black,
                StrokeThickness = 1,
                Fill = color,
                Opacity = 0.8
            };

            if (isAmbiguous)
            {
                square.Fill = GetPatternedTileBrush(color);
            }
            
            square.ToolTip = mod.IdWithMotif;

            if (mod.ChemicalFormula != null)
            {
                square.ToolTip += "\n" + mod.ChemicalFormula.Formula;
            }

            Canvas.SetLeft(square, loc.X);
            Canvas.SetTop(square, loc.Y);
            Panel.SetZIndex(square, 0);
            cav.Children.Add(square);
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

        /// <summary>
        /// Creates a brush to paint half-filled mod annotations, to denote mod localization ambiguity
        /// </summary>
        private static TileBrush GetPatternedTileBrush(SolidColorBrush brushColor)
        {
            PolyLineSegment triangleLinesSegment = new PolyLineSegment();
            triangleLinesSegment.Points.Add(new Point(50, 0));
            triangleLinesSegment.Points.Add(new Point(0, 50));
            PathFigure triangleFigure = new PathFigure();
            triangleFigure.IsClosed = true;
            triangleFigure.StartPoint = new Point(0, 0);
            triangleFigure.Segments.Add(triangleLinesSegment);
            PathGeometry triangleGeometry = new PathGeometry();
            triangleGeometry.Figures.Add(triangleFigure);

            GeometryDrawing triangleDrawing = new GeometryDrawing();
            triangleDrawing.Geometry = triangleGeometry;
            triangleDrawing.Brush = brushColor;
            Pen trianglePen = new Pen(Brushes.Black, 2);
            triangleDrawing.Pen = trianglePen;
            trianglePen.MiterLimit = 0;
            triangleDrawing.Freeze();
            
            DrawingBrush tileBrushWithTiling = new DrawingBrush();
            tileBrushWithTiling.Drawing = triangleDrawing;
            tileBrushWithTiling.TileMode = TileMode.None;
            
            //tileBrushWithTiling.Viewport = new Rect(0, 0, 0.5, 0.5);

            return tileBrushWithTiling;
        }
    }
}