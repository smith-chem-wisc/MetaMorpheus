using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace MetaMorpheusGUI
{
    class GlycanStructureAnnotation
    {
        public static void DrawGlycan(Canvas canvas, string glycanStruct, double scale)
        {
            GlycanNode rootNode = GlycanNode.ReadGlycan2GlycanNode(glycanStruct);

            Queue<GlycanNode> que = new Queue<GlycanNode>();
            GlycanNode temp = null;
            GlycanNode last;

            if (rootNode != null)
            {
                que.Enqueue(rootNode);
            }
            while(que.Count != 0)
            {
                last = temp;
                temp = que.Dequeue();
                temp.canChange = true;

                //Add function here              
                temp.point = GetInitLocation(temp, scale);
                ChangeLocations(temp, last, scale);


                if (temp.LeftChild != null)
                {
                    que.Enqueue(temp.LeftChild);
                }
                if (temp.RightChild != null)
                {
                    que.Enqueue(temp.RightChild);
                }
            }

            //Change rootNode location if 'F' at level 1.
            ChangeFLocation(rootNode, scale);


            //Iterate again to draw the Lines
            if (rootNode != null)
            {
                que.Enqueue(rootNode);
            }
            while (que.Count != 0)
            {
                temp = que.Dequeue();

                //Add function here              
                DrawLine(canvas, temp);
                DrawShape(canvas, temp);

                if (temp.LeftChild != null)
                {
                    que.Enqueue(temp.LeftChild);
                }
                if (temp.RightChild != null)
                {
                    que.Enqueue(temp.RightChild);
                }
            }

            DrawRootLine(canvas, rootNode);
        }

        public static void ChangeFLocation(GlycanNode rootNode, double scale)
        {
            if (rootNode.LeftChild != null && rootNode.LeftChild.Value == 'F' && rootNode.LeftChild.LeftChild == null && rootNode.RightChild != null)
            {             
                rootNode.point = new Point(rootNode.RightChild.point.X, rootNode.point.Y);
                rootNode.LeftChild.point = new Point(rootNode.RightChild.point.X + scale, rootNode.point.Y);
            }
            if (rootNode.RightChild != null && rootNode.RightChild.Value == 'F' && rootNode.RightChild.LeftChild == null && rootNode.LeftChild != null)
            {
                rootNode.point = new Point(rootNode.LeftChild.point.X, rootNode.point.Y);
                rootNode.RightChild.point = new Point(rootNode.LeftChild.point.X + scale, rootNode.point.Y);
            }
        }

        public static Point GetInitLocation(GlycanNode glycanNode, double scale)
        {
            Point point = new Point(20, 20);
            if (glycanNode.Father != null)
            {
                if (glycanNode == glycanNode.Father.LeftChild && glycanNode.Father.RightChild == null)
                {
                    point = new Point(glycanNode.Father.point.X, glycanNode.Father.point.Y - scale);
                }
                if (glycanNode == glycanNode.Father.LeftChild && glycanNode.Father.RightChild != null)
                {
                    point = new Point(glycanNode.Father.point.X - scale/2, glycanNode.Father.point.Y - scale);
                }
                if (glycanNode == glycanNode.Father.RightChild)
                {
                    point = new Point(glycanNode.Father.point.X + scale/2, glycanNode.Father.point.Y - scale);
                }              
            }

            return point;
        }

        public static void ChangeLocations(GlycanNode newNode, GlycanNode lastNode, double scale)
        {

            if (newNode.point.X < 20)
            {
                var rootNode = GobackRootNode(newNode);
                ChangeAllNodeLocation(rootNode, scale/2, 0);
            }
            if (newNode.point.Y < 20)
            {
                var rootNode = GobackRootNode(newNode);
                ChangeAllNodeLocation(rootNode, 0, scale);
            }

            if (lastNode!= null && lastNode.point.Y == newNode.point.Y && newNode.point.X - lastNode.point.X  < scale )
            {
                newNode.point = new Point(lastNode.point.X + scale, newNode.point.Y);              

                Point newNodeFatherPoint = CalNewFatherLocation(newNode, scale);               

                double changedX = newNodeFatherPoint.X - newNode.Father.point.X;

                ChangeFatherAndRightChildrenLocation(newNode.Father, changedX);
            }          
        }

        private static void ChangeFatherAndRightChildrenLocation(GlycanNode glycanNode, double changedX)
        {
            glycanNode.point = new Point(glycanNode.point.X + changedX, glycanNode.point.Y);

            if (glycanNode.Father!=null)
            {
                if (glycanNode == glycanNode.Father.LeftChild)
                {
                    ChangeFatherAndRightChildrenLocation(glycanNode.Father, changedX);
                    if( glycanNode.Father.RightChild != null)
                    {
                        ChangeAllNodeLocation(glycanNode.Father.RightChild, changedX, 0);
                    }
                }


                if (glycanNode == glycanNode.Father.RightChild)
                {
                    changedX = changedX / 2;
                    ChangeFatherAndRightChildrenLocation(glycanNode.Father, changedX);
                }
            }
        }

        private static void ChangeAllNodeLocation(GlycanNode glycanNode, double changedX, double changedY)
        {
            if (glycanNode != null)
            {
                if (glycanNode.canChange)
                {
                    glycanNode.point = new Point(glycanNode.point.X + changedX, glycanNode.point.Y + changedY);
                }              
                ChangeAllNodeLocation(glycanNode.LeftChild, changedX, changedY);
                ChangeAllNodeLocation(glycanNode.RightChild, changedX, changedY);
            }

        }

        public static GlycanNode GobackRootNode(GlycanNode node)
        {
            while (node.Father != null)
            {
                node = node.Father;
            }
            return node;
        }

        private static Point CalNewFatherLocation(GlycanNode glycanNode, double scale)
        {
            Point point = new Point(20, 20);
            if (glycanNode.Father != null)
            {
                if (glycanNode == glycanNode.Father.LeftChild && glycanNode.Father.RightChild == null)
                {
                    point = new Point(glycanNode.point.X, glycanNode.point.Y + scale);
                }
                if (glycanNode == glycanNode.Father.LeftChild && glycanNode.Father.RightChild != null)
                {
                    point = new Point(glycanNode.point.X + scale/2, glycanNode.point.Y + scale);
                }
                if (glycanNode == glycanNode.Father.RightChild)
                {
                    point = new Point(glycanNode.point.X - scale/2, glycanNode.point.Y + scale);
                }
            }
            return point;
        }

        public static void DrawLine(Canvas canvas, GlycanNode glycanNode)
        {
            if (glycanNode.Father!=null)
            {
                Line line = new Line();
                line.X1 = glycanNode.Father.point.X;
                line.Y1 = glycanNode.Father.point.Y;
                line.X2 = glycanNode.point.X;
                line.Y2 = glycanNode.point.Y;
                line.Stroke = Brushes.Black;
                line.StrokeThickness = 2;
                Panel.SetZIndex(line, 1);
                canvas.Children.Add(line);
            }
        }

        public static void DrawShape(Canvas canvas, GlycanNode glycanNode)
        {
            SolidColorBrush colorH = new SolidColorBrush(color: Color.FromArgb(255, 0, 166, 81));
            SolidColorBrush colorN = new SolidColorBrush(color: Color.FromArgb(255, 0, 144, 188));
            SolidColorBrush colorA = new SolidColorBrush(color: Color.FromArgb(255, 165, 67, 153));
            SolidColorBrush colorG = new SolidColorBrush(color: Color.FromArgb(255, 143, 204, 233));
            SolidColorBrush colorF = new SolidColorBrush(color: Color.FromArgb(255, 237, 28, 36));
            var point = glycanNode.point;
            PointCollection myPointCollection = new PointCollection();
            Polygon polygon = new Polygon();

            switch (glycanNode.Value)
            {
                case 'H':
                    Ellipse circle = new Ellipse()
                    {
                        Width = 20,
                        Height = 20,
                        Stroke = Brushes.Black,
                        StrokeThickness = 0.5,
                        Fill = colorH
                    };
                    circle.Margin = new Thickness(point.X-10, point.Y-10, 0, 0);
                    Panel.SetZIndex(circle, 2);
                    canvas.Children.Add(circle);
                    break;

                case 'N':
                    myPointCollection = new PointCollection()
                    {
                        new Point(point.X - 10, point.Y + 10), new Point(point.X + 10, point.Y + 10), new Point(point.X + 10, point.Y - 10), new Point(point.X - 10, point.Y - 10), new Point(point.X - 10, point.Y + 10)
                    };
                    polygon = new Polygon()
                    {
                        Points = myPointCollection,
                        Stroke = Brushes.Black,
                        StrokeThickness = 0.5,
                        Fill = colorN
                    };
                    Panel.SetZIndex(polygon, 2);
                    canvas.Children.Add(polygon);
                    break;

                case 'A':
                    myPointCollection = new PointCollection()
                    {
                        new Point(point.X - 10, point.Y ), new Point(point.X, point.Y + 10), new Point(point.X + 10, point.Y), new Point(point.X, point.Y - 10), new Point(point.X - 10, point.Y )
                    };
                    polygon = new Polygon()
                    {
                        Points = myPointCollection,
                        Stroke = Brushes.Black,
                        StrokeThickness = 0.5,
                        Fill = colorA
                    };
                    Panel.SetZIndex(polygon, 2);
                    canvas.Children.Add(polygon);
                    break;

                case 'G':
                    myPointCollection = new PointCollection()
                    {
                        new Point(point.X - 10, point.Y ), new Point(point.X, point.Y + 10), new Point(point.X + 10, point.Y), new Point(point.X, point.Y - 10), new Point(point.X - 10, point.Y )
                    };
                    polygon = new Polygon()
                    {
                        Points = myPointCollection,
                        Stroke = Brushes.Black,
                        StrokeThickness = 0.5,
                        Fill = colorG                
                    };
                    Panel.SetZIndex(polygon, 2);
                    canvas.Children.Add(polygon);
                    break;

                case 'F':
                    myPointCollection = new PointCollection()
                    {
                        new Point(point.X - 10, point.Y + 10), new Point(point.X + 10, point.Y + 10), new Point(point.X, point.Y - 8), new Point(point.X - 10, point.Y + 10)
                    };
                    polygon = new Polygon()
                    {
                        Points = myPointCollection,
                        Stroke = Brushes.Black,
                        StrokeThickness = 0.5,
                        Fill = colorF
                    };
                    Panel.SetZIndex(polygon, 2);
                    canvas.Children.Add(polygon);
                    break;
            }
        }

        public static void DrawRootLine(Canvas canvas, GlycanNode glycanNode)
        {
            Line line = new Line();
            line.X1 = glycanNode.point.X;
            line.Y1 = glycanNode.point.Y;
            line.X2 = glycanNode.point.X;
            line.Y2 = glycanNode.point.Y + 50;
            line.Stroke = Brushes.Black;
            line.StrokeThickness = 2;
            Panel.SetZIndex(line, 1);
            canvas.Children.Add(line);

            Line line2 = new Line();
            line2.X1 = glycanNode.point.X - 50;
            line2.Y1 = glycanNode.point.Y + 50;
            line2.X2 = glycanNode.point.X + 50;
            line2.Y2 = glycanNode.point.Y + 50;
            line2.Stroke = Brushes.Black;
            line2.StrokeThickness = 4;
            Panel.SetZIndex(line2, 1);

            canvas.Children.Add(line2);
        }
    }

    class GlycanNode 
    {
        private int? level;

        public GlycanNode(char v, int l)
        {
            Value = v;
            LeftChild = null;
            RightChild = null;
            Father = null;
            level = l;
            canChange = false;
        }

        public char Value { get; }
        public GlycanNode Father { get; private set; }
        public GlycanNode LeftChild { get; private set; }
        public GlycanNode RightChild { get; private set; }
        public int Level { get { return level.Value; } }

        public Point point { get; set; }
        public bool canChange { get; set; }

        public static GlycanNode ReadGlycan2GlycanNode (string theGlycanStruct)
        {
            int level = 0;
            GlycanNode curr = new GlycanNode(theGlycanStruct[1], level);
            for (int i = 2; i < theGlycanStruct.Length - 1; i++)
            {
                if (theGlycanStruct[i] == '(')
                {
                    continue;
                }
                if (theGlycanStruct[i] == ')')
                {
                    curr = curr.Father;
                    level--;
                }
                else
                {
                    level++;
                    if (curr.LeftChild == null)
                    {
                        curr.LeftChild = new GlycanNode(theGlycanStruct[i], level);
                        curr.LeftChild.Father = curr;
                        curr = curr.LeftChild;
                    }
                    else
                    {
                        curr.RightChild = new GlycanNode(theGlycanStruct[i], level);
                        curr.RightChild.Father = curr;
                        curr = curr.RightChild;
                    }
                }
            }
            return curr;
        }

    }
}

