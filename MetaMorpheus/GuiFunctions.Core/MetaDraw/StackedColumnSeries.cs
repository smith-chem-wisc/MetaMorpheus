using System;
using System.Collections.Generic;
using OxyPlot;
using OxyPlot.Series;

namespace GuiFunctions
{
    /// <summary>
    /// Vertical stacked bars against a CategoryAxis - what ColumnSeries did before OxyPlot 2.1
    /// removed it. Nothing in 2.2 does both: BarSeries stacks but only renders horizontally (it draws
    /// nothing at all against a bottom CategoryAxis), and LinearBarSeries renders vertically but has no
    /// IsStacked. RectangleBarSeries draws arbitrary rectangles, so the stacking is computed here.
    ///
    /// Stack state is shared through the accumulator passed in, so every series in one plot stacks on
    /// the same running totals per category, the way IsStacked did across sibling series.
    /// </summary>
    public class StackedColumnSeries : RectangleBarSeries
    {
        private readonly Dictionary<int, double> _stackTops;
        private readonly double _baseValue;

        public StackedColumnSeries(Dictionary<int, double> stackTops, double baseValue)
        {
            _stackTops = stackTops;
            _baseValue = baseValue;
            StrokeThickness = 0;
        }

        /// <summary>Where the first segment in a category starts. 0.1 under a log y-axis, else 0.</summary>
        public double BaseValue => _baseValue;

        /// <summary>
        /// Gap between adjacent bars, as a fraction of the category slot. Must be set before AddItem,
        /// and must match the CategoryAxis this series is drawn against: BarSeriesBase sized bars as
        /// BarWidth / (1 + GapWidth), and RectangleBarSeries ignores the axis, so the same arithmetic
        /// is reproduced here to keep the bars the width they were.
        /// </summary>
        public double GapWidth { get; set; } = 0.3;

        private double HalfBarWidth => 0.5 / (1 + GapWidth);

        /// <summary>
        /// Widens the reported x-range to whole category slots. A CategoryAxis used to take its range
        /// from the categorized series it was drawn against, giving -0.5 to N-0.5; RectangleBarSeries is
        /// a plain XYAxisSeries, so without this the axis ends at the last bar's right edge instead of
        /// the end of its slot, which shifts every bar and every gridline.
        /// </summary>
        protected override void UpdateMaxMin()
        {
            base.UpdateMaxMin();

            if (Items.Count == 0)
                return;

            MinX = Math.Min(MinX, MinCategoryIndex - 0.5);
            MaxX = Math.Max(MaxX, MaxCategoryIndex + 0.5);
        }

        private int MinCategoryIndex { get; set; } = int.MaxValue;
        private int MaxCategoryIndex { get; set; } = int.MinValue;

        /// <summary>
        /// Adds one segment on top of whatever earlier series already stacked at this category.
        /// </summary>
        public void AddItem(double value, int categoryIndex, string bin, int total, string group = null)
        {
            double y0 = _stackTops.TryGetValue(categoryIndex, out double top) ? top : _baseValue;
            double y1 = y0 + value;
            _stackTops[categoryIndex] = y1;

            MinCategoryIndex = Math.Min(MinCategoryIndex, categoryIndex);
            MaxCategoryIndex = Math.Max(MaxCategoryIndex, categoryIndex);

            Items.Add(new StackedColumnItem
            {
                X0 = categoryIndex - HalfBarWidth,
                X1 = categoryIndex + HalfBarWidth,
                Y0 = y0,
                Y1 = y1,
                Color = FillColor,
                bin = bin,
                total = total,
                group = group,
            });
        }
    }

    /// <summary>
    /// Carries the fields the tracker strings reference. Lower-case names are deliberate - the tracker
    /// format strings bind by property name and were written against HistItem : BarItem.
    /// </summary>
    public class StackedColumnItem : RectangleBarItem
    {
        public int total { get; set; }
        public string bin { get; set; }
        public string group { get; set; }

        /// <summary>This segment's own magnitude, i.e. what BarItem.Value held before the rewrite.</summary>
        public double Value => Y1 - Y0;

        /// <summary>Index on the category axis, i.e. what BarItem.CategoryIndex held.</summary>
        public int CategoryIndex => (int)Math.Round((X0 + X1) / 2);
    }
}
