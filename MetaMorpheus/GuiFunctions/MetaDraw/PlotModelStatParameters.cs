namespace GuiFunctions.MetaDraw
{
    /// <summary>
    /// Model class for Data Visualization plot parameters
    /// </summary>
    public class PlotModelStatParameters
    {
        /// <summary>
        /// When true, uses logarithmic scale for Y-axis in histograms
        /// </summary>
        public bool UseLogScaleYAxis { get; set; } = false;

        /// <summary>
        /// Property to group PSMs by for nested X-axis plotting
        /// Options: "None", "Notch", "Precursor Charge", "File Name", "Ambiguity Level", "Missed Cleavages"
        /// </summary>
        public string GroupingProperty { get; set; } = "None";

        /// <summary>
        /// Minimum percentage of total for data points to be displayed (0-100)
        /// Applied to all histogram plots
        /// </summary>
        public double MinRelativeCutoff { get; set; } = 0.0;

        /// <summary>
        /// Maximum percentage of total for data points to be displayed (0-100)
        /// Applied to all histogram plots
        /// </summary>
        public double MaxRelativeCutoff { get; set; } = 100.0;

        /// <summary>
        /// When true, normalizes histogram counts to file totals
        /// </summary>
        public bool NormalizeHistogramToFile { get; set; } = false;

        /// <summary>
        /// When true, only displays PSMs that pass the current filters
        /// </summary>
        public bool DisplayFilteredOnly { get; set; } = true;

        public PlotModelStatParameters Clone()
        {
            return new PlotModelStatParameters
            {
                UseLogScaleYAxis = this.UseLogScaleYAxis,
                GroupingProperty = this.GroupingProperty,
                MinRelativeCutoff = this.MinRelativeCutoff,
                MaxRelativeCutoff = this.MaxRelativeCutoff,
                NormalizeHistogramToFile = this.NormalizeHistogramToFile,
                DisplayFilteredOnly = this.DisplayFilteredOnly
            };
        }

        public bool Equals(PlotModelStatParameters other)
        {
            if (other == null) return false;
            return UseLogScaleYAxis == other.UseLogScaleYAxis
                   && GroupingProperty == other.GroupingProperty
                   && System.Math.Abs(MinRelativeCutoff - other.MinRelativeCutoff) < 0.0001
                   && System.Math.Abs(MaxRelativeCutoff - other.MaxRelativeCutoff) < 0.0001
                   && NormalizeHistogramToFile == other.NormalizeHistogramToFile
                   && DisplayFilteredOnly == other.DisplayFilteredOnly;
        }
    }
}
