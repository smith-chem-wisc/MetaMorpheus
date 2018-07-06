namespace TaskLayer
{
    public enum TaskType
    {
        /// <summary>
        /// Search tasks
        /// </summary>
        Search,

        /// <summary>
        /// Global PTM Discovery tasks
        /// </summary>
        Gptmd,

        /// <summary>
        /// Calibration tasks
        /// </summary>
        Calibrate,

        /// <summary>
        /// Crosslinked peptide search tasks
        /// </summary>
        XLSearch,

        /// <summary>
        /// Posttranslational spliced peptide search tasks
        /// </summary>
        Neo
    }
}