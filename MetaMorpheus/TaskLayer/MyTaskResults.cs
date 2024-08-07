﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Easy.Common.Extensions;

namespace TaskLayer
{
    public class MyTaskResults
    {
        public List<string> NewSpectra; // calibration writes new calibrated spectra
        public List<DbForTask> NewDatabases; // gptmd writes new annotated databases
        public List<string> NewFileSpecificTomls; // calibration writes suggested ppm tolerances
        public TimeSpan Time;

        private readonly List<string> resultTexts;

        private readonly StringBuilder TaskSummaryText = new StringBuilder();
        private readonly StringBuilder PsmPeptideProteinSummaryText = new StringBuilder();

        internal MyTaskResults(MetaMorpheusTask s)
        {
            resultTexts = new List<string>();
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("Time to run task: " + Time);
            sb.AppendLine();
            sb.AppendLine("--------------------------------------------------");
            if ((NewSpectra != null && NewSpectra.Any()) || (NewDatabases != null && NewDatabases.Any()))
            {
                sb.AppendLine();
                sb.AppendLine("New files:");
                if (NewSpectra != null && NewSpectra.Any())
                {
                    sb.AppendLine("New spectra: ");
                    sb.AppendLine();
                    sb.AppendLine(string.Join(Environment.NewLine + "\t", NewSpectra));
                }
                if (NewDatabases != null && NewDatabases.Any())
                {
                    sb.AppendLine("New databases: ");
                    sb.AppendLine(string.Join(Environment.NewLine + "\t", NewDatabases.Select(b => b.FilePath)).ToString());
                }
                sb.AppendLine();
                sb.AppendLine("--------------------------------------------------");
            }
            sb.AppendLine();
            sb.AppendLine(PsmPeptideProteinSummaryText.ToString());
            sb.AppendLine(TaskSummaryText.ToString());
            sb.AppendLine();
            sb.AppendLine("--------------------------------------------------");
            sb.AppendLine();
            sb.AppendLine("Engine Results:");
            sb.AppendLine();
            foreach (var ok in resultTexts)
            {
                sb.AppendLine(ok);
                sb.AppendLine();
            }
            return sb.ToString();
        }

        internal void AddResultText(string resultsText)
        {
            resultTexts.Add(resultsText);
        }

        internal void AddPsmPeptideProteinSummaryText(string targetTextString)
        {
            PsmPeptideProteinSummaryText.Append(targetTextString);
        }

        internal void AddTaskSummaryText(string niceTextString)
        {
            TaskSummaryText.AppendLine(niceTextString);
        }
    }
}