using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TaskLayer
{
    public class MyTaskResults
    {
        public List<string> NewSpectra; // calibration writes new calibrated spectra
        public List<DbForTask> NewDatabases; // gptmd writes new annotated databases
        public List<string> NewFileSpecificTomls; // calibration writes suggested ppm tolerances
        public TimeSpan Time;

        private readonly List<string> resultTexts;

        private readonly StringBuilder niceText = new StringBuilder();

        internal MyTaskResults(MetaMorpheusTask s)
        {
            resultTexts = new List<string>();
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("Time to run task: " + Time);
            sb.AppendLine();
            sb.AppendLine();
            sb.AppendLine("--------------------------------------------------");
            if ((NewSpectra != null && NewSpectra.Any()) || (NewDatabases != null && NewDatabases.Any()))
            {
                sb.AppendLine();
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
                sb.AppendLine();
                sb.AppendLine("--------------------------------------------------");
            }
            sb.AppendLine();
            sb.AppendLine();
            sb.AppendLine(niceText.ToString());
            sb.AppendLine();
            sb.AppendLine();
            sb.AppendLine("--------------------------------------------------");
            sb.AppendLine();
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

        internal void AddNiceText(string niceTextString)
        {
            niceText.AppendLine(niceTextString);
        }
    }
}