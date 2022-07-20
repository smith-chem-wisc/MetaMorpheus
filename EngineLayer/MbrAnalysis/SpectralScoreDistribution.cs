using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;
using System.IO;

namespace EngineLayer.MbrAnalysis
{
    public abstract class SpectralScoreDistribution
    {
        public List<SpectralComparison> SpectralComparisonList { get; private set; }
        public string OutputFolder { get; set; }
        public string DistributionType { get; private set; }

        protected SpectralScoreDistribution(string outputFolder)
        {
            OutputFolder = outputFolder; 
            System.IO.Directory.CreateDirectory(outputFolder);
            SpectralComparisonList = new();
        }

        public void AddComparison(string spectraSource, string librarySpectraSource, string fullSequence, SpectralSimilarity similarity)
        {
            SpectralComparisonList.Add(new SpectralComparison(spectraSource, librarySpectraSource, fullSequence, similarity));
        }

        public void WriteResults()
        {
            string fileName = DistributionType + "Distribution.tsv";
            string outputFile = Path.Combine(OutputFolder, fileName);
            using (var output = new StreamWriter(outputFile))
            {
                output.WriteLine(SpectralScoreDistribution.ScoreDistributionHeader);
                foreach (SpectralComparison sc in SpectralComparisonList)
                {
                    output.WriteLine(sc.ToString());
                }
            }

        }

        public static string ScoreDistributionHeader 
        {
            get
            {
                StringBuilder sb = new();
                sb.Append("Spectra File");
                sb.Append('\t');
                sb.Append("Library Spectra File");
                sb.Append('\t');
                sb.Append("Full Sequence");
                sb.Append('\t');
                // Cosine Angle and Spectral Contrast are hard coded and they shouldn't be. This should be extensible any similarity measurement
                sb.Append("Cosine Angle");
                sb.Append('\t');
                sb.Append("Spectral Contrast");
                return sb.ToString();
            }

        }

    }
}
