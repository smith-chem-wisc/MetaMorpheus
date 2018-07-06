using System.Text;

namespace EngineLayer.Neo
{
    public class PsmTsvLine
    {
        public PsmTsvLine(string[] line, int scanNumber, double score, string baseSequence, string fullSequence, string accession, string proteinName, string geneName, string dct, string target, string decoy, string q)
        {
            Line = line;
            ScanNumber = scanNumber;
            Score = score;
            BaseSequence = baseSequence;
            FullSequence = fullSequence;
            Accession = accession;
            ProteinName = proteinName;
            GeneName = geneName;
            DCT = dct;
            Target = target;
            Decoy = decoy;
            Q = q;
            NeoType = null;
        }

        public string[] Line { get; set; }
        public int ScanNumber { get; set; }
        public double Score { get; set; }
        public string BaseSequence { get; set; }
        public string FullSequence { get; set; }
        public string Accession { get; set; }
        public string ProteinName { get; set; }
        public string GeneName { get; set; }
        public string DCT { get; set; }
        public string Target { get; set; }
        public string Decoy { get; set; }
        public string Q { get; set; }
        public NeoType? NeoType { get; set; }

        public PsmTsvLine AggregateLine(PsmTsvLine secondary)
        {
            BaseSequence += "|" + secondary.BaseSequence;
            FullSequence += "|" + secondary.FullSequence;
            Accession += "|" + secondary.Accession;
            ProteinName += "|" + secondary.ProteinName;
            GeneName += "|" + secondary.GeneName;
            return this;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < Line.Length; i++)
            {
                if (i == ImportPsmtsv.BaseIndex)
                    sb.Append(BaseSequence + '\t');
                else if (i == ImportPsmtsv.FullIndex)
                    sb.Append(FullSequence + '\t');
                else if (i == ImportPsmtsv.AccessionIndex)
                    sb.Append(Accession + '\t');
                else if (i == ImportPsmtsv.ProteinIndex)
                    sb.Append(ProteinName + '\t');
                else if (i == ImportPsmtsv.GeneIndex)
                    sb.Append(GeneName + '\t');
                else if (i == ImportPsmtsv.TargetIndex)
                    sb.Append(Target + '\t');
                else if (i == ImportPsmtsv.DecoyIndex)
                    sb.Append(Decoy + '\t');
                else if (i == ImportPsmtsv.QIndex)
                    sb.Append(Q + '\t');
                else
                    sb.Append(Line[i] + '\t');
            }
            if (NeoType != null)
                sb.Append(NeoType.ToString());
            return sb.ToString();
        }
    }
}