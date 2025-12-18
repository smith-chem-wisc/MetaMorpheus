using System.Text;
using FlashLFQ;

namespace EngineLayer.Quantification;

public sealed class QuantificationEngineResults(MetaMorpheusEngine s) : MetaMorpheusEngineResults(s)
{
    public FlashLfqResults FlashLfqResults { get; set; }

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());

        // TODO: Idk man, write something useful here for the results.txt. Maybe some summary statistics from FlashLFQ?

        return sb.ToString();
    }
}