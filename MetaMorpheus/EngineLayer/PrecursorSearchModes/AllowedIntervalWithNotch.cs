namespace EngineLayer;

public class AllowedIntervalWithNotch(double minimum, double maximum, int notch)
{
    public readonly double Minimum = minimum;
    public readonly double Maximum = maximum;
    public readonly int Notch = notch;

    public bool Contains(double value)
    {
        return value >= Minimum && value <= Maximum;
    }
}