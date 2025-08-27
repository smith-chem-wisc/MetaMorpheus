namespace EngineLayer;

public class AllowedIntervalWithNotch(double minimum, double maximum, double notch)
{
    public readonly double Minimum = minimum;
    public readonly double Maximum = maximum;
    public readonly double Notch = notch;

    public bool Contains(double value)
    {
        return value >= Minimum && value <= Maximum;
    }
}