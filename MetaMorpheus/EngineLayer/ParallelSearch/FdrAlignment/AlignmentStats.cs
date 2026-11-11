#nullable enable

namespace EngineLayer.ParallelSearch.FdrAlignment;

public readonly record struct AlignmentStats(int AlignedCount, int ClampedHighCount, int ClampedLowCount)
{
    public static AlignmentStats Empty => new(0, 0, 0);

    public static AlignmentStats operator +(AlignmentStats left, AlignmentStats right)
    {
        return new AlignmentStats(
            left.AlignedCount + right.AlignedCount,
            left.ClampedHighCount + right.ClampedHighCount,
            left.ClampedLowCount + right.ClampedLowCount);
    }
}
