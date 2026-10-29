using System.Globalization;

namespace GuiFunctions.Util;

public static class RetentionTimeRangeParser
{
    public static bool TryParse(string minimumText, string maximumText, out double minimum, out double maximum, out string errorMessage)
    {
        minimum = 0;
        maximum = double.MaxValue;
        errorMessage = null;

        if (!string.IsNullOrWhiteSpace(minimumText)
            && (!double.TryParse(minimumText, NumberStyles.Float, CultureInfo.InvariantCulture, out minimum)
                || !double.IsFinite(minimum) || minimum < 0))
        {
            errorMessage = "Minimum retention time must be a non-negative finite number of minutes.";
            return false;
        }

        if (!string.IsNullOrWhiteSpace(maximumText)
            && (!double.TryParse(maximumText, NumberStyles.Float, CultureInfo.InvariantCulture, out maximum)
                || !double.IsFinite(maximum) || maximum < 0))
        {
            errorMessage = "Maximum retention time must be a non-negative finite number of minutes or blank.";
            return false;
        }

        if (maximum < minimum)
        {
            errorMessage = "Maximum retention time must be greater than or equal to minimum retention time.";
            return false;
        }

        return true;
    }
}
