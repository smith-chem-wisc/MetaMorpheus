#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.TypeConversion;

namespace TaskLayer.ParallelSearch.Util.Converters;

public class SemiColonDelimitedToDoubleArrayTypeConverter : DefaultTypeConverter
{
    public static int RoundingPlaces = 4;
    public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
    {
        var splits = text.Split(';');
        var toReturn = splits.Where(p => p != "");
        return toReturn.Select(p => double.Parse(p, CultureInfo.InvariantCulture)).ToArray();
    }

    public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
    {
        if (value is not IEnumerable<double> list)
            return string.Empty;
        return string.Join(';', list.Select(p => Math.Round(p, RoundingPlaces).ToString(CultureInfo.InvariantCulture)));
    }
}