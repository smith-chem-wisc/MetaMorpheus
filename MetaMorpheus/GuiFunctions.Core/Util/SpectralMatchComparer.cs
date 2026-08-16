using System;
using System.Collections.Generic;
using Readers;

namespace GuiFunctions.Util;

public class SpectralMatchComparer : IEqualityComparer<SpectrumMatchFromTsv>
{
    public static SpectralMatchComparer Instance = new();
    public bool Equals(SpectrumMatchFromTsv x, SpectrumMatchFromTsv y)
    {
        if (ReferenceEquals(x, y)) return true;
        if (x is null) return false;
        if (y is null) return false;
        if (x.GetType() != y.GetType()) return false;
        return x.FullSequence == y.FullSequence && x.Ms2ScanNumber == y.Ms2ScanNumber && x.FileNameWithoutExtension == y.FileNameWithoutExtension && x.PrecursorScanNum == y.PrecursorScanNum;
    }

    public int GetHashCode(SpectrumMatchFromTsv obj)
    {
        return HashCode.Combine(obj.FullSequence, obj.Ms2ScanNumber, obj.FileNameWithoutExtension, obj.PrecursorScanNum);
    }
}
