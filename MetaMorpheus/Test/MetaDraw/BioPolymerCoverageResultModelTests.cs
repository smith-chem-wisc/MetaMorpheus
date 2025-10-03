using System.Collections.Generic;
using NUnit.Framework;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using Readers;

namespace Test.MetaDraw;

public class DummySpectralmatch : SpectrumMatchFromTsv
{
    public DummySpectralmatch(double qValue = 0.05, string ambiguity = "3", string filePath = "file1.psmtsv",  string accession = "ACC", string baseSeq = "ABC", string missedCleavages = "0", string startAndEnd = "[1 to 3]") 
        : base("20100723_Velos1_TaGe_SA_Gamg_2-calib-averaged\t37933\t165.90394\t200.00000\t505630.47644\t37925\t3.00000\t151816.23307\t1194.56112\t3580.66153\t44.328\t39.328\t0\tGGGGGGGGGGGLGGGLGNVLGGLISGAGGGGGGGGGGGGGGGGGGGGTAMR\tGGGGGGGGGGGLGGGLGNVLGGLISGAGGGGGGGGGGGGGGGGGGGGTAMR\tGGGGGGGGGGGLGGGLGNVLGGLISGAGGGGGGGGGGGGGGGGGGGGTAMR\t1\t5\t\t\t\t0\t0\t3580.65338\t0.00815\t2.28000\tP04632\tCalpain small subunit 1\tprimary:CAPNS1, synonym:CAPN4, synonym:CAPNS\tHomo sapiens\t\t\tN\tN\tfull\t[10 to 60]\tK\tI\t \tT\t[y1+1, y2+1, y3+1, y4+1, y5+1, y6+1, y7+1, y8+1, y9+1, y10+1, y11+1, y12+1, y13+1, y14+1, y15+1, y16+1, y17+1, y18+1, y19+1, y20+1, y21+1, y22+1, y23+1, y24+1, y25+1, y26+1, y27+1, y28+1];[b2+1, b3+1, b4+1, b5+1, b6+1, b7+1, b8+1, b9+1, b10+1, b11+1, b12+1, b13+1, b14+1, b15+1, b16+1, b18+1]\t[y1+1:175.11889, y2+1:306.15954, y3+1:377.19654, y4+1:478.24453, y5+1:535.26597, y6+1:592.28604, y7+1:649.30788, y8+1:706.32920, y9+1:763.35010, y10+1:820.37208, y11+1:877.39447, y12+1:934.41413, y13+1:991.43748, y14+1:1048.45701, y15+1:1105.47696, y16+1:1162.50121, y17+1:1219.52228, y18+1:1276.53952, y19+1:1333.56720, y20+1:1390.58477, y21+1:1447.60893, y22+1:1504.63266, y23+1:1561.64614, y24+1:1618.67260, y25+1:1689.70263, y26+1:1746.72044, y27+1:1833.76364, y28+1:1946.83155];[b2+1:115.04997, b3+1:172.07162, b4+1:229.09313, b5+1:286.11466, b6+1:343.13611, b7+1:400.15770, b8+1:457.17939, b9+1:514.20071, b10+1:571.22127, b11+1:628.24464, b12+1:741.32787, b13+1:798.34733, b14+1:855.36264, b15+1:912.38839, b16+1:1025.47895, b18+1:1196.54361]\t[y1+1:-0.00006, y2+1:0.00010, y3+1:-0.00001, y4+1:0.00030, y5+1:0.00028, y6+1:-0.00112, y7+1:-0.00074, y8+1:-0.00088, y9+1:-0.00145, y10+1:-0.00093, y11+1:-0.00000, y12+1:-0.00181, y13+1:0.00008, y14+1:-0.00186, y15+1:-0.00337, y16+1:-0.00058, y17+1:-0.00097, y18+1:-0.00520, y19+1:0.00101, y20+1:-0.00288, y21+1:-0.00018, y22+1:0.00208, y23+1:-0.00590, y24+1:-0.00090, y25+1:-0.00799, y26+1:-0.01164, y27+1:-0.00047, y28+1:-0.01663];[b2+1:-0.00023, b3+1:-0.00004, b4+1:-0.00000, b5+1:0.00007, b6+1:0.00005, b7+1:0.00018, b8+1:0.00040, b9+1:0.00026, b10+1:-0.00064, b11+1:0.00126, b12+1:0.00043, b13+1:-0.00157, b14+1:-0.00773, b15+1:-0.00344, b16+1:0.00305, b18+1:0.00333]\t[y1+1:-0.37, y2+1:0.33, y3+1:-0.02, y4+1:0.62, y5+1:0.52, y6+1:-1.89, y7+1:-1.14, y8+1:-1.25, y9+1:-1.90, y10+1:-1.14, y11+1:-0.00, y12+1:-1.94, y13+1:0.08, y14+1:-1.78, y15+1:-3.05, y16+1:-0.50, y17+1:-0.80, y18+1:-4.08, y19+1:0.76, y20+1:-2.07, y21+1:-0.12, y22+1:1.38, y23+1:-3.78, y24+1:-0.56, y25+1:-4.73, y26+1:-6.67, y27+1:-0.26, y28+1:-8.54];[b2+1:-2.03, b3+1:-0.26, b4+1:-0.02, b5+1:0.24, b6+1:0.14, b7+1:0.45, b8+1:0.89, b9+1:0.51, b10+1:-1.13, b11+1:2.01, b12+1:0.58, b13+1:-1.97, b14+1:-9.04, b15+1:-3.78, b16+1:2.98, b18+1:2.78]\t[y1+1:9775, y2+1:1858, y3+1:3049, y4+1:1076, y5+1:1476, y6+1:2977, y7+1:3453, y8+1:3289, y9+1:3932, y10+1:4650, y11+1:4428, y12+1:4487, y13+1:4638, y14+1:4270, y15+1:4701, y16+1:5381, y17+1:3626, y18+1:3552, y19+1:2672, y20+1:2941, y21+1:2807, y22+1:2851, y23+1:2541, y24+1:6153, y25+1:1475, y26+1:2602, y27+1:5708, y28+1:888];[b2+1:2509, b3+1:8489, b4+1:9065, b5+1:8360, b6+1:7417, b7+1:4225, b8+1:4917, b9+1:3361, b10+1:2780, b11+1:3262, b12+1:2070, b13+1:1145, b14+1:908, b15+1:3497, b16+1:932, b18+1:1424]\t44\t0.5968\t \t \t1\t0\t0.000019\t1\t0\t0.000019\t0.00015681982040405273\t5.9E-05\r", ['\t'], HeaderToIndex)
    {
        QValue = qValue;
        AmbiguityLevel = ambiguity;
        FileNameWithoutExtension = filePath;
        Accession = accession;
        BaseSeq = baseSeq;
        MissedCleavage = missedCleavages;
        StartAndEndResiduesInParentSequence = startAndEnd;
    }

    // Header to ensure constructor runs. 
    public static Dictionary<string, int> HeaderToIndex = new Dictionary<string, int>
    {
        { "File Name", 0 },
        { "Scan Number", 1 },
        { "Scan Retention Time", 2 },
        { "Num Experimental Peaks", 3 },
        { "Total Ion Current", 4 },
        { "Precursor Scan Number", 5 },
        { "Precursor Charge", 6 },
        { "Precursor Intensity", 7 },
        { "Precursor MZ", 8 },
        { "Precursor Mass", 9 },
        { "Score", 10 },
        { "Delta Score", 11 },
        { "Notch", 12 },
        { "Base Sequence", 13 },
        { "Full Sequence", 14 },
        { "Essential Sequence", 15 },
        { "Ambiguity Level", 16 },
        { "PSM Count (unambiguous, <0.01 q-value)", 17 },
        { "Mods", 18 },
        { "Mods Chemical Formulas", 19 },
        { "Mods Combined Chemical Formula", 20 },
        { "Num Variable Mods", 21 },
        { "Missed Cleavages", 22 },
        { "Monoisotopic Mass", 23 },
        { "Mass Diff (Da)", 24 },
        { "Mass Diff (ppm)", 25 },
        { "Accession", 26 },
        { "Name", 27 },
        { "Gene Name", 28 },
        { "Organism Name", 29 },
        { "Identified Sequence Variations", 30 },
        { "Splice Sites", 31 },
        { "Contaminant", 32 },
        { "Decoy", 33 },
        { "Description", 34 },
        { "Start and End Residues In Full Sequence", 35 },
        { "Previous Residue", 36 },
        { "Next Residue", 37 },
        { "Theoreticals Searched", 38 },
        { "Decoy/Contaminant/Target", 39 },
        { "Matched Ion Series", 40 },
        { "Matched Ion Mass-To-Charge Ratios", 41 },
        { "Matched Ion Mass Diff (Da)", 42 },
        { "Matched Ion Mass Diff (Ppm)", 43 },
        { "Matched Ion Intensities", 44 },
        { "Matched Ion Counts", 45 },
        { "Normalized Spectral Angle", 46 },
        { "Localized Scores", 47 },
        { "Improvement Possible", 48 },
        { "Cumulative Target", 49 },
        { "Cumulative Decoy", 50 },
        { "QValue", 51 },
        { "Cumulative Target Notch", 52 },
        { "Cumulative Decoy Notch", 53 },
        { "QValue Notch", 54 },
        { "PEP", 55 },
        { "PEP_QValue", 56 },
        { "Intersecting Sequence Variations", -1 },
        { "1/K0", -1 },
        { "Intersectings", -1 },

    };

    // Dont use these, they will crash. 
    public DummySpectralmatch(string line, char[] split, Dictionary<string, int> parsedHeader) : base(line, split, parsedHeader)
    { 
    }

    public DummySpectralmatch(SpectrumMatchFromTsv psm, string fullSequence, int index = 0, string baseSequence = "") : base(psm, fullSequence, index, baseSequence)
    {
    }
}

[TestFixture]
public class BioPolymerCoverageResultModelTests
{
    [Test]
    public void Constructor_SetsProperties()
    {
        var match = new DummySpectralmatch();
        var model = new BioPolymerCoverageResultModel(match, "SEQ", 2, 5, BioPolymerCoverageType.UniqueMissedCleavage);

        Assert.That(model.Match, Is.EqualTo(match));
        Assert.That(model.BaseSequence, Is.EqualTo("SEQ"));
        Assert.That(model.Start, Is.EqualTo(2));
        Assert.That(model.End, Is.EqualTo(5));
        Assert.That(model.CoverageType, Is.EqualTo(BioPolymerCoverageType.UniqueMissedCleavage));
    }

    [Test]
    public void Properties_DelegateToMatch()
    {
        var match = new DummySpectralmatch();
        var model = new BioPolymerCoverageResultModel(match, "SEQ", 1, 3, BioPolymerCoverageType.Shared);

        Assert.That(model.QValue, Is.EqualTo(0.05));
        Assert.That(model.Ambiguity, Is.EqualTo("3"));
        Assert.That(model.FileName, Is.EqualTo("file1.psmtsv"));
    }
}