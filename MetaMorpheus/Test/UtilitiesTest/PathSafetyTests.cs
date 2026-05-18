using NUnit.Framework;
using System;
using System.IO;
using EngineLayer.Util;
using System.Reflection;

namespace Test.UtilitiesTest
{
    [TestFixture]
    public class PathSafetyTests
    {
        [Test]
        public void AppendsRequiredEnding_WhenMissing()
        {
            string dir = @"C:\work\out";
            string input = Path.Combine(dir, "result");
            string safe = PathSafety.MakeSafeOutputPath(input, ".mzML");

            Assert.That(safe, Does.EndWith(".mzML"));
            Assert.That(Path.GetFileName(safe), Is.EqualTo("result.mzML"));
        }

        [Test]
        public void RespectsEnding_CaseInsensitive()
        {
            string dir = @"C:\work\out";
            string input = Path.Combine(dir, "result.MZml");
            string safe = PathSafety.MakeSafeOutputPath(input, ".mzml");

            NUnit.Framework.Assert.That(safe.ToLowerInvariant(), Does.EndWith(".mzml"));
        }

        [Test]
        public void ReplacesInvalidChars_InFilename()
        {
            string dir = @"C:\work\out";
            string input = Path.Combine(dir, "re:sul*t?.mzML"); // invalid chars
            string safe = PathSafety.MakeSafeOutputPath(input, ".mzML");

            NUnit.Framework.Assert.That(Path.GetFileName(safe), Is.EqualTo("re_sul_t_.mzML"));
        }

        [Test]
        public void AvoidsReservedDeviceNames()
        {
            string dir = @"C:\work\out";
            string input = Path.Combine(dir, "CON.mzML");
            string safe = PathSafety.MakeSafeOutputPath(input, ".mzML");

            NUnit.Framework.Assert.That(Path.GetFileName(safe), Is.EqualTo("_CON.mzML"));
        }

        [Test]
        public void TrimsBaseName_ToRespectFileNameLimit()
        {
            string dir = @"C:\work\out";
            string veryLongBase = new string('a', 300);
            string input = Path.Combine(dir, veryLongBase + ".mzML");

            string safe = PathSafety.MakeSafeOutputPath(input, ".mzML", maxPath: 260, maxFileName: 50);

            NUnit.Framework.Assert.That(Path.GetFileName(safe).Length, Is.LessThanOrEqualTo(50));
            NUnit.Framework.Assert.That(safe, Does.EndWith(".mzML"));
        }

        [Test]
        public void TrimsBaseName_ToRespectTotalPathLimit()
        {
            // Make a long directory + base name that would exceed maxPath unless trimmed
            string longDir = @"C:\verylong\" + new string('d', 120);
            string baseName = new string('b', 120);
            string input = Path.Combine(longDir, baseName) + ".mzML";

            string safe = PathSafety.MakeSafeOutputPath(input, ".mzML", maxPath: 200, maxFileName: 180);

            NUnit.Framework.Assert.That(safe.Length, Is.LessThanOrEqualTo(200));
            NUnit.Framework.Assert.That(safe, Does.EndWith(".mzML"));
            // Base is trimmed, but directory remains intact
            NUnit.Framework.Assert.That(safe.StartsWith(longDir, StringComparison.Ordinal), Is.True);
        }

        [Test]
        public void Throws_WhenDirectoryAloneExceedsPathLimit()
        {
            // Directory is already too long; cannot be fixed by trimming file name.
            string longDir = @"C:\" + new string('x', 300);
            string input = Path.Combine(longDir, "file.mzML");

            NUnit.Framework.Assert.That(
                () => PathSafety.MakeSafeOutputPath(input, ".mzML", maxPath: 260, maxFileName: 255),
                Throws.TypeOf<PathTooLongException>());
        }

        [Test]
        public void SuppliesDefaultFileName_WhenOnlyDirectoryProvided()
        {
            string dirOnly = @"C:\work\onlydir\";
            string safe = PathSafety.MakeSafeOutputPath(dirOnly, ".mzML");

            NUnit.Framework.Assert.That(Path.GetDirectoryName(safe), Is.EqualTo(Path.GetFullPath(dirOnly).TrimEnd(Path.DirectorySeparatorChar)));
            NUnit.Framework.Assert.That(Path.GetFileName(safe), Is.EqualTo("output.mzML"));
        }

        [Test]
        public void Throws_OnNullOrEmptyInputs()
        {
            NUnit.Framework.Assert.That(() => PathSafety.MakeSafeOutputPath(null, ".mzML"), Throws.TypeOf<ArgumentException>());
            NUnit.Framework.Assert.That(() => PathSafety.MakeSafeOutputPath("C:\\a\\b", ""), Throws.TypeOf<ArgumentException>());
        }

        [Test]
        public void LeavesAlreadyValidPath_Alone()
        {
            string input = @"C:\work\out\good_name.mzML";
            string safe = PathSafety.MakeSafeOutputPath(input, ".mzML");

            NUnit.Framework.Assert.That(safe, Is.EqualTo(input));
        }

        [Test]
        public void Throws_OnWhitespacePath()
        {
            Assert.That(
                () => PathSafety.MakeSafeOutputPath("   \t", ".mzML"),
                Throws.TypeOf<ArgumentException>().With.Message.Contains("Path is null or whitespace"));
        }

        [Test]
        public void Throws_OnNullRequiredEnding()
        {
            Assert.That(
                () => PathSafety.MakeSafeOutputPath(@"C:\out\file", null),
                Throws.TypeOf<ArgumentException>().With.Message.Contains("Required ending is null or whitespace"));
        }

        [Test]
        public void Throws_OnWhitespaceRequiredEnding()
        {
            Assert.That(
                () => PathSafety.MakeSafeOutputPath(@"C:\out\file", "   "),
                Throws.TypeOf<ArgumentException>().With.Message.Contains("Required ending is null or whitespace"));
        }

        [Test]
        public void Throws_WhenDirectoryAlreadyExceedsMaxPath()
        {
            // Construct a directory longer than the maxPath by itself.
            string tooLongDir = @"C:\" + new string('x', 400);
            string input = Path.Combine(tooLongDir, "file.mzML");

            Assert.That(
                () => PathSafety.MakeSafeOutputPath(input, ".mzML", maxPath: 260, maxFileName: 255),
                Throws.TypeOf<PathTooLongException>().With.Message.Contains("Directory portion exceeds"));
        }

        [Test]
        public void Throws_WhenDirectoryOnlyAndTooLong()
        {
            // Ends with a separator -> treated as directory-only.
            string tooLongDir = @"C:\" + new string('x', 400) + Path.DirectorySeparatorChar;
            Assert.That(
                () => PathSafety.MakeSafeOutputPath(tooLongDir, ".mzML", maxPath: 260, maxFileName: 255),
                Throws.TypeOf<PathTooLongException>().With.Message.Contain("exceeds"));
        }

        [Test]
        public void Throws_WhenEvenMinBasePlusEndingExceedsMaxPath()
        {
            // Make dir such that: dirWithSep + (1-char base) + ".mzML" > maxPath
            string dir = @"C:\verylong\" + new string('d', 40);
            string dirWithSep = dir.EndsWith(Path.DirectorySeparatorChar.ToString(), StringComparison.Ordinal)
                ? dir
                : dir + Path.DirectorySeparatorChar;

            string ending = ".mzML";
            int minPossibleFileLen = 1 + ending.Length; // 1-char base + ending

            // Choose maxPath just below the minimum required to fit (force failure by 1 char)
            int maxPath = dirWithSep.Length + minPossibleFileLen - 1;

            string input = Path.Combine(dir, new string('b', 100)) + ending;

            Assert.That(
                () => PathSafety.MakeSafeOutputPath(input, ending, maxPath: maxPath, maxFileName: 255),
                Throws.TypeOf<PathTooLongException>().With.Message.Contains("exceeds maximum path length"));
        }

        [Test]
        public void TrimsBaseName_ToTinyFileNameLimit_ButThrowsIfTotalPathStillTooLong()
        {
            // Directory part leaves very limited room for filename within maxPath.
            string dir = @"C:\fixed\" + new string('q', 120);
            string dirWithSep = dir.EndsWith(Path.DirectorySeparatorChar.ToString(), StringComparison.Ordinal)
                ? dir
                : dir + Path.DirectorySeparatorChar;

            string ending = ".mzML";
            string longBase = new string('z', 300);
            string input = Path.Combine(dir, longBase) + ending;

            // Make filename limit tiny; method will try to keep at least 1 char of base.
            int maxFileName = 6; // smaller than base+ending; 1-char base + 5-char ".mzML" would just fit
                                 // But set total path to be too small even for that 1+ending after directory
            int maxPath = dirWithSep.Length + maxFileName - 1; // 1 less than needed overall

            Assert.That(
                () => PathSafety.MakeSafeOutputPath(input, ending, maxPath: maxPath, maxFileName: maxFileName),
                Throws.TypeOf<PathTooLongException>());
        }

        [Test]
        public void Handles_PathEndingWithDirectorySeparator_ButThrowsIfTooLong()
        {
            string longDir = @"C:\" + new string('a', 300) + Path.DirectorySeparatorChar;
            // With no filename, method will synthesize "output.mzML"; directory alone is already too long.
            Assert.That(
                () => PathSafety.MakeSafeOutputPath(longDir, ".mzML", maxPath: 260, maxFileName: 255),
                Throws.TypeOf<PathTooLongException>());
        }

        [Test]
        public void Throws_WhenMaxPathIsSmallerThanJustEndingAndOneCharFile()
        {
            // Minimal viable file length = 1 + ending length
            string dir = @"C:\tiny";
            string ending = ".mzML";
            string dirWithSep = dir.EndsWith(Path.DirectorySeparatorChar.ToString(), StringComparison.Ordinal)
                ? dir
                : dir + Path.DirectorySeparatorChar;

            // Set maxPath to be smaller than dirWithSep + (1+ending)
            int maxPath = dirWithSep.Length + (1 + ending.Length) - 1;

            string input = Path.Combine(dir, "base") + ending;

            Assert.That(
                () => PathSafety.MakeSafeOutputPath(input, ending, maxPath: maxPath, maxFileName: 255),
                Throws.TypeOf<PathTooLongException>());
        }

        [Test]
        public void EnsureDirWithSeparator_ReturnsEmptyString_WhenDirectoryIsNull()
        {
            var method = typeof(PathSafety).GetMethod("EnsureDirWithSeparator", BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(method, Is.Not.Null);

            var result = (string)method.Invoke(null, new object[] { null });
            Assert.That(result, Is.EqualTo(string.Empty));
        }

        [Test]
        public void EnsureDirWithSeparator_ReturnsEmptyString_WhenDirectoryIsEmpty()
        {
            var method = typeof(PathSafety).GetMethod("EnsureDirWithSeparator", BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(method, Is.Not.Null);

            var result = (string)method.Invoke(null, new object[] { string.Empty });
            Assert.That(result, Is.EqualTo(string.Empty));
        }

        [Test]
        public void CombineDirectoryAndFile_ReturnsFileName_WhenDirectoryIsNull()
        {
            var method = typeof(PathSafety).GetMethod("CombineDirectoryAndFile", BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(method, Is.Not.Null);

            var result = (string)method.Invoke(null, new object[] { null, "file.txt" });
            Assert.That(result, Is.EqualTo("file.txt"));
        }

        [Test]
        public void CombineDirectoryAndFile_ReturnsFileName_WhenDirectoryIsEmpty()
        {
            var method = typeof(PathSafety).GetMethod("CombineDirectoryAndFile", BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(method, Is.Not.Null);

            var result = (string)method.Invoke(null, new object[] { string.Empty, "file.txt" });
            Assert.That(result, Is.EqualTo("file.txt"));
        }

        [Test]
        public void ToSafeOutputPath_LeavesAlreadyValidPath_Alone()
        {
            string input = @"C:\work\out\good_name.mzML"; string safe = input.ToSafeOutputPath(".mzML");
            Assert.That(safe, Is.EqualTo(input));
        }

        // Extension method: throws on whitespace path
        [Test]
        public void ToSafeOutputPath_Throws_OnWhitespacePath() 
        { 
            Assert.That( () => "   \t".ToSafeOutputPath(".mzML"), Throws.TypeOf<ArgumentException>().With.Message.Contains("Path is null or whitespace")); 
        }

    }
}