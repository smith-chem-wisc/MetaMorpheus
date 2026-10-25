using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using Readers;
using System.Reflection.Metadata.Ecma335;

namespace TaskLayer
{
    public class MyFileManager
    {
        private readonly bool DisposeOfFileWhenDone;
        private readonly ConcurrentDictionary<string, MsDataFile> MyMsDataFiles = new ConcurrentDictionary<string, MsDataFile>(); // Dictionary to hold loaded MsDataFiles, with the file path as the key.

        public MyFileManager(bool disposeOfFileWhenDone)
        {
            DisposeOfFileWhenDone = disposeOfFileWhenDone;
        }

        public static event EventHandler<StringEventArgs> WarnHandler;

        public bool SeeIfOpen(string path)
        {
            return MyMsDataFiles.TryGetValue(path, out var file);
        }

        public MsDataFile LoadFile(string origDataFile, CommonParameters commonParameters)
        {
            FilteringParams filter = new FilteringParams(
                commonParameters.NumberOfPeaksToKeepPerWindow,
                commonParameters.MinimumAllowedIntensityRatioToBasePeak,
                commonParameters.WindowWidthThomsons,
                commonParameters.NumberOfWindows,
                commonParameters.NormalizePeaksAccrossAllWindows,
                commonParameters.TrimMs1Peaks,
                commonParameters.TrimMsMsPeaks);

            if (commonParameters.DissociationType == DissociationType.LowCID || commonParameters.MS2ChildScanDissociationType == DissociationType.LowCID || commonParameters.MS3ChildScanDissociationType == DissociationType.LowCID)
            {
                filter = null;
            }

            if (MyMsDataFiles.TryGetValue(origDataFile, out MsDataFile value))
                return value;

            return MyMsDataFiles.GetOrAdd(origDataFile,
                MsDataFileReader.GetDataFile(origDataFile).LoadAllStaticData(filter, commonParameters.MaxThreadsToUsePerFile));
        }

        /// <summary>
        /// Pairs each scan of a loaded file with the same scan in a second read of that file, by position: one reader
        /// returns a file's scans in the same order however it filters their peaks. A scan number alone does not identify
        /// a scan, because an MGF can repeat a SCANS= value and a static load keeps every such block. A position whose two
        /// scans have different scan numbers is left unpaired.
        /// </summary>
        /// <returns>The stored scan for each paired loaded scan, keyed by reference.</returns>
        public static Dictionary<MsDataScan, MsDataScan> PairWithStoredScans(MsDataFile loaded, MsDataFile stored)
        {
            List<MsDataScan> loadedScans = loaded.GetAllScansList();
            List<MsDataScan> storedScans = stored.GetAllScansList();
            var storedScanFor = new Dictionary<MsDataScan, MsDataScan>(ReferenceEqualityComparer.Instance);
            for (int i = 0; i < Math.Min(loadedScans.Count, storedScans.Count); i++)
            {
                if (loadedScans[i].OneBasedScanNumber == storedScans[i].OneBasedScanNumber)
                {
                    storedScanFor.Add(loadedScans[i], storedScans[i]);
                }
            }
            return storedScanFor;
        }

        /// <summary>
        /// Sets each scan's reporter ion intensities from the scan they are read from (<see cref="Ms2ScanWithSpecificMass.ReporterIonScan"/>)
        /// as <paramref name="filePath"/> stores it. Those scans are read one at a time through the reader's dynamic connection,
        /// so no second copy of the whole file is held beside <paramref name="loaded"/>. A file that cannot be read so is read
        /// whole and paired by position (<see cref="PairWithStoredScans"/>): a format other than .raw, .mzML and .mgf; an MGF
        /// that repeats a SCANS= value, which its scan index refuses; and a file whose read by scan number returns another scan
        /// than the loaded one, as an mzML does whose native IDs do not number its spectra from 1.
        /// </summary>
        /// <returns>True when the scans were read one at a time, false when the file was read whole.</returns>
        public static bool SetReporterIonIntensitiesFromStoredScans(string filePath, MsDataFile loaded,
            IReadOnlyCollection<Ms2ScanWithSpecificMass> scans, IsobaricMassTag massTag, int maxThreads)
        {
            if (TrySetReporterIonIntensitiesScanByScan(filePath, scans, massTag))
            {
                return true;
            }

            var storedScanFor = PairWithStoredScans(loaded, MsDataFileReader.GetDataFile(filePath).LoadAllStaticData(null, maxThreads));
            foreach (var scan in scans)
            {
                scan.SetIsobaricMassTagReporterIonIntensities(massTag, storedScanFor);
            }
            return false;
        }

        /// <summary>
        /// Reads each reporter ion scan by its scan number and keeps it only if it is the loaded scan: the same scan number,
        /// MS order, native ID and retention time. At the first one that is not, it stops and returns false, and the caller
        /// reads the file whole, which sets every scan again.
        /// </summary>
        private static bool TrySetReporterIonIntensitiesScanByScan(string filePath, IEnumerable<Ms2ScanWithSpecificMass> scans, IsobaricMassTag massTag)
        {
            MsDataFile stored = MsDataFileReader.GetDataFile(filePath);
            if (stored is not (ThermoRawFileReader or Mzml or Mgf))
            {
                return false;
            }

            try
            {
                stored.InitiateDynamicConnection();
            }
            catch (MzLibException)
            {
                return false; // an MGF that repeats a SCANS= value
            }

            try
            {
                var storedScanFor = new Dictionary<MsDataScan, MsDataScan>(ReferenceEqualityComparer.Instance);
                foreach (var scansOfOneSpectrum in scans.GroupBy(scan => scan.ReporterIonScan, ReferenceEqualityComparer.Instance))
                {
                    MsDataScan loadedScan = (MsDataScan)scansOfOneSpectrum.Key;
                    MsDataScan storedScan = stored.GetOneBasedScanFromDynamicConnection(loadedScan.OneBasedScanNumber);
                    if (storedScan == null || !Identity(storedScan).Equals(Identity(loadedScan)))
                    {
                        return false;
                    }

                    // one spectrum at a time, so only the scan being read is held
                    storedScanFor.Clear();
                    storedScanFor.Add(loadedScan, storedScan);
                    foreach (var scan in scansOfOneSpectrum)
                    {
                        scan.SetIsobaricMassTagReporterIonIntensities(massTag, storedScanFor);
                    }
                }
                return true;
            }
            catch (MzLibException)
            {
                return false; // an MGF asked for a scan number its index does not hold
            }
            finally
            {
                stored.CloseDynamicConnection();
            }
        }

        private static (int ScanNumber, int MsnOrder, string NativeId, double RetentionTime) Identity(MsDataScan scan) =>
            (scan.OneBasedScanNumber, scan.MsnOrder, scan.NativeId, scan.RetentionTime);

        internal void DoneWithFile(string origDataFile)
        {
            if (DisposeOfFileWhenDone)
            {
                // This would only return false if the file was not in the dictionary
                MyMsDataFiles.TryRemove(origDataFile, out var file);
            }
        }

        private void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
        }

        public static void CompressDirectory(DirectoryInfo directorySelected)
        {
            foreach (FileInfo fileToCompress in directorySelected.GetFiles())
            {
                CompressFile(fileToCompress);
            }
        }

        public static void CompressFile(FileInfo fileToCompress)
        {
            using (FileStream originalFileStream = fileToCompress.OpenRead())
            {
                if ((File.GetAttributes(fileToCompress.FullName) &
                   FileAttributes.Hidden) != FileAttributes.Hidden & fileToCompress.Extension != ".gz")
                {
                    using (FileStream compressedFileStream = File.Create(fileToCompress.FullName + ".gz"))
                    {
                        using (GZipStream compressionStream = new GZipStream(compressedFileStream,
                           CompressionMode.Compress))
                        {
                            originalFileStream.CopyTo(compressionStream);
                        }
                    }
                }
            }

            if (File.Exists(fileToCompress.FullName))
            {
                File.Delete(fileToCompress.FullName);
            }
        }

        public static void DecompressDirectory(DirectoryInfo directorySelected)
        {
            foreach (FileInfo fileToDecompress in directorySelected.GetFiles())
            {
                DecompressFile(fileToDecompress);
            }
        }

        public static void DecompressFile(FileInfo fileToDecompress)
        {
            using (FileStream originalFileStream = fileToDecompress.OpenRead())
            {
                string currentFileName = fileToDecompress.FullName;
                string newFileName = currentFileName.Remove(currentFileName.Length - fileToDecompress.Extension.Length);

                using (FileStream decompressedFileStream = File.Create(newFileName))
                {
                    using (GZipStream decompressionStream = new GZipStream(originalFileStream, CompressionMode.Decompress))
                    {
                        decompressionStream.CopyTo(decompressedFileStream);
                        Console.WriteLine("Decompressed: {0}", fileToDecompress.Name);
                    }
                }
            }

            if (File.Exists(fileToDecompress.FullName))
            {
                File.Delete(fileToDecompress.FullName);
            }
        }
    }
}