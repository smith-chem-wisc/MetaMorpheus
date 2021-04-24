using System;
using System.Collections.Generic;
using System.Net;
using System.IO;
using System.Text.RegularExpressions;
using System.Runtime.InteropServices;

namespace EngineLayer
{
    public class DotnetVersion
    {
        public string VersionFetchedFromWeb { get; private set; }

        /// <summary>
        /// Gets the windows dotnet runtime version from a readme text string.
        /// </summary>
        /// <param name="readme"></param>
        /// <returns></returns>
        public static string GetDotnetVersionFromReadme(string readme)
        {
            return Regex.Match(readme, @"runtime-desktop-(\d\.\d\.\d)").Groups[1].ToString();
        }

        /// <summary>
        /// Determines if the present dotnet version is the same as the latest version on  the web based on readmes.
        /// </summary>
        /// <param name="readmeUrl"></param>
        /// <returns></returns>
        public bool IsSameAsLatestWebVersion(string readmeUrl)
        {
            return IsSameAsVersion(VersionFetchedFromWeb ?? GetDotnetVersionFromWeb(readmeUrl));
        }

        public bool VersionOnWebIsAvailble(string readmeUrl)
        {
            var versionPieces = GetDotnetVersionFromWeb(readmeUrl).Split('.');
            string majorVersion = $"{versionPieces[0]}.{versionPieces[1]}";
            // get the list of runtimes available on machine
            // check that the major version is in that list
            // check that runtime and desktop are availble (possibly with Process.Start("dotnet --list-runtimes") and capturing output
            throw new NotImplementedException();
        }

        /// <summary>
        /// Determines if the present dotnet version is the same as a version fetched from a readme.
        /// </summary>
        /// <param name="readmeVersion"></param>
        /// <returns></returns>
        public static bool IsSameAsVersion(string readmeVersion)
        {
            var versionPieces = readmeVersion.Split('.');
            string majorVersion = $"{versionPieces[0]}.{versionPieces[1]}";
            return RuntimeInformation.FrameworkDescription.Contains(majorVersion);
        }

        /// <summary>
        /// Gets the dotnet version from a MetaMorpheus readme URL.
        /// </summary>
        /// <param name="readmeUrl"></param>
        /// <returns></returns>
        private string GetDotnetVersionFromWeb(string readmeUrl)
        {
            using var client = new WebClient();
            var readme = new StreamReader(client.OpenRead(readmeUrl)).ReadToEnd();
            VersionFetchedFromWeb = GetDotnetVersionFromReadme(readme);
            return VersionFetchedFromWeb;
        }
    }
}
