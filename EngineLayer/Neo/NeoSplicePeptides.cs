using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Neo
{
    public static class NeoSplicePeptides
    {
        public static void SplicePeptides(string NFilePath, string CFilePath)
        {
            string header = "";
            string[] NLines = (System.IO.File.ReadAllLines(NFilePath));
            string[] CLines = (System.IO.File.ReadAllLines(CFilePath));
        }
    }
}
