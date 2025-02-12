using System.IO;
using TaskLayer;

namespace EngineLayer
{
    public static class Logger
    {
        static bool isAttached;
        static StreamWriter logWriter;

        public static void Attach(string path)
        {
            logWriter = new StreamWriter(path);
            isAttached = true;
        }

        internal static void Log(object sender, StringEventArgs s)
        {
            if (!isAttached) return;

            logWriter.WriteLine(s.S);
            logWriter.Flush(); // Ensure the log is written to the file
        }

        internal static void StartingTaskHander(object sender, SingleTaskEventArgs s)
        {
            if (!isAttached) return;

            logWriter.WriteLine();
            logWriter.WriteLine($"=======================   {s.DisplayName}    =======================");
            logWriter.WriteLine();
            logWriter.Flush(); // Ensure the log is written to the file
        }

        internal static void Detach()
        {
            if (!isAttached) return;
            logWriter.Close();
            isAttached = false;
        }
    }
}
