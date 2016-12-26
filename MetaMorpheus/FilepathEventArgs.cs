using System;

namespace MetaMorpheus
{
    public class FilepathEventArgs : EventArgs
    {
        public string Filepath { get; private set; }

        public FilepathEventArgs(string filepath)
        {
            Filepath = filepath;
        }
    }
}