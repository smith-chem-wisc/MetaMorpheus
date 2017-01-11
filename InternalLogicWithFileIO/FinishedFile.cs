namespace InternalLogicWithFileIO
{
    public class FinishedFile
    {
        public string filepath { get; private set; }

        public FinishedFile(string filepath)
        {
            this.filepath = filepath;
        }
    }
}