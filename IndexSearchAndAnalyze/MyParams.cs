namespace IndexSearchAndAnalyze
{
    public abstract class MyParams
    {
        public AllTasksParams allTasksParams;

        public MyParams(AllTasksParams allTasksParams)
        {
            this.allTasksParams = allTasksParams;
        }

        internal abstract void Validate();
    }
}