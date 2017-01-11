namespace InternalLogic
{
    public abstract class MyParams
    {
        public AllTasksParams allTasksParams;
        protected MyParams(AllTasksParams allTasksParams)
        {
            this.allTasksParams = allTasksParams;
        }
        internal abstract void Validate();
    }
}