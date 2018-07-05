using System;
using System.Collections.Generic;

namespace TaskLayer
{
    public class XmlForTaskListEventArgs : EventArgs
    {
        public List<DbForTask> newDatabases;

        public XmlForTaskListEventArgs(List<DbForTask> newDatabases)
        {
            this.newDatabases = newDatabases;
        }
    }
}