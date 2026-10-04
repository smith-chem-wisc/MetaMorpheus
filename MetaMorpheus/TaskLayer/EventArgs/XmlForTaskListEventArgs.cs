using System;
using System.Collections.Generic;
using EngineLayer.DatabaseLoading;

namespace TaskLayer
{
    public class XmlForTaskListEventArgs : EventArgs
    {
        public List<DbForTask> NewDatabases;

        public XmlForTaskListEventArgs(List<DbForTask> newDatabases)
        {
            NewDatabases = newDatabases;
        }
    }
}