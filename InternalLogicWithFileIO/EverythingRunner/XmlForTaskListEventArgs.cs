using System;
using System.Collections.Generic;

namespace InternalLogicTaskLayer
{
    public class XmlForTaskListEventArgs : EventArgs
    {

        #region Public Fields

        public List<XmlForTask> newDatabases;

        #endregion Public Fields

        #region Public Constructors

        public XmlForTaskListEventArgs(List<XmlForTask> newDatabases)
        {
            this.newDatabases = newDatabases;
        }

        #endregion Public Constructors

    }
}