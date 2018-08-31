using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{

    public class Node
    {
        internal Char value;
        internal Node lChild;
        internal Node rChild;
        internal Node father;

        public Node(Char v, Node l, Node r)
        {
            value = v;
            lChild = l;
            rChild = r;
        }

        public Node(Char v)
        {
            value = v;
            lChild = null;
            rChild = null;
        }

        public override string ToString()
        {
            return value.ToString();
        }
    }
}