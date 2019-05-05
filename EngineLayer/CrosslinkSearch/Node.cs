
namespace EngineLayer
{

    public class Node
    {
        internal char value;
        internal Node lChild;
        internal Node rChild;
        internal Node mChild;
        internal Node father;
        internal int? level;

        public Node(char v, Node l, Node r, Node m)
        {
            value = v;
            lChild = l;
            rChild = r;
            mChild = m;
            level = null;
        }

        public Node(char v)
        {
            value = v;
            lChild = null;
            rChild = null;
            mChild = null;
            level = null;
        }

        public Node(char v, int l)
        {
            value = v;
            lChild = null;
            rChild = null;
            mChild = null;
            level = l;
        }

        public char Value { get { return value; } }
        public Node Father { get { return father; } }
        public Node LeftChild { get { return lChild; } }
        public Node RightChild { get { return rChild; } }
        public Node MiddleChild { get { return mChild; } }
        public int Level { get { return level.Value; } }

        public override string ToString()
        {
            return value.ToString();
        }
    }
}