
namespace EngineLayer
{

    public class Node
    {
        public Node(char v, Node l, Node r, Node m)
        {
            Value = v;
            LeftChild = l;
            RightChild = r;
            MiddleChild = m;
            Level = null;
        }

        public Node(char v)
        {
            Value = v;
            LeftChild = null;
            RightChild = null;
            MiddleChild = null;
            Level = null;
        }

        public Node(char v, int l)
        {
            Value = v;
            LeftChild = null;
            RightChild = null;
            MiddleChild = null;
            Level = l;
        }

        public char Value { get; set; }
        public Node Father { get; set; }
        public Node LeftChild { get; set; }
        public Node RightChild { get; set; }
        public Node MiddleChild { get; set; }
        public int? Level { get; set; }

        public override string ToString()
        {
            return Value.ToString();
        }
    }
}