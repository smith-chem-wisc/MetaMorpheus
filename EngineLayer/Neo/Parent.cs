namespace EngineLayer.Neo
{
    public class Parent
    {
        public Parent(string id, string seq)
        {
            ID = id;
            Seq = seq;
        }

        public string ID { get; set; }
        public string Seq { get; set; }
    }
}