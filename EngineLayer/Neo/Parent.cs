namespace EngineLayer.Neo
{
    public class Parent
    {
        public Parent(string id, string seq)
        {
            this.id = id;
            this.seq = seq;
        }

        public string id { get; set; }
        public string seq { get; set; }
    }
}