namespace EngineLayer.Neo
{
    public class Parent
    {
        #region Public Constructors

        public Parent(string id, string seq)
        {
            this.id = id;
            this.seq = seq;
        }

        #endregion Public Constructors

        #region Public Properties

        public string id { get; set; }
        public string seq { get; set; }

        #endregion Public Properties
    }
}