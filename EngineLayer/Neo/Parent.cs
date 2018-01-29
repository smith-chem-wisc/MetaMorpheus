using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Neo
{
    public class Parent
    {
        public string id { get; set; }
        public string seq { get; set; }

        public Parent(string id, string seq)
        {
            this.id = id;
            this.seq = seq;
        }
    }
}
