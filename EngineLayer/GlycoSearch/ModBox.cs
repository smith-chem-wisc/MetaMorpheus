namespace EngineLayer
{
    public class ModBox
    {
        public ModBox(int[] ids)
        {
            ModIds = ids;
            NumberOfMods = ids.Length;
            TargetDecoy = true;
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public double Mass { get; set; }

        public double DecoyMass { get; set; }

        public bool TargetDecoy { get; set; }

    }
}
