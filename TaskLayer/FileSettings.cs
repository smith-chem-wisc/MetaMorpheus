using System;
using EngineLayer;
using Nett;
using System.Linq.Expressions;

namespace TaskLayer

{
    public class FileSettings
    {


        public FileSettings()
        {
            // Set default values here:
            Protease = GlobalEngineLevelSettings.ProteaseDictionary["trypsin"];
        }


        public Protease Protease { get; set; }


    }
}