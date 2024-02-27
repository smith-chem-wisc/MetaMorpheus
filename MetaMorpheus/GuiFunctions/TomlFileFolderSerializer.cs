using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Xml.Linq;
using Newtonsoft.Json;
using Proteomics;
using Nett;
using TaskLayer;
using EngineLayer;

namespace GuiFunctions

{
    /// <summary>
    /// This class is about the implementation of inputs and outputs. 
    /// </summary>
    /// 
    public static class TomlFileFolderSerializer
    {
        public static string PathToCheck = Path.Combine(GlobalVariables.DataDir, "DefaultParameters");
        /// <summary>
        /// Create a file that contains the serialized object. 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="path">The path the file will be created</param>
        /// <param name="objectToSerialize">The object to serialize</param>
        public static void Serialize<T>(string path, T objectToSerialize)
        {
            Toml.WriteFile<T>(objectToSerialize, path, MetaMorpheusTask.tomlConfig);
        }

        /// <summary>
        /// Return a toml deserialized object
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="path">The path of file that will be deserialized</param>
        /// <returns>The deserialized object</returns>
        /// <exception cref="ArgumentNullException"></exception>
        public static T Deserialize<T>(string path)
        {
            T toml = Toml.ReadFile<T>(path, MetaMorpheusTask.tomlConfig);
            return toml;
        }

        /// <summary>
        /// Check if the type exist
        /// </summary>
        /// <param name="type">The type to check</param>
        /// <returns>True if it exist, false if it does not</returns>
        public static bool CheckforTypeFolder(Type type)
        {
            string path = Path.Combine(PathToCheck, type.ToString());
            return Directory.Exists(path);
        }

        /// <summary>
        /// Return the path of the folder (type)
        /// </summary>
        /// <param name="folder">The type to check</param>
        /// <param name="name">The name of the type</param>
        /// <returns></returns>
        public static string GetFilePath(Type folder, string name)
        {
            string txtFilePath = Path.Combine(folder.FullName, $"{name}.toml");
            string path = Path.Combine(PathToCheck, txtFilePath);
            return path;
        }

        /// <summary>
        /// Create the type's folder if it does not exist, save it if it exists. 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="name">name of the object to save</param>
        /// <param name="o">object to save</param>
        public static void Save<T>(string name, T o)
        {
            if (!CheckforTypeFolder(o.GetType()))
            {
                string pathOfFolder = Path.Combine(PathToCheck, o.GetType().ToString());
                Directory.CreateDirectory(pathOfFolder);
            }
            string path = GetFilePath(o.GetType(), name);
            Serialize<T>(path, o);
        }

        /// <summary>
        /// Return an Dictionary of all deserialized file
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="folder">The folder to check</param>
        /// <returns>Dictionary of all deserialized file</returns>
        public static Dictionary<string, T> LoadAllOfTypeT<T>(T obj)
        {
            Dictionary<string, T> dictionary = new Dictionary<string, T>();
            Type folder = obj.GetType();
            if (CheckforTypeFolder(folder))
            {
                string pathOfFolder = Path.Combine(PathToCheck, folder.ToString());
                string[] contentsInFolder = Directory.GetFiles(pathOfFolder);
                foreach (string path in contentsInFolder)
                {
                   dictionary.Add(Path.GetFileNameWithoutExtension(path), Deserialize<T>(path));
                }                
            }
            return dictionary;
        }


        /// <summary>
        /// Delete file
        /// </summary>
        /// <typeparam name="T"></typeparam>m 
        /// <param name="type">type of file to delete</param>
        /// <param name="name">name of the file</param>
        public static void Delete(Type type, string name)
        {
            if (CheckforTypeFolder(type))
            {
                string filePath = GetFilePath(type, name);
                if (File.Exists(filePath))
                {
                    File.Delete(filePath);
                }
            }
        }

        /// <summary>
        /// Returns the names of each saved object of Type passed
        /// </summary>
        /// <param name="type">Type of object to load names for</param>
        /// <returns></returns>
        public static string[] LoadAllNamesOfType(Type type)
        {
            string path = Path.Combine(PathToCheck, type.ToString());
            var folders = Directory.GetFiles(path)
                .Select(p => Path.GetFileNameWithoutExtension(p)).ToArray();
            return folders;
        }     
        
       
    }
}