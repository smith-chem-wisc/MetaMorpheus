using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public static class ByteSerializer
    {
        #region Single Object
        public static byte[] ObjectToByteArray(object obj)
        {
            if (obj == null)
                return null;
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            using (System.IO.MemoryStream ms = new System.IO.MemoryStream())
            {
                bf.Serialize(ms, obj);
                return ms.ToArray();
            }
        }

        public static object ByteArrayToObject(byte[] arrBytes)
        {
            if (arrBytes == null)
                return null;
            using (System.IO.MemoryStream memStream = new System.IO.MemoryStream())
            {
                System.Runtime.Serialization.Formatters.Binary.BinaryFormatter binForm = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
                memStream.Write(arrBytes, 0, arrBytes.Length);
                memStream.Seek(0, System.IO.SeekOrigin.Begin);
                object obj = binForm.Deserialize(memStream);
                return obj;
            }
        }

        #endregion

        #region SingleObjectToFile

        public static void ObjectToByteArrayFile<T>(this T obj, string fileName)
        {
            if (obj == null)
                return;
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            using (System.IO.FileStream fs = new System.IO.FileStream(fileName, System.IO.FileMode.Create))
            {
                bf.Serialize(fs, obj);
            }
        }

        public static T ByteArrayFileToObject<T>(string fileName)
        {
            if (fileName == null)
                return default;
            using (System.IO.FileStream fs = new System.IO.FileStream(fileName, System.IO.FileMode.Open))
            {
                System.Runtime.Serialization.Formatters.Binary.BinaryFormatter binForm = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
                T obj = (T)binForm.Deserialize(fs);
                return obj;
            }
        }

        #endregion

        #region Object Collection in Memory

        public static byte[] ObjectArrayToByteArray(object[] obj)
        {
            if (obj == null)
                return null;
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            using (System.IO.MemoryStream ms = new System.IO.MemoryStream())
            {
                bf.Serialize(ms, obj);
                return ms.ToArray();
            }
        }

        public static object[] ByteArrayToObjectArray(byte[] arrBytes)
        {
            if (arrBytes == null)
                return null;
            using (System.IO.MemoryStream memStream = new System.IO.MemoryStream())
            {
                System.Runtime.Serialization.Formatters.Binary.BinaryFormatter binForm = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
                memStream.Write(arrBytes, 0, arrBytes.Length);
                memStream.Seek(0, System.IO.SeekOrigin.Begin);
                object[] obj = (object[])binForm.Deserialize(memStream);
                return obj;
            }
        }

        #endregion

        #region Object Collection To File

        public static void ObjectArrayToByteArrayFile<T>(this T[] obj, string fileName)
        {
            if (obj == null)
                return;
            System.Runtime.Serialization.Formatters.Binary.BinaryFormatter bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            using (System.IO.FileStream fs = new System.IO.FileStream(fileName, System.IO.FileMode.Create))
            {
                bf.Serialize(fs, obj);
            }
        }

        public static T[] ByteArrayFileToObjectArray<T>(string fileName)
        {
            if (fileName == null)
                return null;
            using (System.IO.FileStream fs = new System.IO.FileStream(fileName, System.IO.FileMode.Open))
            {
                System.Runtime.Serialization.Formatters.Binary.BinaryFormatter binForm = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
                T[] obj = (T[])binForm.Deserialize(fs);
                return obj;
            }
        }

        #endregion
    }
}
