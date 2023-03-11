using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using FlashLFQ;
using UsefulProteomicsDatabases.Generated;

namespace TaskLayer
{
    public class ExtendedWriter
    {
        public OrderedDictionary HeaderOrderedDictionary { get; private set; }
        /// <summary>
        /// It's difficult to enumerate through an OrderedDictionary directly. Instead,
        /// keys and values can be stored in arrays
        /// </summary>
        private string[] _headerFields;
        /// <summary>
        /// This contains the position of each field in the header.
        /// Fields that occur in the base TabSeparatedHeader have positive values, corresponding
        /// to their position in the original TabSeparatedHeader.Additional fields have position -1
        /// </summary>
        private int[] _headerPositions;
        public string ExtendedHeader => String.Join('\t', _headerFields).Trim();

        public HashSet<string> AdditionalFields;
        public HashSet<string> OriginalFields;
        public Dictionary<Object, Dictionary<string, string>> AdditionalInfoDictionary { get; }
        private Dictionary<string, string> _fieldInfoDictionaryReference = new();

        /// <summary>
        /// ExtendedWriter is a simple way to extend the results of any object that is written to a tsv.
        /// Pass in a list of fields and their desired column position.
        /// ExtendedWrite will create a modified header (Extended Header)
        /// And produce corresponding string for every itemToWrite via the WriteExtendedString method
        /// </summary>
        /// <param name="exampleItem"> Currently accepts any item, but should be modified to accept some interfaced (e.g. IWriteable) </param>
        /// <param name="tabSeparatedHeader"> The tab separated header associated with the example item </param>
        /// <param name="tabSeparatedHeader"> The names of the additional </param>
        public ExtendedWriter(IEnumerable<Object> itemsToWrite, string tabSeparatedHeader, List<(string, int)> newFieldsWithPositions)
        {
            if (!(itemsToWrite.First() is ChromatographicPeak)) throw new NotImplementedException();

            WriteHeaderDictionary(tabSeparatedHeader);
            WriteAdditionalFields(newFieldsWithPositions);
            WriteHeaderArrays();

            AdditionalInfoDictionary = new Dictionary<Object, Dictionary<string, string>>();
            foreach (var peak in itemsToWrite)
            {
                AdditionalInfoDictionary.TryAdd(peak, new Dictionary<string, string>(_fieldInfoDictionaryReference));
            }
        }

        /// <summary>
        /// Instantiates OriginalFields and HeaderOrderedDictionary,
        /// where the name of each header field are keys and the zero indexed position of each header field are values
        /// </summary>
        /// <param name="tabSeparatedHeader"> The tab separated header corresponding to the ToString method of the object </param>
        private void WriteHeaderDictionary(string tabSeparatedHeader)
        {
            string[] headerFields = tabSeparatedHeader.Split('\t').Where(s => s.IsNotNullOrEmptyOrWhiteSpace()).ToArray();

            HeaderOrderedDictionary = new();
            OriginalFields = new HashSet<string>(headerFields);
            for (int position = 0; position < headerFields.Length; position++)
            {
                HeaderOrderedDictionary.Add(headerFields[position], position);
            }
        }

        /// <summary>
        /// Adds new fields to the HeaderDictionary and AdditionalFields.
        /// Instantiates _fieldInfoDictionaryReference
        /// </summary>
        /// <param name="newFieldsWithPositions"></param>
        private void WriteAdditionalFields(List<(string, int)> newFieldsWithPositions)
        {
            _fieldInfoDictionaryReference = new();
            AdditionalFields = new();
            foreach (var fieldPositionPair in newFieldsWithPositions)
            {
                try
                {
                    AddField(fieldPositionPair.Item1, fieldPositionPair.Item2);
                    _fieldInfoDictionaryReference.Add(fieldPositionPair.Item1, "");
                }
                catch (ArgumentException e) { }
            }
        }

        /// <summary>
        /// Adds new fields to the HeaderOrderedDictionary at the appropriate position, as long as the field name doesn't already exist
        /// </summary>
        private void AddField(string fieldName, int position)
        {
            if (AdditionalFields.Contains(fieldName))
                throw new ArgumentException("A field with this name has already been added");
            if (OriginalFields.Contains(fieldName))
                throw new ArgumentException("A field with this name exists in the original header");

            HeaderOrderedDictionary.Insert(position, fieldName, -1);
            AdditionalFields.Add(fieldName);
        }

        /// <summary>
        /// Instantiates the _headerFields and _headerPositions arrays and populates them
        /// with the contents of the HeaderOrderedDictionary. Must be called AFTER WriteAdditionalFields
        /// </summary>
        private void WriteHeaderArrays()
        {
            _headerFields = new string[HeaderOrderedDictionary.Count];
            _headerPositions = new int[HeaderOrderedDictionary.Count];
            HeaderOrderedDictionary.Keys.CopyTo(_headerFields, 0);
            HeaderOrderedDictionary.Values.CopyTo(_headerPositions, 0);
        }

        /// <summary>
        /// Adds additional information to an object
        /// </summary>
        /// <param name="associatedPeak">The object associated with the additional information</param>
        /// <param name="field">Name of the field the information will be associate with</param>
        /// <param name="info">Additional information</param>
        /// <exception cref="ArgumentException"></exception>
        public void AddInfo(Object associatedPeak, string field, string info)
        {
            if (!AdditionalInfoDictionary.ContainsKey(associatedPeak))
                throw new ArgumentException("Associated peak not found");

            if (!AdditionalInfoDictionary[associatedPeak].ContainsKey(field))
                throw new ArgumentException("Field not found");

            AdditionalInfoDictionary[associatedPeak][field] = info;
        }

        /// <summary>
        /// Writes an object to a tab separated string corresponding to the ExtendedHeader
        /// </summary>
        /// <param name="peak"> Object to be written </param>
        /// <returns> A tab separated string </returns>
        /// <exception cref="ArgumentException"></exception>
        public string WriteExtendedString(Object peak)
        {
            if (!AdditionalInfoDictionary.ContainsKey(peak))
                throw new ArgumentException("Object not found");

            string[] originalInfo = peak.ToString().Split('\t');
            string[] extendedInfo = new string[HeaderOrderedDictionary.Count];
            int positionOffset = 0;
            for (int i = 0; i < extendedInfo.Length; i++)
            {
                if (_headerPositions[i] < 0)
                {
                    positionOffset++;
                    extendedInfo[i] = AdditionalInfoDictionary[peak][_headerFields[i]];
                }
                else
                {
                    extendedInfo[i] = originalInfo[i - positionOffset];
                }
            }

            return String.Join('\t', extendedInfo);
        }

        
    }
}
