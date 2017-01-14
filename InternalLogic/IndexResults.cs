using System.Collections.Generic;
using System.Text;

namespace InternalLogicEngineLayer
{
	public class IndexResults : MyResults
	{

		#region Public Constructors

		public IndexResults(List<CompactPeptide> peptideIndex, Dictionary<float, List<int>> fragmentIndexDict, IndexEngine indexParams) : base(indexParams)
		{
			this.peptideIndex = peptideIndex;
			this.fragmentIndexDict = fragmentIndexDict;
		}

		#endregion Public Constructors

		#region Public Properties

		public Dictionary<float, List<int>> fragmentIndexDict { get; private set; }
		public List<CompactPeptide> peptideIndex { get; private set; }

		#endregion Public Properties

		#region Protected Methods

		protected override string GetStringForOutput()
		{
			var sb = new StringBuilder();
			sb.AppendLine("\t\tfragmentIndexDict.Count: " + fragmentIndexDict.Count);
			sb.Append("\t\tpeptideIndex.Count: " + peptideIndex.Count);
			return sb.ToString();
		}

		#endregion Protected Methods

	}
}