using System.Collections.Generic;
using System.Text;

namespace InternalLogicEngineLayer
{
	public class ModernSearchResults : MyResults
	{

		#region Public Constructors

		public ModernSearchResults(List<ModernSpectrumMatch>[] newPsms, ModernSearchEngine s) : base(s)
		{
			this.newPsms = newPsms;
		}

		#endregion Public Constructors

		#region Public Properties

		public List<ModernSpectrumMatch>[] newPsms { get; private set; }

		#endregion Public Properties

		#region Protected Methods

		protected override string GetStringForOutput()
		{
			var sb = new StringBuilder();
			return sb.ToString();
		}

		#endregion Protected Methods

	}
}