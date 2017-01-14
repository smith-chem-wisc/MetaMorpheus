using System.Text;

namespace InternalLogicTaskLayer
{
	class MySearchTaskResults : MyTaskResults
	{
		public MySearchTaskResults(MyTaskEngine s) : base(s)
		{
		}

		protected override string GetStringForOutput()
		{
			var sb = new StringBuilder();
			return sb.ToString();
		}
	}
}