namespace EngineLayer.Calibration
{
    public interface IHasInputsAndOutputs
    {

        #region Public Properties

        double[] Inputs { get; }
        double Label { get; }

        #endregion Public Properties

    }
}