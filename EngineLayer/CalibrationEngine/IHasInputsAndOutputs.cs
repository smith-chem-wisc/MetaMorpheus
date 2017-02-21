namespace EngineLayer.Calibration
{
    public interface IHasInputsAndOutputs
    {

        #region Public Properties

        double[] inputs { get; }
        double label { get; }

        #endregion Public Properties

    }
}