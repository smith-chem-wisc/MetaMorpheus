using System;
using System.Data.Entity.Core.Objects.DataClasses;

namespace TaskLayer
{
    public enum MassDiffAcceptorType
    {
        Exact,
        OneMM,
        TwoMM,
        ThreeMM,
        ModOpen,
        Open,
        Custom,
        PlusOrMinusThreeMM,
        MostAbundant_Exact,
        MostAbundant_PlusMinusOne,
        MostAbundant_PlusMinusTwo
    }

    public static class MassDiffAcceptorTypeExtensions
    {
        public static bool IsMostAbundant(this MassDiffAcceptorType type)
        {
            switch (type)
            {
                case MassDiffAcceptorType.Exact:
                case MassDiffAcceptorType.OneMM:
                case MassDiffAcceptorType.TwoMM:
                case MassDiffAcceptorType.ThreeMM:
                case MassDiffAcceptorType.ModOpen:
                case MassDiffAcceptorType.Open:
                case MassDiffAcceptorType.Custom:
                case MassDiffAcceptorType.PlusOrMinusThreeMM:
                    return false;

                case MassDiffAcceptorType.MostAbundant_Exact:
                case MassDiffAcceptorType.MostAbundant_PlusMinusOne:
                case MassDiffAcceptorType.MostAbundant_PlusMinusTwo:
                    return true;

                    default:
                    throw new ArgumentException($"Unknown MassDiffAcceptorType: {type}");
            }
        }
    }
}
