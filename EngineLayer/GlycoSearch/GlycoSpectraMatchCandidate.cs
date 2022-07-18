using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.GlycoSearch
{
    public struct GlycoSpectraMatchCandidate
    {
        public GlycoSpectraMatchCandidate(int pepId, double score, double pScore, double xcorrScore, GlycoType gType, int rank, int[] modPos, string[] modMotifs, int[] nPos, List<int> glycanBoxIds)
        {
            PepId = pepId;
            Score = score;
            PScore = pScore;
            XcorrScore = xcorrScore;
            GType = gType;
            Rank = rank;
            ModPos = modPos;
            ModMotifs = modMotifs;
            NPos = nPos;
            GlycanBoxIds = glycanBoxIds;
        }

        public int PepId { get; }

        public double Score { get; }

        public double PScore { get; }

        public double XcorrScore { get; }
        public GlycoType GType { get; }

        public int Rank { get; }

        public int[] ModPos { get; }

        public string[] ModMotifs { get; }

        public int[] NPos { get; }

        public List<int> GlycanBoxIds{get;}



    }
}
