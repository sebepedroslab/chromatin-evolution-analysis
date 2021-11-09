SELECT MasterProteinAccessions, Sequence, Modifications, ptmRSBestSiteProbabilities FROM TargetPsms WHERE ptmRSBestSiteProbabilities IS NOT NULL AND MasterProteinAccessions IS NOT NULL;
