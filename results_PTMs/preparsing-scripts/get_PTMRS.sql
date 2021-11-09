SELECT MasterProteinAccessions, Sequence, ptmRSBestSiteProbabilities FROM TargetPeptideGroups WHERE ptmRSBestSiteProbabilities IS NOT "" AND MasterProteinAccessions IS NOT NULL;
