# load libraries
library(stringr)
library(GenomicRanges)
library(Rsamtools)

lists_genes = list(
  "SETMAR" =                     c("Hsap_ENSP00000373354.3"),
  "Chromo: rve fusions, HG3"   = c("Rhidel_EIE76626","Auraur_scaffold1257.g13.t1","Exapal_XP_020910819.1","Notgen_g14498.t1","Ptyfla_pfl_40v0_9_20150316_1g16730.t1"),
  "Chromo: RVT fusions, HG2.1" = c("Auraur_scaffold921.g6.t1","Mlei_ML1025_g4_t1","Sros_EGD82126","Cowc_KJE97446","Cavfas_EGG14977","Hetalb_XP_020435614.1","Cvel_21008.t1-p1"),
  "Chromo: RVT fusions, HG6.0" = c("Mlei_ML1025_g4_t1"),
  "Chromo: CHD1/2"             = c("Morvir_scaffold102.g15.t1"),
  "PHD: PYGO1/2-like" =          c("Auraur_scaffold1329.g6.t1","Nvec_NVE15670","Morvir_scaffold1.g18.t2","Spur_XP_030841878.1"),
  "PHD: Tc2 fusions, HG9.0" =    c("Apla_XP_022090421.1","Notgen_g20558.t1"),
  "PHD: RVT fusions, HG9.1" =    c("Crei_PNW79026"),
  "PHD: MULE fusions, HG9.1" =   c("Apla_XP_022111688.1","Adig_XP_015759311.1","Auraur_scaffold1574.g2.t1"),
  "MBT: PHF20L1" =               c("Tcas_XP_972506.1"),
  "Bromo" =                      c("Dpul_EFX85712"),
  "TUDOR" =                      c("Morvir_scaffold240.g28.t1","Morvir_scaffold21.g164.t3","Auraur_scaffold72.g1.t1")
)




bam_fn = "../../fusions-genome-alignments/work/"

pdf("results_TEfusions/plot_pileups.pdf", height = 8, width = 12)
layout(mat = matrix(1:30, ncol = 6, byrow = TRUE))
par(mar=c(3, 1, 3, 1))
for (set in names(lists_genes)) {
  
  genes = lists_genes[[set]]
  
  for (gei in genes) {

    # define species
    spi = stringr::str_split(gei, "_", simplify = T)[,1]
    if (spi == "Auraur") { spg = "Aaur" } else { spg = spi }
    message(gei)
    
    # load tx to gene mapping
    tx2ge = read.table(sprintf("results_TEfusions/data/fusions.%s.diamond.csv", spi))

    # get transcript
    txi = tx2ge [ tx2ge[,1] == gei, 2 ] [ 1 ]
    
    # load alignments
    gei_gr = GenomicRanges::GRanges(sprintf("%s:1-10000", txi))
    gei_param = Rsamtools::ScanBamParam(which = gei_gr, what = scanBamWhat())
    bai = Rsamtools::scanBam(file = sprintf("%s/%s/mapping_bwa_merge/merge.bwa.bam", bam_fn, spg), param = gei_param)
    
    # alignment coordinates
    bai_start = unlist(lapply(1:length(bai),  function(b) bai[[b]]$pos ))
    bai_width = unlist(lapply(1:length(bai),  function(b) bai[[b]]$qwidth ))
    bai_stop  = bai_start + bai_width

    # keep non-NA reads    
    keep_reads = !is.na(bai_start) & !is.na(bai_stop)
    bai_start  = bai_start [ keep_reads ]
    bai_stop   = bai_stop [ keep_reads ]
    
    # granges
    bai_gr = GenomicRanges::GRanges(seqnames = txi, IRanges::IRanges(start = bai_start, end = bai_stop))
    
    # read domain coordinates
    arq = read.table(sprintf("results_TEfusions/data/fusions.%s.Pfamscan.seqs.csv", spi))
    arq = arq [ arq[,1] == gei, ]
    dom_start = arq[,2] * 3
    dom_stop  = arq[,3] * 3
    dom_names = arq[,5]
    
    # # granges of domain boundaries
    # dom_gr = GenomicRanges::GRanges(seqnames = txi, IRanges::IRanges(start = sort(c(dom_start,dom_stop)), width = 1))
    # 
    # # overlap between domain boundaries and reads
    # ovs = GenomicRanges::findOverlaps(bai_gr, dom_gr)
    # bai_f_gr = bai_gr [ ovs@from ]
    # bai_f_start = GenomicRanges::start(bai_f_gr)
    # bai_f_stop  = GenomicRanges::end(bai_f_gr)
    
    # granges of interdomain regions
    dom_gr = GenomicRanges::GRanges(
      seqnames = txi,
      IRanges::IRanges(
        start = dom_stop [ 1:length(dom_stop) - 1 ],
        end   = dom_start [ 2:length(dom_stop) ]))

    # overlap between interdomain regions and reads
    ovs = GenomicRanges::findOverlaps(bai_gr, dom_gr)
    bai_f_gr = bai_gr [ ovs@from ]
    bai_f_start = GenomicRanges::start(bai_f_gr)
    bai_f_stop  = GenomicRanges::end(bai_f_gr)
    
    # plot empty
    plot(
      x=NA, y=NA, 
      xlim = c(0,max(c(bai_stop,2000), na.rm = TRUE)), 
      ylim = c(-60,1000), 
      xlab = "", ylab = "", yaxt = "n", 
      frame.plot = FALSE,
      main = sprintf("%s\n%s",gei,set), las = 2,
      cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7)
    
    # plot reads
    for (r in 1:length(bai_f_start)) {
      ry = 1000/ length(bai_f_start) * r
      segments(x0=bai_f_start[r], x1=bai_f_stop[r], y0=ry, y1=ry, col = "lightblue3", lwd = 2, lend = 1)
    }
    
    # plot domains
    segments(x0=dom_start, x1=dom_stop, y0=-30, y1=-30, col = "hotpink3", lwd = 4, lend = 1)
    text(x = dom_start, y = -60, dom_names, cex = 0.7, col = "hotpink4", pos = 4)
    abline(v=c(dom_start, dom_stop), lty = 2, col = "hotpink1")
    
    
    
  }
  
}

dev.off()