#!/usr/bin/env Rscript
##
## Build SVclone's three input files for one chrX simulation condition.
##
## SVclone was designed around FACETS, which gives integer allele-specific copy number per segment.
## We have no FACETS here and would not trust it if we did: it derives copy number from heterozygous
## germline SNPs, and a male X has none. That is the same fact the whole hemizygous estimator exists
## to work around. SVCFit is given DNAcopy segments instead, so out of fairness SVclone gets the
## same DNAcopy segments, converted to the format it can read.
##
## WHAT THE CONVERSION COSTS SVCLONE, stated plainly because it is not neutral. DNAcopy gives a
## continuous tumour-to-normal depth ratio, cn_bar, with no allelic split. SVclone's cnv.txt wants
## integers. Rounding is therefore forced by SVclone's input format -- it is not a choice made to
## disadvantage it, and SVCFit deliberately does NOT round (see METHODS-DRAFT.md: cn_bar "is a
## continuous measurement and is not rounded to an integer copy number").
##
## THE FLOOR. round(cn_bar) is 0 for any segment below 0.5, and a segment at zero total copies is
## something SVclone may drop outright, which would silently shrink its denominator. Measured over
## the 45 REP=0 conditions: 291 of 4188 segments (6.95%) fall below 0.5, covering 25.7 Mb of
## 6976 Mb (0.37% of bases). The minimum cn_bar anywhere is 0.187, so nothing is near zero and the
## floor is effectively binary -- 0 or 1, with no useful value between. Those segments sit entirely
## in e1 and e2 (e4 has NONE, because its deletions lie on amplified backgrounds where cn_bar =
## kappa - s stays above 0.5) and concentrate at high purity: 12.96% at p60, 15.56% at p80.
##
## Default FLOOR = 1, so no locus is lost. The alternative loses segments specifically at high
## purity, where both methods should do best, distorting SVclone's purity profile in a direction
## that cannot be reasoned about afterwards. Every clamped segment is logged so the count is a
## reported number rather than a hidden one. Pass --floor 0 to run the sensitivity check.
##
## DELIBERATE DEVIATION FROM visor_replicates/scripts/helper/make_input.R: that script keeps only
## FILTER == "PASS". SVCFit's chrX path (13_chrx_score_svcfit.R:94) reads the svtyper VCF and does
## not filter on FILTER at all, so a PASS filter here would hand SVclone 32 of 42 calls against
## SVCFit's 42 and lose it detections it never had a chance at. Both methods get the same call set,
## from the same file.
##
## Usage:
##   20_chrx_svclone_input.R --vcf <svt_*.vcf> --seg <chrx_seg_*.csv> --sample <name> \
##                           --purity <0-1> --out <dir> [--floor 1] [--ploidy 2]

suppressMessages(library(dplyr))

## Arguments are parsed with base R rather than optparse ON PURPOSE. optparse is installed in the
## `visor` environment only -- not in `svclone`, which is the environment this runs under, and not
## in `svcfit`. visor_replicates/scripts/helper/make_input.R does `library(optparse)` and is invoked
## by the SVclone stage. Depending on a third environment just to read six flags
## is not worth the coupling.
.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = NULL) {
  i <- match(flag, .args)
  if (is.na(i)) return(default)
  if (i == length(.args)) stop(flag, " given with no value", call. = FALSE)
  .args[i + 1L]
}
opt <- list(vcf    = .get("--vcf"),    seg    = .get("--seg"),
            sample = .get("--sample"), purity = as.numeric(.get("--purity", NA)),
            out    = .get("--out"),    floor  = as.integer(.get("--floor", "1")),
            ploidy = as.numeric(.get("--ploidy", "2")))
if (any(!.args[startsWith(.args, "--")] %in%
        c("--vcf","--seg","--sample","--purity","--out","--floor","--ploidy")))
  stop("unrecognised flag: ",
       paste(setdiff(.args[startsWith(.args, "--")],
                     c("--vcf","--seg","--sample","--purity","--out","--floor","--ploidy")),
             collapse = ", "), call. = FALSE)

for (r in c("vcf", "seg", "sample", "purity", "out"))
  if (is.null(opt[[r]])) stop("--", r, " is required", call. = FALSE)
if (!file.exists(opt$vcf)) stop("VCF not found: ", opt$vcf, call. = FALSE)
if (!file.exists(opt$seg)) stop("segment file not found: ", opt$seg, call. = FALSE)
if (!is.finite(opt$purity) || opt$purity <= 0 || opt$purity > 1)
  stop("--purity must be a fraction in (0,1]; got ", opt$purity, call. = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

## ---- 1. SV breakpoints -------------------------------------------------------------------------
## Column derivation and the INV position swap are transcribed from make_input.R so the two
## benchmarks build SVclone's input the same way. The FILTER == "PASS" step is deliberately absent;
## see the header.
raw <- read.table(opt$vcf, quote = "\"", stringsAsFactors = FALSE)
colnames(raw) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","normal","tumor")[seq_len(ncol(raw))]

sv <- raw %>%
  rename(chr1 = CHROM, pos1 = POS) %>%
  mutate(chr2  = gsub(".*(chr[0-9XYM]+).*", "\\1", ALT),
         chr2  = ifelse(grepl("chr", chr2), chr2, chr1),
         pos2  = gsub(".*:(\\d+).*", "\\1", ALT),
         pos2  = ifelse(grepl("^\\d+$", pos2), pos2, gsub(".*END=(\\d+).*", "\\1", INFO)),
         pos2  = suppressWarnings(as.integer(pos2)),
         class = gsub(".*SVTYPE=(\\w+).*", "\\1", INFO),
         pp    = pos1,
         pos1  = ifelse(class == "INV", pos2, pos1),
         pos2  = ifelse(class == "INV", pp,   pos2),
         chr1  = ifelse(grepl("^chr", chr1), chr1, paste0("chr", chr1)),
         chr2  = ifelse(grepl("^chr", chr2), chr2, paste0("chr", chr2))) %>%
  filter(class != "INS") %>%
  select(chr1, pos1, chr2, pos2)

n_drop_na <- sum(is.na(sv$pos2))
if (n_drop_na) {
  ## An unparseable mate coordinate is dropped rather than guessed, and counted so it cannot pass
  ## silently as a detection SVclone was never offered.
  message(sprintf("[%s] %d SV(s) dropped: mate position unparseable from ALT and INFO",
                  opt$sample, n_drop_na))
  sv <- sv[!is.na(sv$pos2), , drop = FALSE]
}
write.table(sv, file.path(opt$out, paste0("svc_", opt$sample, "_simple.txt")),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

## ---- 2. purity / ploidy ------------------------------------------------------------------------
## Purity is the TRUE value from the condition name. SVclone takes it as a parameter; SVCFit has no
## purity argument at all (calc_svcf's formals contain none) and infers everything from depth and
## reads. The asymmetry favours SVclone and matches the submitted autosomal benchmark, which also
## passed true purity. Say so wherever these results are reported.
write.table(data.frame(sample = opt$sample, purity = opt$purity, ploidy = opt$ploidy),
            file.path(opt$out, paste0("spp_", opt$sample, ".txt")),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

## ---- 3. CNV ------------------------------------------------------------------------------------
seg <- read.csv(opt$seg, stringsAsFactors = FALSE)
need <- c("chrom", "start", "end", "cn_bar")
if (!all(need %in% names(seg)))
  stop("segment file must have columns ", paste(need, collapse = ", "), "; got: ",
       paste(names(seg), collapse = ", "), call. = FALSE)

t_round <- round(seg$cn_bar)
t_total <- pmax(opt$floor, t_round)
clamped <- which(t_round < opt$floor)

cnv <- data.frame(chrom   = ifelse(grepl("^chr", seg$chrom), seg$chrom, paste0("chr", seg$chrom)),
                  start   = seg$start,
                  end     = seg$end,
                  ## NORMAL copy number on a hemizygous chromosome is one copy, no minor allele --
                  ## not the 2/1 make_input.R hardcodes for autosomes. Getting this wrong would
                  ## halve every cellular fraction SVclone reports.
                  n_major = 1L,
                  n_minor = 0L,
                  t_total = as.integer(t_total),
                  ## No allelic split exists to report on a single-copy chromosome.
                  t_minor = 0L,
                  stringsAsFactors = FALSE)

## A ONE-SEGMENT FILE IS SPLIT IN HALF. SVclone's load_cnvs() probes the file with delimiter='\t'
## first; a comma-separated file collapses to a single column, which is what routes it to the
## intended "caveman csv" re-read. But the probe consumes line 1 as a HEADER, so a file with exactly
## one segment leaves zero rows, hits the `len(cnv_df) == 0` branch, returns chr='', and every SV
## then fails to match with "No SV chroms match the CNV input!".
##
## One segment is the CORRECT DNAcopy answer for a copy-neutral chromosome with no detectable
## breakpoint -- e.g. chrX 10001-156030001 at cn_bar 1.001 -- and it arises in roughly one or two of
## the 45 conditions in most replicates, concentrated at low purity. Letting SVclone fail there would
## score it as undetected on precisely the conditions where the copy-number picture is simplest,
## which is a parser artefact being recorded as a method failure.
##
## The split is at the midpoint into two segments carrying the SAME copy number, so it adds no
## information and changes no value SVclone can read -- it only makes the file parseable. Logged, so
## the count of affected conditions is reportable.
if (nrow(cnv) == 1L) {
  mid <- as.integer((cnv$start[1] + cnv$end[1]) %/% 2)
  cnv <- rbind(transform(cnv[1, ], end = mid), transform(cnv[1, ], start = mid + 1L))
  message(sprintf("[%s] single segment split at %d into two of identical copy number, ",
                  opt$sample, mid),
          "so SVclone's loader can read the file; see the note in this script")
}

write.table(cnv, file.path(opt$out, "cnv.txt"),
            row.names = TRUE, col.names = FALSE, sep = ",", quote = FALSE)

## ---- clamp log ---------------------------------------------------------------------------------
if (length(clamped)) {
  log <- data.frame(sample = opt$sample, chrom = cnv$chrom[clamped],
                    start = seg$start[clamped], end = seg$end[clamped],
                    cn_bar = seg$cn_bar[clamped], rounded = t_round[clamped],
                    assigned = as.integer(t_total)[clamped], stringsAsFactors = FALSE)
  write.table(log, file.path(opt$out, paste0("clamped_", opt$sample, ".tsv")),
              row.names = FALSE, sep = "\t", quote = FALSE)
}
message(sprintf("[%s] %d SV(s), %d segment(s), %d clamped to floor=%d (%.2f%% of segments, %.2f%% of bases)",
                opt$sample, nrow(sv), nrow(cnv), length(clamped), opt$floor,
                100 * length(clamped) / max(1, nrow(cnv)),
                100 * sum(as.numeric(seg$end - seg$start)[clamped]) /
                      max(1, sum(as.numeric(seg$end - seg$start)))))
