# Normalize the two sample columns in a matched normal/tumor VCF.
#
# Required variables:
#   normal_sample  unique sample ID for column 10
#   tumor_sample   unique sample ID for column 11
# Optional variables:
#   emit_metadata  write the normalization metadata line (default: 1)

BEGIN {
    FS = OFS = "\t"
    header_count = 0
    record_count = 0
    status = 0

    if (emit_metadata == "") {
        emit_metadata = 1
    }

    if (normal_sample !~ /^[A-Za-z0-9_.-]+$/) {
        fail("invalid normal sample ID: " normal_sample, 20)
    }
    if (tumor_sample !~ /^[A-Za-z0-9_.-]+$/) {
        fail("invalid tumor sample ID: " tumor_sample, 21)
    }
    if (normal_sample == tumor_sample) {
        fail("normal and tumor sample IDs must differ", 22)
    }
}

function fail(message, code) {
    print "ERROR: normalize_two_sample_vcf: " message > "/dev/stderr"
    status = code
    exit code
}

/^##/ {
    print
    next
}

/^#CHROM/ {
    header_count++
    if (header_count != 1) {
        fail("multiple #CHROM headers", 23)
    }
    if (NF != 11) {
        fail("expected exactly two sample columns; found " (NF - 9), 24)
    }

    original_normal = $10
    original_tumor = $11
    if (emit_metadata) {
        print "##SVCFitSampleNameNormalization=<Normal=" normal_sample ",Tumor=" tumor_sample ">"
    }
    $10 = normal_sample
    $11 = tumor_sample
    print
    next
}

/^#/ {
    print
    next
}

{
    if (header_count != 1) {
        fail("variant record encountered before #CHROM header", 25)
    }
    if (NF != 11) {
        fail("record " (record_count + 1) " has " NF " fields; expected 11", 26)
    }
    record_count++
    print $0
}

END {
    if (status != 0) {
        exit status
    }
    if (header_count != 1) {
        print "ERROR: normalize_two_sample_vcf: missing #CHROM header" > "/dev/stderr"
        exit 27
    }
    print "Normalized VCF samples: " original_normal " -> " normal_sample \
          "; " original_tumor " -> " tumor_sample \
          "; records=" record_count > "/dev/stderr"
}
