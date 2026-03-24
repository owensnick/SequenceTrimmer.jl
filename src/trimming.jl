

"""
    Read a block of fastq reads into Vector{String}[]
"""
function read_block!(io, data::Vector{Vector{String}}, blocksize=length(data))
    r = 0
    for i = 1:blocksize
        eof(io) && break
        data[i][1] = readline(io)
        data[i][2] = readline(io)
        data[i][3] = readline(io)
        data[i][4] = readline(io)
        r += 1
    end
    r
end


trimdata(data, n) = @view data[1:min(length(data), n)]

"""
    Write a block of fastq reads from Vector{String}[]
    Writes trimmed sequence and quality string
"""
function write_block(io, data, rn, trim_data)
    for i = 1:rn
        ind = 1:trim_data[2, i]
        write(io, data[i][1], "\n")
        write(io, trimdata(data[i][2], trim_data[2, i]), "\n")
        write(io, data[i][3], "\n")
        write(io, trimdata(data[i][4], trime_data[2, i]), "\n")
    end
end


"""
    Count the number of Ambiguous bases in a sequence
    i.e. count number of N's
"""
count_n(seq) = count(isambiguous, seq)



"""
    approxsuffix(seq, pattern, max_mismatches)
    Find approximate suffix match of pattern in seq
"""
@inline function approxsuffix(seq, pattern, max_mismatches)
    n = length(seq)
    m = length(pattern)

    # scan from right to left (like findprev)
    for endpos in n:-1:m
        start = endpos - m + 1
        mismatches = 0

        @inbounds for i = 1:m
            if seq[start+i-1] != pattern[i]
                mismatches += 1
                mismatches > max_mismatches && break
            end
        end

        if mismatches ≤ max_mismatches
            return start:endpos
        end
    end

    return 1:0
end


"""
    Find reverse complement alignment of two sequences
    find 
"""
function align_rc(seq1, seq2, align_length=20, max_distance=1, max_n=5)
    rseq2 = reverse_complement(seq2)

    (seq1 == rseq2) && return 1:0
    rc2 = rseq2[(end-align_length+1):end]
    (count_n(rc2) >= max_n) && return 1:0

    cand_ind = approxsuffix(seq1, rc2, max_distance)

    isempty(cand_ind) && return 1:0
    (cand_ind.stop == length(seq1)) && return 1:0

    return 1:cand_ind.stop
end




"""
    Increment a pwm with a sequence
"""
function add_pwm!(PWM, seq)
    for i = 1:min(size(PWM, 2), length(seq))
        PWM[5, i] += 1
        if seq[i] == DNA_A
            PWM[1, i] += 1
        elseif seq[i] == DNA_C
            PWM[2, i] += 1
        elseif seq[i] == DNA_G
            PWM[3, i] += 1
        elseif seq[i] == DNA_T
            PWM[4, i] += 1
        end
    end
end


"""
    Estimate over respresented sequences in the portion surrounding a reverse complement match between read 1 and 2
"""
function estimate_libtype(fq1, fq2, num_reads=5, rl=50, sl=0)

    f1 = GzipDecompressorStream(open(fq1))
    f2 = GzipDecompressorStream(open(fq2))

    PWM_L = zeros(5, rl)
    PWM_R = zeros(5, rl)


    nr = 0
    tls = 0
    while !eof(f1) && !eof(f2) && (nr != (num_reads + sl))
        id1 = readline(f1)
        read1 = readline(f1)
        qid1 = readline(f1)
        qs1 = readline(f1)

        id2 = readline(f2)
        read2 = readline(f2)
        qid2 = readline(f2)
        qs2 = readline(f2)

        nr += 1
        (nr <= sl) && continue

        seq1 = LongDNA{4}(read1)
        seq2 = LongDNA{4}(read2)
        ind = align_rc(seq1, seq2)

        if ind != 1:0
            tls += 1
            adapter1 = seq1[(ind.stop+1):end]
            adapter2 = seq2[(ind.stop+1):end]

            add_pwm!(PWM_L, adapter1)
            add_pwm!(PWM_R, adapter2)
        end
    end

    close(f1)
    close(f2)

    L_consen, L_adapter_type, L_ic = pwm_summary(PWM_L, "[TFR]\tLeft  $(tls) sequences :")
    R_consen, R_adapter_type, R_ic = pwm_summary(PWM_R, "[TFR]\tRight $(tls) sequences :")

    L_consen, R_consen
end

"""
    Get consensus base from a vector of probabilites
"""
function getbase(v, seq=[DNA_A, DNA_C, DNA_G, DNA_T], thresh=0.25)
    m, mi = findmax(v)
    if m < thresh
        return DNA_Gap
    else
        seq[mi]
    end
end

"""
    Get consensus from PWM
"""
consensus(PWM, l=size(PWM, 2)) = LongDNA{4}([getbase(view(PWM, 1:4, i)) for i = 1:l])

"""
    Determine is seq `a` is an instance of seq `b` upto mismatches
"""
function isinstance(a, b, mm=3)
    tmm = 0
    for i = 1:min(length(a), length(b))
        tmm += ifelse(a[i] == b[i], 0, 1)
    end
    tmm ≤ mm
end

"""
    classify_consensus sequences
"""
function classify_consensus(cl)
    isinstance(cl, dna"AGATCGGAAGAGC") && return :Illumina
    isinstance(cl, dna"CTGTCTCTTATACACATCT") && return :Nextera
    return :Unknown
end

"""
    Print pwm summary
"""
function pwm_summary(PWM, label="Adapt")
    consen = consensus(PWM)
    adapter_type = classify_consensus(consen)
    ic = info_content(PWM)
    println(label, " consensus: ", consen, "\tType: ", adapter_type, "\tIC: ", ic, " bits")
    consen, adapter_type, ic
end

"""
    Calculate information content of a PWM in bits
"""
function info_content(PWM)
    ts = PWM[5, :]
    en = 3.0 ./ (2 * ts * log(2))
    F = PWM[1:4, :] ./ ts'
    H = -sum(F .* log2.(F), dims=1)
    R = log2(4) .- (H + en')
    sum(filter(r -> !isnan(r) && !isinf(r), R))
end

"""
    Find the adapter location within a sequence upto a mismatch threshold
"""
function find_adapt_location(seq1, seq2, L_match, R_match, thresh=3 / 13, tail=7)
    n = min(length(seq1), length(seq2))
    amm = 0
    bmm = 0

    min_p = min(length(L_match), length(R_match))

    for i = 1:(n-tail)
        amm = 0
        bmm = 0
        loop_len = min(min_p, n - i + 1)


        try
            for j = 1:loop_len
                amm += ifelse((seq1[i+j-1] != L_match[j]), 1, 0)
                bmm += ifelse((seq2[i+j-1] != R_match[j]), 1, 0)
            end
        catch err
            @show seq1
            @show seq2
            @show n
            @show L_match
            @show R_match
            throw(err)
        end

        if 0.5 * (amm + bmm) / loop_len < thresh
            return i - 1, amm, bmm
        end
    end
    return 0, 0, 0
end

"""
    Determine if seqA ≈ rc(seqB) upto mismatches
"""
function rcmm(seqA, seqB)
    rb = reverse_complement(seqB)
    mm = mismatches(seqA, rb) - count(isambiguous, seqA, rb)
    max(0, mm)
end


"""
    Trim Reads function
    Takes pairs of reads and trimm


    ## Trimming procedure
    1. Attempt to find L_match and R_Match adapters in R1 and R2 respectively
    2. Check remaining sequence is RC
    3. If adapter discovery fails attempt to find RC region
    4. Check that remaining region contains adapters

"""
function trim_reads(readR1, readR2, L_match, R_match, align_length, max_distance, mismatch_rate, max_n)

    # Convert string to DNA sequence
    seq1 = LongDNA{4}(readR1[2])
    seq2 = LongDNA{4}(readR2[2])

    ### Find the adapter location
    adapt_location, r1mm, r2mm = find_adapt_location(seq1, seq2, L_match, R_match, mismatch_rate)
    trim_ind = 1:adapt_location

    if adapt_location > 0 && (rcmm(seq1[trim_ind], seq2[trim_ind])) / length(trim_ind) < mismatch_rate
        return 1, adapt_location
    else
        ### find RC
        trim_ind = align_rc(seq1, seq2, align_length, max_distance, max_n)

        if !isempty(trim_ind)
            adaptind = (trim_ind.stop+1):length(seq1)

            mm1 = mismatches(seq1[adaptind], L_match)
            mm2 = mismatches(seq1[adaptind], R_match)

            compare_length = min(length(adaptind), length(L_match))
            if (mm1 + mm2) / (2 * compare_length) < mismatch_rate
                return 2, trim_ind.stop
            end
        end
    end
    0, length(seq1)
end

"""
    Threaded fastq trimmer

"""
function trim_fastq_threads(fq_R1, fq_R2, align_length=20, max_distance=1, mismatch_rate=1 / 5, max_n=5, maxreads=-1, blocksize=100, adapt_length=13, method=:estimate, est_consenus_reads=100000)


    # Set output
    outfile_R1 = replace(fq_R1, ".fastq.gz" => ".trim.fastq.gz")
    outfile_R2 = replace(fq_R2, ".fastq.gz" => ".trim.fastq.gz")
    outfile_R1 = replace(outfile_R1, ".fq.gz" => ".trim.fq.gz")
    outfile_R2 = replace(outfile_R2, ".fq.gz" => ".trim.fq.gz")

    (fq_R1 == outfile_R1) && error("In/Out Equal: $fq_R1 $outfile_R1")
    (fq_R2 == outfile_R2) && error("In/Out Equal: $fq_R2 $outfile_R2")



    ### Print parameters
    println("[TFR]\tFrag align trim       :\t$fq_R1, $fq_R2")
    println("[TFR]\talign_length          :\t", align_length)
    println("[TFR]\tmax_distance          :\t", max_distance)
    println("[TFR]\tmismatch_rate         :\t", mismatch_rate)
    println("[TFR]\tmax_n                 :\t", max_n)
    println("[TFR]\tmaxreads              :\t", maxreads)
    println("[TFR]\tblocksize             :\t", blocksize)
    println("[TFR]\tnum_threads           :\t", Threads.nthreads())
    println("[TFR]\test_consensus_reads   :\t", est_consenus_reads)
    println("[TFR]\tadapt_length          :\t", adapt_length)


    ### Find consensus sequences to trim from
    if method == :estimate # Estimate the libtype from overepresented sequences
        println("[TFR]\tEstimating adapters   :    ")
        L_consen, R_consen = estimate_libtype(fq_R1, fq_R2, est_consenus_reads)
    elseif method == :Illumina
        println("[TFR]\tAdapter Search        :\t", method)
        L_consen = dna"AGATCGGAAGAGC"
        R_consen = dna"AGATCGGAAGAGC"
    elseif method == :Nextera
        println("[TFR]\tAdapter Search        :\t", method)
        L_consen = dna"CTGTCTCTTATACACATCT"
        R_consen = dna"CTGTCTCTTATACACATCT"
    else
        println("[TFR]\tEstimating adapters   :    ")
        L_consen, R_consen = estimate_libtype(fq_R1, fq_R2, est_consenus_reads)
    end


    L_match = L_consen[1:min(length(L_consen), adapt_length)]
    R_match = R_consen[1:min(length(R_consen), adapt_length)]

    println("[TFR]\tL_match               :\t", L_match)
    println("[TFR]\tR_match               :\t", R_match)
    println("[TFR]\tWriting               :\t$outfile_R1, $outfile_R2")
    flush(stdout) ### Flushes IO so can read the output parameters whilst trimming
    starttime = time()


    #### Open readers and writers
    reader_r1 = open(fq_R1) |> GzipDecompressorStream
    reader_r2 = open(fq_R2) |> GzipDecompressorStream

    writer_r1 = open(outfile_R1, "w") |> GzipCompressorStream
    writer_r2 = open(outfile_R2, "w") |> GzipCompressorStream

    #### Datastructure to read fastq files and trim them
    datablock_R1 = [Vector{String}(undef, 4) for i = 1:blocksize]
    datablock_R2 = [Vector{String}(undef, 4) for i = 1:blocksize]

    ### Trimming stats
    temp_trim_data = Vector{Tuple{Int,Int}}(undef, blocksize)
    trim_data = zeros(Int, 2, blocksize)
    trim_stats = counter(Tuple{Int,Int})

    total_reads = 0
    max_trim_length = 0

    while true
        ### Read block of reads
        rn1 = read_block!(reader_r1, datablock_R1, blocksize)
        rn2 = read_block!(reader_r2, datablock_R2, blocksize)
        (rn1 != rn2) && error("FASTQ mismatch")
        (rn1 == 0) && break

        ### Threaded trimming of the block
        Threads.@threads for i = 1:rn1
            trim_data[1, i], trim_data[2, i] = trim_reads(datablock_R1[i], datablock_R2[i], L_match, R_match, align_length, max_distance, mismatch_rate, max_n)
        end

        #### Accumulate trim stats
        for i = 1:rn1
            if trim_data[1, i] > 0
                max_trim_length = max(trim_data[2, i], max_trim_length)
            end
            push!(trim_stats, (trim_data[1, i], trim_data[2, i]))
        end

        #### Write Reads
        write_block(writer_r1, datablock_R1, rn1, trim_data)
        write_block(writer_r2, datablock_R2, rn1, trim_data)

        total_reads += rn1
        (maxreads != -1) && (total_reads > maxreads) && break
    end

    ### Close readers and writers
    close(reader_r1)
    close(reader_r2)
    close(writer_r1)
    close(writer_r2)


    ### Print stats
    total_trimmmed = 0

    F = zeros(Int, 3, max_trim_length)

    for ((trim_type, frag), trim_count) ∈ trim_stats
        (trim_type == 0) && continue
        total_trimmmed += trim_count

        F[1, frag] += trim_count
        F[trim_type+1, frag] += trim_count
    end


    println("[TFR]\tComplete in           :\t", time() - starttime, " seconds")
    println("[TFR]\ttotal_reads           :\t", total_reads)
    println("[TFR]\ttotal_trimmed         :\t", total_trimmmed, "\t", total_trimmmed / total_reads)
    println("[TFR]\tTrimmed frag dist     :\t")


    for i = 1:max_trim_length
        println(i, "\t", F[1, i], "\t", F[2, i], "\t", F[3, i])
    end

end