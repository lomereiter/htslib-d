module htslib.cram;
import dstep.cram;

import bio.bam.read, bio.bam.utils.value, bio.bam.tagvalue, bio.core.base,
    bio.bam.bai.bin, bio.bam.abstractreader, bio.bam.referenceinfo,
    bio.sam.header, bio.bam.reference;
import std.string, std.exception, std.array, std.range;

class CramException : Exception {
    this(string msg) { super(msg); }
}

auto zeroChecked(alias func, T...)(string err_msg, auto ref T params) {
    int ret = func(params);
    if (!ret)
        throw new CramException(err_msg);
}

struct CramBamReadRange {
    private {
        dstep.cram.cram_fd* _fd;
	IBamSamReader _reader;
    }

    this(cram_fd* fd, CramReader reader) {
        _fd = fd;
	_reader = reader;
        popFront();
    }

    bio.bam.read.BamRead front;
    bool empty;

    void popFront() {
        auto cr = cram_get_seq(_fd);
        if (cr is null) {
            if (_fd.err == 0)
                empty = true;
            else
                throw new CramException("Failed to read CRAM record");
        }
	front = cram2bam(_fd.header, _fd, _fd.ctr.slice, cr,
			 _fd.ctr.curr_rec - 1);
	front.associateWithReader(_reader);
    }

    private BamRead cram2bam(SAM_hdr *bfd, cram_fd *fd, cram_slice *s,
			     cram_record *cr, int rec) {
	char[1024] name_buf = void;
	char* name;
	int name_len;
	if (cr.name_len > 0) {
	    name = cast(char*)(s.name_blk.data + cr.name);
	    name_len = cr.name_len;
	} else {
	    auto p = name_buf.ptr;
	    name = p;
	    auto prefix = fd.prefix;
	    while (*prefix != '\0')
		*p++ = *prefix++;
	    bio.core.utils.format.write(p, ':');
	    bio.core.utils.format.write(p, s.id);
	    bio.core.utils.format.write(p, ':');
	    if (cr.mate_line >= 0 && cr.mate_line < rec)
		bio.core.utils.format.write(p, cr.mate_line);
	    else
		bio.core.utils.format.write(p, rec);
	    name_len = cast(int)(p - name);
	}

	if (cr.rg < -1 || cr.rg >= bfd.nrg)
	    throw new Exception("Read group id is out of range");

	// extra 4 bytes account for RGZ and zero
	auto rg_len = (cr.rg != -1) ? (bfd.rg[cr.rg].name_len + 4) : 0;

	enforce(s.seqs_blk.data !is null && s.qual_blk.data !is null);

	auto name_offset = 32;
	auto cigar_offset = name_offset + name_len + 1;
	auto seq_offset = cigar_offset + cr.ncigar * uint.sizeof;
	auto qual_offset = seq_offset + (cr.len + 1) / 2;
	auto tag_offset = qual_offset + cr.len;

	ubyte[] raw_data = allocate(cr.aux_size + rg_len + tag_offset);
	uint bin_mq_nl = (reg2bin(cr.apos - 1, cr.aend) << 16) |
	    (cr.mqual << 8) | (name_len + 1);
	uint flag_nc = (cr.flags << 16) | cr.ncigar;

	*cast(int*)(raw_data.ptr) = cr.ref_id;
	*cast(int*)(raw_data.ptr + 4) = cr.apos - 1;
        *cast(uint*)(raw_data.ptr + 8) = bin_mq_nl;
	*cast(uint*)(raw_data.ptr + 12) = flag_nc;
	*cast(uint*)(raw_data.ptr + 16) = cast(uint)cr.len;
	*cast(int*)(raw_data.ptr + 20) = cr.mate_ref_id;
	*cast(int*)(raw_data.ptr + 24) = cr.mate_pos - 1;
	*cast(int*)(raw_data.ptr + 28) = cr.tlen;

	raw_data[name_offset .. name_offset+name_len] = (cast(ubyte[])name[0 .. name_len])[];
	auto cigar_ptr = cast(uint*)(raw_data.ptr + cigar_offset);
	cigar_ptr[0 .. cr.ncigar] = (s.cigar + cr.cigar)[0 .. cr.ncigar];
	auto src = cast(string)((s.seqs_blk.data + cr.seq)[0 .. cr.len]);
	auto dst = raw_data.ptr + seq_offset;
	foreach (k; 0 .. cr.len / 2) {
	    auto b1 = Base16(src[2*k]).internal_code;
	    auto b2 = Base16(src[2*k+1]).internal_code;
	    dst[k] = cast(ubyte)((b1 << 4) | b2);
	}
	if ((cr.len & 1) != 0)
	    dst[cr.len / 2] = cast(ubyte)(Base16(src[$-1]).internal_code << 4);

	auto quals = (s.qual_blk.data + cr.qual)[0 .. cr.len];
	raw_data[qual_offset .. $][0 .. cr.len] = quals[];

	auto tag_data = raw_data[tag_offset .. $][0 .. cr.aux_size];
	tag_data[] = (s.aux_blk.data + cr.aux)[0 .. cr.aux_size];

	if (cr.rg != -1) {
	    auto rg_name_str = bfd.rg[cr.rg].name[0 .. rg_len - 4];
	    auto rg_name = Value(cast(string)rg_name_str);
	    (tag_data.ptr + cr.aux_size)[0 .. 2] = cast(ubyte[])"RG";
	    emplaceValue(tag_data.ptr + cr.aux_size + 2, rg_name);
	}

	return BamRead(raw_data);
    }

    private {
	size_t _used;
	ubyte[] _block;
	enum _capa = 65536;
	ubyte[] allocate(size_t size) {
	    if (size + _used <= _capa && _block !is null) {
		auto result = _block[_used .. $][0 .. size];
		_used += size;
		return result;
	    } else {
		_block = uninitializedArray!(ubyte[])(_capa);
		_used = size;
		return _block[0 .. size];
	    }
	}
    }
}

class CramReader : IBamSamReader {
    private {
        string _fn;
        string _mode;
        dstep.cram.cram_fd* _fd;

	SamHeader _header; // (RE-)parsed lazily

	SamHeader parseHeader() {
	    auto ptr = sam_hdr_str(_fd.header);
	    auto len = sam_hdr_length(_fd.header);
	    // FIXME: avoid reparsing header; turn it into an interface
	    // with two implementations?
	    return new SamHeader(cast(string)(ptr[0 .. len]));
	}

	ReferenceSequenceInfo[] _reference_sequences;
	int[string] _reference_sequence_dict;
    }

    this(string filename) {
        _fn = filename;
        _fd = cram_open(toStringz(filename), "r");
        if (_fd == null) {
            throw new CramException("can't open file " ~ filename);
        }

	_reference_sequences.length = _fd.header.nref;
	foreach (k; 0 .. _fd.header.nref) {
	    auto seq = _fd.header.ref_[k];
	    _reference_sequences[k] = ReferenceSequenceInfo(seq.name.to!string,
							    seq.len);
	    _reference_sequence_dict[_reference_sequences[k].name] = cast(int)k;
	}
    }

    ///
    bool hasReference(string reference) {
	return null != (reference in _reference_sequence_dict);
    }

    ///
    string filename() @property const {
	return _fn;
    }

    ///
    bio.bam.reference.ReferenceSequence opIndex(string ref_name) {
	enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " is not present in the header");
	auto ref_id = _reference_sequence_dict[ref_name];
	return ReferenceSequence(null, ref_id, _reference_sequences[ref_id]);
    }

    SamHeader header() @property {
	if (_header is null)
	    _header = parseHeader();
	return _header;
    }

    const(ReferenceSequenceInfo)[] reference_sequences() @property const {
	return _reference_sequences;
    }

    /// no-op
    void assumeSequentialProcessing() {}

    void close() {
        zeroChecked!cram_close("failed to close CRAM file descriptor",
                _fd);
    }

    CramBamReadRange reads() {
        auto fd = _fn == "-" ? _fd : cram_open(toStringz(_fn), "r");
	cram_set_option(fd, cram_option.CRAM_OPT_DECODE_MD);
	cram_set_option(fd, cram_option.CRAM_OPT_NTHREADS, 3);
        return typeof(return)(fd, this);
    }

    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property {
	return inputRangeObject(reads());
    }
}

import std.stdio, std.range;
void main() {
    auto cram = new CramReader("/home/lomereiter/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam.cram");
    writefln("%(%s\n%)",cram.reads.take(1_000_000));
}
