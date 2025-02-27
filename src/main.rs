use anyhow::Result;
use boxcar::Vec as V;
use clap::Parser;
use core::array::from_fn;
use mimalloc::MiMalloc;
use niffler::send::from_path;
use packed_seq::{PackedSeqVec, Seq, SeqVec};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rayon::{current_num_threads, ThreadPoolBuilder};
use regex::bytes::{Regex, RegexBuilder};
use rustc_hash::FxBuildHasher;
use seq_io::fasta;
use seq_io_parallel::{MinimalRefRecord, ParallelProcessor, ParallelReader};
use simd_minimizers::minimizer_and_superkmer_positions;
use std::collections::HashSet;
use std::path::Path;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

type KT = u64;
type SKT = u128; // together as one

const SHARD_BASES: usize = 5;
const SHARDS: usize = 1 << (2 * SHARD_BASES);
const SKLEN_BITS: usize = 6;
const SKLEN_MASK: SKT = (1 << SKLEN_BITS) - 1;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA, possibly compressed)
    #[arg(short, long)]
    input: String,
    /// K-mer size (up to 32)
    #[arg(short)]
    k: usize,
    /// Minimizer size
    #[arg(short, default_value_t = 21)]
    m: usize,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

#[derive(Clone)]
pub struct SuperkmerCollector<'a> {
    k: usize,
    m: usize,
    buckets: &'a [V<SKT>; SHARDS],
    match_n: &'a Regex,
    match_newline: &'a Regex,
    min_pos_vec: Vec<u32>,
    sk_pos_vec: Vec<u32>,
}

impl ParallelProcessor for SuperkmerCollector<'_> {
    fn process_record<'a, Rf: MinimalRefRecord<'a>>(&mut self, record: Rf) -> Result<()> {
        let w = self.k - self.m + 1;
        for raw_seq in self
            .match_n
            .split(record.ref_seq())
            .filter(|&s| s.len() >= self.k)
        {
            let mut packed_seq = PackedSeqVec::default();
            for line in self.match_newline.split(raw_seq) {
                if !line.is_empty() {
                    packed_seq.push_ascii(line);
                }
            }
            let len = packed_seq.len();
            if len >= self.k {
                self.min_pos_vec.clear();
                self.min_pos_vec.reserve(len * 5 / 2 / (w + 1));
                self.sk_pos_vec.clear();
                self.sk_pos_vec.reserve(len * 5 / 2 / (w + 1));
                minimizer_and_superkmer_positions(
                    packed_seq.as_slice(),
                    self.m,
                    w,
                    &mut self.min_pos_vec,
                    &mut self.sk_pos_vec,
                );
                self.min_pos_vec.push(u32::MAX);
                self.sk_pos_vec.push((len - (self.k - 1)) as u32);
                let mut min_pos = self.min_pos_vec[0];
                let mut sk_pos = self.sk_pos_vec[0];
                for (&next_min_pos, &next_sk_pos) in
                    self.min_pos_vec.iter().zip(self.sk_pos_vec.iter()).skip(1)
                {
                    let shard_range = (min_pos as usize)..(min_pos as usize + SHARD_BASES);
                    let shard = packed_seq.slice(shard_range).to_word();
                    let sk_range = (sk_pos as usize)..((next_sk_pos as usize) + self.k - 1);
                    let sk_mid = (sk_range.start + sk_range.end) / 2;
                    let left = packed_seq.slice(sk_range.start..sk_mid).to_word() as SKT;
                    let right = packed_seq.slice(sk_mid..sk_range.end).to_word() as SKT;
                    let skmer = (((right << (2 * (sk_mid - sk_range.start))) | left) << SKLEN_BITS)
                        | (sk_range.len() as SKT); // little-endian order
                    self.buckets[shard].push(skmer);
                    min_pos = next_min_pos;
                    sk_pos = next_sk_pos;
                }
            }
        }
        Ok(())
    }
}

fn collect_superkmers<P: AsRef<Path>>(
    k: usize,
    m: usize,
    path: P,
    threads: usize,
) -> [V<SKT>; SHARDS] {
    let match_n = RegexBuilder::new(r"[N]+")
        .case_insensitive(true)
        .unicode(false)
        .build()
        .unwrap();
    let match_newline = RegexBuilder::new(r"[\r\n]+")
        .unicode(false)
        .build()
        .unwrap();
    let buckets = from_fn(|_| V::with_capacity(1 << 14));
    let processor = SuperkmerCollector {
        k,
        m,
        buckets: &buckets,
        match_n: &match_n,
        match_newline: &match_newline,
        min_pos_vec: vec![],
        sk_pos_vec: vec![],
    };
    let (reader, _) = from_path(path).expect("Failed to open input file");
    let reader = fasta::Reader::new(reader);
    reader.process_parallel(processor, threads).unwrap();
    buckets
}

fn main() {
    let args = Args::parse();
    let k = args.k;
    assert!(k <= 32);
    let m = args.m;
    assert!(m <= k);
    let w = k - m + 1;
    let path = args.input;
    let threads = if let Some(t) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
        t
    } else {
        current_num_threads()
    };
    let buckets = collect_superkmers(k, m, path, threads);
    let kmer_mask = (1u128 << (2 * k)) - 1;
    let count: usize = buckets
        .into_par_iter()
        .map(|v| {
            let mut set =
                HashSet::with_capacity_and_hasher(v.count() * (w + 1) * 3 / 5, FxBuildHasher);
            for skmer in v {
                let len = (skmer & SKLEN_MASK) as usize;
                let skmer = skmer >> SKLEN_BITS;
                for i in 0..(len - k + 1) {
                    let kmer = ((skmer >> (2 * i)) & kmer_mask) as KT; // start with low bits
                    set.insert(kmer);
                }
            }
            set.len()
        })
        .sum();
    println!("{count} distinct {k}-mers");
}
