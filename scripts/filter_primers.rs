//! Split primer-assigned from primerless reads (upstream primers.smk
//! filter_primerless_reads).
//!
//! Port of the upstream workflow/scripts/filter_primers.rs (MIT,
//! snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
//!   filter_primers.rs INPUT.bam --primers PRIMERS.bam --primerless PRIMERLESS.bam
//!
//! The snakemake script API is replaced by positional/flag arguments; the
//! stderr redirect (snakemake.redirect_stderr) is dropped because the rule
//! shell captures stderr with `2> {log}`. Reads carrying the `ra` (assigned)
//! tag, the `ma` (multiple assignment) tag, or a secondary alignment to an
//! assigned primary go to the primers output; everything else is primerless.
//!
//! ```cargo
//! cargo-features = ["edition2021"]
//! [dependencies]
//! indexmap = "1.8"
//! noodles = { version = "0.18.0", features = ["bam", "sam", "bgzf"] }
//! ```

use indexmap::IndexMap;
use noodles::bam::{Reader, Record, Writer};
use noodles::bgzf::writer;
use noodles::sam::{header::ReferenceSequence, record::data::field::Tag, Header};
use std::collections::HashSet;
use std::error::Error;
use std::ffi::CString;
use std::fs::File;
use std::i32;
use std::io::BufWriter;
use std::str::FromStr;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut input_path: Option<&str> = None;
    let mut primers_path: Option<&str> = None;
    let mut primerless_path: Option<&str> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--primers" => {
                primers_path = Some(&args[i + 1]);
                i += 2;
            }
            "--primerless" => {
                primerless_path = Some(&args[i + 1]);
                i += 2;
            }
            other => {
                input_path = Some(other);
                i += 1;
            }
        }
    }
    let input_path = input_path.expect("missing INPUT.bam argument");
    let primers_path = primers_path.expect("missing --primers argument");
    let primerless_path = primerless_path.expect("missing --primerless argument");

    let mut input = File::open(input_path).map(Reader::new)?;
    let header = input.read_header()?.parse()?;
    let reference_sequences = input.read_reference_sequences()?;

    let mut primerless_writer =
        build_writer(primerless_path, &header, &reference_sequences)?;
    let mut primer_writer = build_writer(primers_path, &header, &reference_sequences)?;

    let ra_tag: Tag = Tag::from_str("ra")?;
    let mut primary_records = HashSet::new();
    for result in input.records() {
        let record = result?;

        let data = record.data();
        match data.get(ra_tag) {
            Some(Ok(_)) => {
                let idx = i32::from(record.reference_sequence_id().unwrap());
                let chr = reference_sequences
                    .get_index(idx as usize)
                    .unwrap()
                    .0
                    .as_bytes();
                let pos: i32 = i32::from(record.position().unwrap());
                primary_records.insert((chr, pos, record.read_name().unwrap().to_owned()));
            }
            _ => continue,
        }
    }
    let mut input = File::open(input_path).map(Reader::new)?;
    input.read_header()?;
    input.read_reference_sequences()?;
    let ma_tag = Tag::from_str("ma")?;
    for result in input.records() {
        let record = result?;
        let data = record.data();
        if data.get(ra_tag).is_some()
            || data.get(ma_tag).is_some()
            || is_secondary_alignment(&record, &primary_records)?
        {
            primer_writer.write_record(&record)?;
        } else {
            primerless_writer.write_record(&record)?
        }
    }
    Ok(())
}

fn is_secondary_alignment(
    record: &Record,
    primary_records: &HashSet<(&[u8], i32, CString)>,
) -> Result<bool, Box<dyn Error>> {
    let data = record.data();
    match data.get(Tag::OtherAlignments) {
        Some(Ok(sa_entry)) => {
            let split_tag = sa_entry
                .value()
                .as_str()
                .unwrap()
                .split(',')
                .collect::<Vec<&str>>();
            let chrom = split_tag[0].as_bytes();
            let pos = split_tag[1].parse::<i32>()?;
            return Ok(primary_records.contains(&(
                chrom,
                pos,
                record.read_name().unwrap().to_owned(),
            )));
        }
        _ => Ok(false),
    }
}

fn build_writer(
    file_path: &str,
    header: &Header,
    reference_sequences: &IndexMap<String, ReferenceSequence>,
) -> Result<Writer<writer::Writer<BufWriter<File>>>, Box<dyn Error>> {
    let mut writer = std::fs::File::create(file_path)
        .map(BufWriter::new)
        .map(Writer::new)?;
    writer.write_header(header)?;
    writer.write_reference_sequences(reference_sequences)?;
    Ok(writer)
}
