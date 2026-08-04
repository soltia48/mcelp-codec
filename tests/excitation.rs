use mcelp::{codebook, pitch};

#[test]
fn pitch_lag_decoding() {
    // Absolute lag, first subframe of male.mcelp frame 1.
    let a = pitch::decode_absolute(0xc1);
    assert_eq!((a.integer, a.frac), (0x54, -1));
    // Relative lag, second subframe of the same half-frame.
    let b = pitch::decode_relative(0x1f, a.integer);
    assert_eq!((b.integer, b.frac), (0x59, -1));
    // The integer-only branch at the top of the range.
    let c = pitch::decode_absolute(200);
    assert_eq!((c.integer, c.frac), (88, 0));
}

#[test]
fn innovation_index_7() {
    let inv = codebook::decode(7);
    assert_eq!(inv.class, 0);
    let want: [u16; 23] = [
        0x18e3, 0x1e7e, 0x0587, 0x187c, 0x039b, 0xf930, 0x00ea, 0xf960, 0xfe7d, 0xf9d2, 0xfca5,
        0xfaf8, 0xfc30, 0xfbb0, 0xfa95, 0xfe37, 0xf91f, 0xfe71, 0xf9ed, 0xfbea, 0xfc66, 0xfe1a,
        0xfe3c,
    ];
    let got: Vec<u16> = inv.code[..23].iter().map(|&v| v as u16).collect();
    assert_eq!(got, want.to_vec());
}

#[path = "data/mod.rs"]
mod data;

#[test]
fn adaptive_codebook_prediction() {
    // The excitation buffer is 0x7279..0x73b2; the subframe being built starts
    // at 0x7313, i.e. 154 words in.
    const AT: usize = 0x7313 - 0x7279;
    let mut exc: Vec<i16> = data::EXC_BEFORE.iter().map(|&v| v as i16).collect();
    let lag = pitch::Lag {
        integer: data::LAG.0,
        frac: data::LAG.1,
    };
    pitch::predict(&mut exc, AT, &lag);
    let got: Vec<u16> = exc[AT..AT + pitch::SUBFRAME]
        .iter()
        .map(|&v| v as u16)
        .collect();
    let want = &data::EXC_AFTER[AT..AT + pitch::SUBFRAME];
    assert_eq!(got, want.to_vec());
}
