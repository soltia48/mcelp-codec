use mcelp::{bitstream, lsp};

/// First frame of examples/male.mcelp, decoded through the LSF/LSP stage.
/// Expected values are the reference's own working memory (0x7bec.. LSF,
/// 0x6f34.. LSP).
#[test]
fn first_frame_lsf_and_lsp() {
    let frame = bitstream::parse_hex_line("27ef9c1000701f00070182000e03e000e000").unwrap();
    let words = bitstream::canonicalize(&frame);
    let p = bitstream::unpack(&words);
    assert_eq!(p.field[0], 0x27);
    assert_eq!(p.field[1], 0x0ef9);

    let mut st = lsp::LsfState::default();
    let idx = lsp::LsfIndices::unpack(p.field[0], p.field[1]);
    let lsf = lsp::decode(&mut st, idx);
    let want_lsf: [i16; 10] = [
        0x0ea3, 0x131a, 0x1ab9, 0x2262, 0x2c31, 0x34c4, 0x44f5, 0x4a36, 0, 0,
    ];
    for i in 0..8 {
        assert_eq!(lsf[i], want_lsf[i], "lsf[{i}]");
    }

    let l = lsp::lsf_to_lsp(&lsf);
    let want_lsp: [u16; 8] = [
        0x72d1, 0x69d9, 0x55e7, 0x3cf5, 0x1826, 0xf605, 0xb96d, 0xa8ef,
    ];
    for i in 0..8 {
        assert_eq!(l[i] as u16, want_lsp[i], "lsp[{i}]");
    }
}

/// The same frame carried on into LSP interpolation and LPC conversion.
#[test]
fn first_frame_lpc() {
    let frame = bitstream::parse_hex_line("27ef9c1000701f00070182000e03e000e000").unwrap();
    let p = bitstream::unpack(&bitstream::canonicalize(&frame));
    let mut st = lsp::LsfState::default();
    let lsf = lsp::decode(&mut st, lsp::LsfIndices::unpack(p.field[0], p.field[1]));
    let want_lsf: [u16; 10] = [
        0x0ea3, 0x131a, 0x1ab9, 0x2262, 0x2c31, 0x34c4, 0x44f5, 0x4a36, 0x514a, 0x54c2,
    ];
    assert_eq!(lsf.map(|v| v as u16), want_lsf);

    let cur = lsp::lsf_to_lsp(&lsf);
    let want_lsp: [u16; 10] = [
        0x72d1, 0x69d9, 0x55e7, 0x3cf5, 0x1826, 0xf605, 0xb96d, 0xa8ef, 0x9679, 0x8f3f,
    ];
    assert_eq!(cur.map(|v| v as u16), want_lsp);

    let mut prev = [0i16; 10];
    prev.copy_from_slice(&mcelp::tables::LSP_MEAN);
    let (mid, _, lpc_mid, lpc_cur) = lsp::interpolate_pair(&prev, &cur);
    let want_mid: [u16; 10] = [
        0x76d1, 0x6ac4, 0x54dd, 0x3910, 0x152e, 0xf1e7, 0xc220, 0xaa8d, 0x9564, 0x8a37,
    ];
    assert_eq!(mid.map(|v| v as u16), want_mid);
    let want_a1: [u16; 11] = [
        0x1000, 0xffa4, 0x00e2, 0x004c, 0x05ac, 0xfe4e, 0x036b, 0x00be, 0x01d4, 0xff0a, 0x0244,
    ];
    let want_a2: [u16; 11] = [
        0x1000, 0xff47, 0x01c9, 0x0095, 0x0b2b, 0xfc69, 0x06ce, 0x0181, 0x03e0, 0xfde5, 0x0488,
    ];
    assert_eq!(lpc_mid.map(|v| v as u16), want_a1);
    assert_eq!(lpc_cur.map(|v| v as u16), want_a2);
}
