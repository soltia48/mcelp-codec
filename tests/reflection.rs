//! Speech.mcelp frame 40, subframe 1.

use mcelp::lsp;

const LPC_TAIL: [u16; 10] = [
    0xffcf, 0x04ae, 0xf6cd, 0x028e, 0xf74d, 0x0694, 0xfedf, 0x05c7, 0xfe79, 0x0233,
];
const REFLECTION: [u16; 10] = [
    0x02f6, 0x01ad, 0xfaa9, 0x0104, 0xfa15, 0x04b0, 0x0064, 0x0547, 0xfe78, 0x0233,
];

#[test]
fn lpc_to_reflection() {
    let mut a = [0i16; 11];
    a[0] = 4096;
    for (i, &v) in LPC_TAIL.iter().enumerate() {
        a[i + 1] = v as i16;
    }
    let k = lsp::reflection_coefficients(&a);
    assert_eq!(
        k.iter().map(|&v| v as u16).collect::<Vec<_>>(),
        REFLECTION.to_vec()
    );
}
