//! Shared SVG pieces: colour tokens for light and dark mode, label formatting, base64.

/// The `<svg>` opening tag and `<style>` block. Colours are role tokens with light and dark
/// values (OS setting, or a `data-theme` attribute on the `<svg>`); the reference palette of the
/// dataviz guide, validated for adjacent-pair and three-slot all-pair use.
pub fn open(width: f64, height: f64, title: &str) -> String {
    const LIGHT: &str = "color-scheme: light; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; \
        --muted: #898781; --grid: #e1e0d9; --axis: #c3c2b7; --s1: #2a78d6; --s2: #eb6834; \
        --s3: #1baf7a; --s4: #eda100; --s5: #e87ba4; --s6: #008300; --s7: #6250d6;";
    const DARK: &str = "color-scheme: dark; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; \
        --muted: #898781; --grid: #2c2c2a; --axis: #383835; --s1: #3987e5; --s2: #d95926; \
        --s3: #199e70; --s4: #c98500; --s5: #d55181; --s6: #008300; --s7: #9085e9;";
    format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" role="img" aria-label="{title}">
<title>{title}</title>
<style>
svg {{ {LIGHT} }}
@media (prefers-color-scheme: dark) {{ svg:where(:not([data-theme="light"])) {{ {DARK} }} }}
svg[data-theme="dark"] {{ {DARK} }}
text {{ font-family: system-ui, -apple-system, "Segoe UI", sans-serif; fill: var(--ink2); font-size: 12px; font-variant-numeric: tabular-nums; }}
.title {{ fill: var(--ink); font-size: 18px; font-weight: 600; font-variant-numeric: normal; }}
.subtitle {{ font-size: 13px; font-variant-numeric: normal; }}
.panel {{ fill: var(--ink); font-size: 14px; font-weight: 600; }}
.tick {{ fill: var(--muted); }}
.grid {{ stroke: var(--grid); stroke-width: 1; }}
.axis {{ stroke: var(--axis); stroke-width: 1; }}
.line {{ fill: none; stroke-width: 2; stroke-linejoin: round; stroke-linecap: round; }}
.ink {{ stroke: var(--ink); }} .ink2 {{ stroke: var(--ink2); }} .muted {{ stroke: var(--muted); }}
.s1 {{ stroke: var(--s1); }} .s2 {{ stroke: var(--s2); }} .s3 {{ stroke: var(--s3); }} .s4 {{ stroke: var(--s4); }}
.s5 {{ stroke: var(--s5); }} .s6 {{ stroke: var(--s6); }} .s7 {{ stroke: var(--s7); }}
</style>
<rect width="{width}" height="{height}" fill="var(--surface)"/>
"#
    )
}

/// `10⁻¹²`-style label for a power of ten.
pub fn pow10(exponent: i32) -> String {
    const SUP: [char; 10] = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
    let digits: String = exponent.unsigned_abs().to_string().chars().map(|c| SUP[c as usize - '0' as usize]).collect();
    format!("10{}{digits}", if exponent < 0 { "⁻" } else { "" })
}

/// Escape text for SVG content and attributes.
pub fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// Standard base64 (RFC 4648) for embedding a PNG as a data URI.
pub fn base64(bytes: &[u8]) -> String {
    const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let n = u32::from(chunk[0]) << 16 | u32::from(*chunk.get(1).unwrap_or(&0)) << 8 | u32::from(*chunk.get(2).unwrap_or(&0));
        for i in 0..4 {
            out.push(if i <= chunk.len() { A[(n >> (18 - 6 * i) & 63) as usize] as char } else { '=' });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn base64_and_labels() {
        assert_eq!(super::base64(b"Man"), "TWFu");
        assert_eq!(super::base64(b"Ma"), "TWE=");
        assert_eq!(super::base64(b"M"), "TQ==");
        assert_eq!(super::pow10(-12), "10⁻¹²");
        assert_eq!(super::pow10(3), "10³");
    }
}
