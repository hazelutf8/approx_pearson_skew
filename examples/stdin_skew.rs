use approx_pearson_skew::pearson_skew_median;
use std::io::Read;

fn main() -> std::io::Result<()> {
    let mut lines = String::new();
    std::io::stdin().read_to_string(&mut lines)?;
    println!("Input from stdin: {:?}", lines);

    let skew = pearson_skew_median(lines.as_bytes()).expect("Zero bytes found");
    println!("Pearson (Median) Skew: {}", skew);
    Ok(())
}
