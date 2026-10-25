use solar::surface::{Parameters, Surface};

pub const INITIAL: &str = include_str!("../../data/gong_20190108.csv");

pub fn map(csv: &str, parameters: Parameters) -> Result<Surface, &'static str> {
    let values = csv
        .lines()
        .filter(|line| !line.starts_with('#'))
        .flat_map(|line| line.split(','))
        .filter(|value| !value.trim().is_empty())
        .map(|value| {
            value
                .trim()
                .parse::<f64>()
                .map_err(|_| "invalid magnetogram")
        })
        .collect::<Result<Vec<_>, _>>()?;
    Surface::new(60, 120, parameters, values)
}
