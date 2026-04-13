//! Handler for the `vareffect init` subcommand.
//!
//! Scaffolds a `vareffect_build.toml` configuration file at a discoverable
//! location so that `setup`, `check`, and other subcommands can find it
//! without an explicit `--config` flag.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::Args;
use console::style;

/// The default `vareffect_build.toml` embedded at compile time.
const DEFAULT_CONFIG: &str = include_str!("../vareffect_build.toml");

/// Arguments for the `init` subcommand.
#[derive(Debug, Args)]
pub struct InitArgs {
    /// Path to write the config file.
    /// Overrides `VAREFFECT_BUILD_CONFIG` env and the XDG default.
    #[arg(long)]
    pub config: Option<PathBuf>,

    /// Overwrite an existing config without confirmation.
    #[arg(long)]
    pub force: bool,

    /// Set the data directory in the generated config.
    /// Overrides `output_dir` and `raw_dir` in the template.
    #[arg(long)]
    pub data_dir: Option<PathBuf>,
}

/// Resolve where to WRITE the config file.
///
/// Precedence: CLI flag > `VAREFFECT_BUILD_CONFIG` env > XDG default.
///
/// Unlike [`crate::config::find_config`] (which searches for an existing
/// file), this function returns a target path even if the file does not
/// exist yet.
fn resolve_config_write_path(cli_override: Option<&Path>) -> PathBuf {
    if let Some(p) = cli_override {
        return p.to_path_buf();
    }
    if let Ok(p) = std::env::var("VAREFFECT_BUILD_CONFIG") {
        return PathBuf::from(p);
    }
    dirs::config_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join("vareffect")
        .join("vareffect_build.toml")
}

/// Apply `--data-dir` overrides to the config template.
///
/// Replaces the default `output_dir` and `raw_dir` values with paths
/// derived from the given data directory.
fn apply_data_dir_override(template: &str, data_dir: &Path) -> String {
    let data_dir_str = data_dir.display().to_string().replace('\\', "/");
    let raw_dir_str = data_dir
        .join("raw")
        .display()
        .to_string()
        .replace('\\', "/");

    template
        .replace(
            "output_dir = \"data/vareffect\"",
            &format!("output_dir = \"{data_dir_str}\""),
        )
        .replace(
            "raw_dir = \"data/raw\"",
            &format!("raw_dir = \"{raw_dir_str}\""),
        )
}

/// Run the `vareffect init` subcommand.
///
/// # Errors
///
/// Returns an error if:
/// - The config file exists and the user declines to overwrite
/// - The parent directory cannot be created (permission denied)
/// - The config file cannot be written
pub fn run(args: &InitArgs) -> Result<()> {
    let config_path = resolve_config_write_path(args.config.as_deref());

    // Check for existing config.
    if config_path.exists() && !args.force {
        let term = console::Term::stdout();
        if !term.is_term() {
            bail!(
                "config already exists at {}. Use --force to overwrite.",
                config_path.display(),
            );
        }
        eprintln!(
            "Config already exists at {}. Overwrite? [y/N] ",
            config_path.display(),
        );
        let answer = term
            .read_line()
            .context("reading confirmation from terminal")?;
        let answer = answer.trim();
        if !answer.eq_ignore_ascii_case("y") {
            bail!("Aborted. Use --force to overwrite.");
        }
    }

    // Create parent directories.
    if let Some(parent) = config_path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("cannot create config directory: {}", parent.display(),))?;
    }

    // Generate config content.
    let content = match &args.data_dir {
        Some(dir) => apply_data_dir_override(DEFAULT_CONFIG, dir),
        None => DEFAULT_CONFIG.to_string(),
    };

    // Write config file.
    fs::write(&config_path, &content)
        .with_context(|| format!("writing config to {}", config_path.display()))?;

    eprintln!(
        "  {}  Config written to {}",
        style("OK").green(),
        config_path.display(),
    );
    eprintln!();
    eprintln!("Next steps:");
    eprintln!(
        "  1. Edit config if needed:  $EDITOR {}",
        config_path.display(),
    );
    eprintln!("  2. Download reference data: vareffect setup");
    eprintln!("  3. Validate setup:         vareffect check");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use tempfile::TempDir;

    #[test]
    fn default_template_is_identical_to_repo_file() {
        // The embedded template must match the repo file exactly.
        assert_eq!(DEFAULT_CONFIG, include_str!("../vareffect_build.toml"));
    }

    #[test]
    fn default_template_parses_as_valid_config() {
        let config: crate::config::VareffectConfig =
            toml::from_str(DEFAULT_CONFIG).expect("default template must parse");
        assert_eq!(config.vareffect.fasta_build, "GRCh38");
    }

    #[test]
    fn init_creates_config_file() {
        let tmp = TempDir::new().unwrap();
        let config_path = tmp.path().join("vareffect/vareffect_build.toml");

        let args = InitArgs {
            config: Some(config_path.clone()),
            force: false,
            data_dir: None,
        };
        run(&args).expect("init should succeed");

        let content = fs::read_to_string(&config_path).unwrap();
        assert_eq!(content, DEFAULT_CONFIG);
    }

    #[test]
    fn init_creates_parent_directories() {
        let tmp = TempDir::new().unwrap();
        let config_path = tmp.path().join("deep/nested/dir/vareffect_build.toml");

        let args = InitArgs {
            config: Some(config_path.clone()),
            force: false,
            data_dir: None,
        };
        run(&args).expect("init should create parent dirs");

        assert!(config_path.exists());
    }

    #[test]
    fn init_with_force_overwrites_existing() {
        let tmp = TempDir::new().unwrap();
        let config_path = tmp.path().join("vareffect_build.toml");
        fs::write(&config_path, "old content").unwrap();

        let args = InitArgs {
            config: Some(config_path.clone()),
            force: true,
            data_dir: None,
        };
        run(&args).expect("init --force should succeed");

        let content = fs::read_to_string(&config_path).unwrap();
        assert_eq!(content, DEFAULT_CONFIG);
    }

    #[test]
    fn data_dir_override_substitutes_paths() {
        let tmp = TempDir::new().unwrap();
        let config_path = tmp.path().join("vareffect_build.toml");

        let args = InitArgs {
            config: Some(config_path.clone()),
            force: false,
            data_dir: Some(PathBuf::from("/custom/data")),
        };
        run(&args).expect("init with --data-dir should succeed");

        let content = fs::read_to_string(&config_path).unwrap();
        assert!(content.contains("output_dir = \"/custom/data\""));
        assert!(content.contains("raw_dir = \"/custom/data/raw\""));
        assert!(!content.contains("output_dir = \"data/vareffect\""));
        assert!(!content.contains("raw_dir = \"data/raw\""));
    }

    #[test]
    fn data_dir_override_produces_valid_config() {
        let result = apply_data_dir_override(DEFAULT_CONFIG, Path::new("/opt/ve"));
        let _config: crate::config::VareffectConfig =
            toml::from_str(&result).expect("overridden template must still parse");
    }
}
